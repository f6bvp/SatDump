#include "recorder_fswitch.h"
#include "imgui/imgui.h"
#include "imgui/implot/implot.h"
#include "imgui/implot/implot_internal.h"
#include "logger.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <ctime>
#include <fstream>

namespace satdump
{
    static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

    void FrequencySwitcher::start(double center_hz, double sr, int fft_size, double fdt, const float *init_display)
    {
        stop();

        {
            std::lock_guard<std::mutex> lk(mtx);
            f0 = center_hz;
            samplerate = sr;
            nbins = fft_size;

            // Use an offset that is a whole number of FFT bins, so that
            // the folding is exact
            dk = std::max<int>(1, std::lround(offset_khz * 1e3 * nbins / samplerate));
            if (dk > nbins / 2 - 1)
                dk = nbins / 2 - 1;
            delta_eff = dk * samplerate / nbins;

            run_dwell = std::max(0.1f, dwell_s);
            run_settle = std::max(0.0f, settle_s);
            run_delay = std::max(0.0f, delay_s);
            calibrating = auto_delay;
            slots.clear();
            measured_delay = -1;
            delay_score_ratio = 0;
            slot_len = std::max(fdt, 0.05);

            frame_dt = fdt;
            events.clear();
            n_frames = 0;
            lag_s = 0;
            speed = 1;
            rate_hist.clear();
            speed_now = 1;
            sum_a.assign(nbins, 0);
            sum_b.assign(nbins, 0);
            n_a = n_b = 0;
            cycles = 0;
            plot_x.clear();
            plot_y.clear();
            disp.assign(nbins, -150.0f);
            // Keep showing the current spectrum until the first frames arrive
            if (init_display != nullptr)
                std::copy(init_display, init_display + nbins, disp.begin());
            disp_a.assign(nbins, 0.0f);
            disp_b.assign(nbins, 0.0f);
        }

        logger->info("Frequency switching: f0 %.0f Hz, offset +/- %.0f Hz (%d bins), dwell %.2f s, settle %.2f s", f0, delta_eff, dk, run_dwell, run_settle);

        running = true;
        thread = std::thread(&FrequencySwitcher::run, this);
    }

    void FrequencySwitcher::stop()
    {
        if (!thread.joinable())
            return;
        {
            std::lock_guard<std::mutex> lk(cv_mtx);
            running = false;
        }
        cv.notify_all();
        thread.join();
        logger->info("Frequency switching stopped after %llu cycles", (unsigned long long)cycles);
    }

    void FrequencySwitcher::reset()
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (auto_delay && running)
        {
            calibrating = true;
            slots.clear();
        }
        std::fill(sum_a.begin(), sum_a.end(), 0);
        std::fill(sum_b.begin(), sum_b.end(), 0);
        n_a = n_b = 0;
        cycles = 0;
        plot_x.clear();
        plot_y.clear();
        unf_x.clear();
        unf_12.clear();
        unf_21.clear();
    }

    void FrequencySwitcher::run()
    {
        int ph = 0;
        while (running)
        {
            set_sdr_frequency(f0 + (ph == 0 ? delta_eff : -delta_eff));

            double t = now_s();
            {
                std::lock_guard<std::mutex> lk(mtx);
                events.push_back({t, ph});
                while (events.size() > 64)
                    events.pop_front();
            }

            // While the delay is measured, B lasts twice as long as A: the
            // pattern is then not symmetric and only one delay matches it
            double dwell = run_dwell;
            {
                std::lock_guard<std::mutex> lk(mtx);
                if (calibrating && ph == 1)
                    dwell = 2 * run_dwell;
            }
            auto until = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(run_settle + dwell));
            {
                std::unique_lock<std::mutex> lk(cv_mtx);
                cv.wait_until(lk, until, [this] { return !running; });
            }

            if (running && ph == 1)
            {
                std::lock_guard<std::mutex> lk(mtx);
                cycles++;
            }
            ph ^= 1;
        }

        set_sdr_frequency(f0);
        std::lock_guard<std::mutex> lk(mtx);
        events.push_back({now_s(), -1});
    }

    void FrequencySwitcher::feed(const float *db, int size)
    {
        if (!running)
            return;

        double now = now_s();
        std::lock_guard<std::mutex> lk(mtx);
        if (size != nbins || frame_dt <= 0)
            return;

        // Date the frame by its arrival. The sample count only measures
        // losses (speed below 100 %)
        if (n_frames == 0)
            t_first = now - frame_dt;
        n_frames++;
        double fend = now;
        double fstart = now - frame_dt;
        lag_s = now - (t_first + n_frames * frame_dt);
        if (now > t_first)
            speed = (n_frames * frame_dt) / (now - t_first);

        // Same, over the last 5 seconds
        rate_hist.push_back({now, n_frames});
        while (rate_hist.size() > 2 && now - rate_hist.front().first > 5.0)
            rate_hist.pop_front();
        if (now - rate_hist.front().first >= 1.0)
            speed_now = (n_frames - rate_hist.front().second) * frame_dt / (now - rate_hist.front().first);

        // Main display, with the current delay
        int ph_disp = phase_of(fstart, fend, run_delay);
        if (ph_disp >= 0)
            update_display(db, ph_disp);

        // Keep short decimated slots, to measure then track the delay
        add_slot(db, fstart, fend);
        if (calibrating)
        {
            // Enough data: a few full (asymmetric) cycles plus one cycle for the delay
            double cycle = 2 * run_settle + 3 * run_dwell;
            if (fend - slots.front().t_start >= (CAL_CYCLES + 1) * cycle)
                estimate_delay();
            return;
        }
        else if (auto_delay)
        {
            double period = run_settle + run_dwell;
            if (fend - slots.front().t_start >= (2 * TRACK_CYCLES + 1) * period)
                track_delay();
        }

        int ph = phase_of(fstart, fend, run_delay);
        if (ph < 0)
            return;

        // dB -> linear power: 10^(x/10) = exp(x * ln(10) / 10)
        std::vector<double> &s = ph == 0 ? sum_a : sum_b;
        for (int k = 0; k < nbins; k++)
            s[k] += std::exp((double)db[k] * 0.23025850929940458);
        if (ph == 0)
            n_a++;
        else
            n_b++;
    }

    void FrequencySwitcher::update_display(const float *db, int ph)
    {
        // mtx must be held. Tuned at f0 + D (A), FFT bin k holds the sky
        // frequency of bin k + d around f0: move it back, and average each
        // position on its own.
        std::vector<float> &buf = ph == 0 ? disp_a : disp_b;
        int shift = ph == 0 ? dk : -dk;
        float r = 1.0f / std::max(1.0f, display_avg);
        for (int j = 0; j < nbins; j++)
        {
            int k = j - shift;
            if (k < 0 || k >= nbins)
                continue;
            float v = std::exp(db[k] * 0.23025850929940458f);
            buf[j] = buf[j] == 0 ? v : buf[j] * (1.0f - r) + v * r;
        }

        // Each bin is seen by both positions, at two places of the IF
        // filter. Weight each position by its own power: the one that sees
        // the bin inside its passband wins, the one on the filter edge fades
        // out, with no step between them.
        for (int j = 0; j < nbins; j++)
        {
            float a = disp_a[j], b = disp_b[j];
            if (a + b > 0)
                disp[j] = 10.0f * std::log10((a * a + b * b) / (a + b));
        }
        on_display(disp.data());
    }

    int FrequencySwitcher::phase_of(double fstart, double fend, double delay)
    {
        // mtx must be held. Time of the signal, as commanded
        fstart -= delay;
        fend -= delay;

        // Find the retune that applies to this frame
        int i = (int)events.size() - 1;
        while (i >= 0 && events[i].first > fstart)
            i--;
        if (i < 0 || events[i].second < 0)
            return -1;
        if (fstart < events[i].first + run_settle)
            return -1;
        double end_limit = (i + 1 < (int)events.size()) ? events[i + 1].first : events[i].first + run_settle + run_dwell;
        if (fend > end_limit)
            return -1;
        return events[i].second;
    }

    void FrequencySwitcher::add_slot(const float *db, double fstart, double fend)
    {
        // mtx must be held
        int ch = std::min(CAL_CHANNELS, nbins);
        if (slots.empty() || fstart >= slots.back().t_start + slot_len)
            slots.push_back({fstart, fend, std::vector<double>(ch, 0.0), 0});
        Slot &sl = slots.back();
        sl.t_end = fend;
        for (int k = 0; k < nbins; k++)
            sl.sum[(int64_t)k * ch / nbins] += std::exp((double)db[k] * 0.23025850929940458);
        sl.count++;
    }

    bool FrequencySwitcher::search_delay(double lo, double hi, double &best_d, double &contrast)
    {
        // mtx must be held. The right delay separates A from B best.
        int ch = slots.empty() ? 0 : (int)slots[0].sum.size();
        double best = -1, sum_scores = 0;
        int n_scores = 0;
        std::vector<std::pair<double, double>> scores; // (delay, score)
        for (double d = lo; d < hi; d += slot_len)
        {
            std::vector<double> a(ch, 0), b(ch, 0);
            int na = 0, nb = 0;
            for (auto &sl : slots)
            {
                int ph = phase_of(sl.t_start, sl.t_end, d);
                if (ph < 0)
                    continue;
                std::vector<double> &t = ph == 0 ? a : b;
                for (int k = 0; k < ch; k++)
                    t[k] += sl.sum[k] / sl.count;
                (ph == 0 ? na : nb)++;
            }
            if (na < 2 || nb < 2)
                continue;
            double score = 0;
            for (int k = 0; k < ch; k++)
            {
                double am = a[k] / na, bm = b[k] / nb;
                if (bm > 0)
                    score += std::pow((am - bm) / bm, 2);
            }
            score *= (double)na * nb / (na + nb); // same noise level for every delay
            sum_scores += score;
            n_scores++;
            scores.push_back({d, score});
            if (score > best)
                best = score;
        }
        if (n_scores == 0)
            return false;

        // Several delays can fit equally (the settle time gives some room):
        // take the middle of that range, to keep a margin on both sides
        int ib = 0;
        while (scores[ib].second < best)
            ib++;
        int l = ib, r = ib;
        while (l > 0 && scores[l - 1].second >= 0.95 * best)
            l--;
        while (r + 1 < (int)scores.size() && scores[r + 1].second >= 0.95 * best)
            r++;
        best_d = 0.5 * (scores[l].first + scores[r].first);
        contrast = best / (sum_scores / n_scores);
        return true;
    }

    void FrequencySwitcher::estimate_delay()
    {
        // mtx must be held. While measuring, B lasts twice as long as A, so
        // one asymmetric cycle holds only one matching delay
        double cycle = 2 * run_settle + 3 * run_dwell;
        double d = 0, contrast = 0;
        // The window of good delays can start slightly below zero: search
        // one full cycle starting a bit before zero, so that it is not cut
        // in two and wrapped to the end of the cycle
        double lo = -std::min(1.0, cycle / 4);
        bool ok = search_delay(lo, lo + cycle, d, contrast);
        d = std::max(0.0, d);

        calibrating = false;
        slots.clear();
        if (!ok)
        {
            logger->warn("Frequency switching: delay measurement failed, using %.2f s", run_delay);
            return;
        }
        measured_delay = d;
        delay_score_ratio = contrast;
        run_delay = d;
        std::fill(sum_a.begin(), sum_a.end(), 0);
        std::fill(sum_b.begin(), sum_b.end(), 0);
        n_a = n_b = 0;
        cycles = 0;
        logger->info("Frequency switching: measured delay %.2f s (contrast %.1f)", measured_delay, delay_score_ratio);
    }

    void FrequencySwitcher::track_delay()
    {
        // mtx must be held. Small corrections only: within half a period
        // the pattern holds a single matching delay
        double period = run_settle + run_dwell;
        double d = 0, contrast = 0;
        // Stay within half a period of the first measurement: drifting by a
        // whole period would swap A and B
        double centre = measured_delay >= 0 ? measured_delay : run_delay;
        double lo = std::max(-run_settle, centre - 0.45 * period), hi = centre + 0.45 * period;
        if (search_delay(lo, hi, d, contrast) && contrast > 1.5)
        {
            d = std::max(0.0, d);
            if (std::fabs(d - run_delay) >= slot_len)
                logger->debug("Frequency switching: delay %.2f -> %.2f s (contrast %.1f)", run_delay, d, contrast);
            run_delay = d;
            delay_score_ratio = contrast;
        }

        // Keep the last cycles only
        double keep_from = slots.back().t_end - 2 * TRACK_CYCLES * period;
        while (!slots.empty() && slots.front().t_end < keep_from)
            slots.erase(slots.begin());
    }

    void FrequencySwitcher::compute()
    {
        // mtx must be held
        plot_x.clear();
        plot_y.clear();
        unf_x.clear();
        unf_12.clear();
        unf_21.clear();
        if (n_a == 0 || n_b == 0 || nbins == 0)
            return;

        const double nan = std::numeric_limits<double>::quiet_NaN();

        std::vector<double> a(nbins), b(nbins), q(nbins);
        for (int k = 0; k < nbins; k++)
        {
            a[k] = sum_a[k] / n_a;
            b[k] = sum_b[k] / n_b;
        }

        // Hide the band edges: bins well below the median level
        std::vector<double> med = b;
        std::nth_element(med.begin(), med.begin() + nbins / 2, med.end());
        double threshold = med[nbins / 2] * std::pow(10.0, -edge_db / 10.0);
        std::vector<bool> valid(nbins);
        std::vector<double> q21(nbins); // o2 - o1, reference A
        for (int k = 0; k < nbins; k++)
        {
            valid[k] = a[k] > threshold && b[k] > threshold;
            q[k] = valid[k] ? (a[k] - b[k]) / b[k] : nan;      // o1 - o2, reference B
            q21[k] = valid[k] ? (b[k] - a[k]) / a[k] : nan; // o2 - o1, reference A
        }

        // FFT bin k is at (k / N - 0.5) * samplerate from the tuned frequency
        std::vector<double> y;
        for (int j = dk; j < nbins - dk; j++)
        {
            plot_x.push_back((f0 + ((double)j / nbins - 0.5) * samplerate) / 1e6);
            y.push_back((valid[j - dk] && valid[j + dk]) ? 0.5 * (q[j - dk] + q21[j + dk]) : nan);
        }
        smooth_scale(y, plot_y);

        // Unfolded, on the receiver axis (bin k shown at f0 + its offset)
        unf_x.clear();
        for (int k = 0; k < nbins; k++)
            unf_x.push_back((f0 + ((double)k / nbins - 0.5) * samplerate) / 1e6);
        smooth_scale(q, unf_12);
        smooth_scale(q21, unf_21);
    }

    void FrequencySwitcher::smooth_scale(const std::vector<double> &y, std::vector<double> &out)
    {
        // Boxcar smoothing ignoring hidden bins (prefix sums), then display units
        const double nan = std::numeric_limits<double>::quiet_NaN();
        int n = y.size();
        std::vector<double> ps(n + 1, 0);
        std::vector<int> pc(n + 1, 0);
        for (int i = 0; i < n; i++)
        {
            bool ok = !std::isnan(y[i]);
            ps[i + 1] = ps[i] + (ok ? y[i] : 0);
            pc[i + 1] = pc[i] + (ok ? 1 : 0);
        }
        int w = std::max(1, smoothing);
        out.resize(n);
        for (int i = 0; i < n; i++)
        {
            if (std::isnan(y[i]))
            {
                out[i] = nan;
                continue;
            }
            int lo = std::max(0, i - w / 2), hi = std::min(n - 1, i + (w - 1) / 2);
            double v = (ps[hi + 1] - ps[lo]) / (pc[hi + 1] - pc[lo]);
            if (use_db) // power ratio in dB; not defined below -100 %
                out[i] = 1 + v > 0 ? 10 * std::log10(1 + v) : nan;
            else
                out[i] = 1e3 * v; // in units of 1e-3 of the system power
        }
    }

    void FrequencySwitcher::draw_status()
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (nbins == 0)
            return;
        ImGui::Text("Offset : +/- %.1f kHz (%d bins)", delta_eff / 1e3, dk);
        ImGui::Text("Cycles : %llu (%.0f s on sky)", (unsigned long long)cycles, cycles * 2 * run_dwell);
        ImGui::Text("Frames : A %llu / B %llu", (unsigned long long)n_a, (unsigned long long)n_b);
        if (calibrating)
            ImGui::TextColored(ImVec4(1, 0.8, 0.2, 1), "Measuring the delay (%d cycles)...", CAL_CYCLES);
        else if (measured_delay >= 0)
            ImGui::Text("Delay : %.2f s tracked (contrast %.1f)", run_delay, delay_score_ratio);
        else
            ImGui::Text("Delay : %.2f s", run_delay);
        ImGui::Text("Stream : %.2f MS/s (last 5 s)", speed_now * samplerate / 1e6);
        ImGui::Text("Samples lost : %.1f %% now, %.1f %% since start", std::max(0.0, (1 - speed_now) * 100), std::max(0.0, (1 - speed) * 100));
        if (running && n_frames > 20 && speed_now < 0.9)
            ImGui::TextColored(ImVec4(1, 0.6, 0.2, 1), "More than 10 %% of the samples are lost");

        // The tracked delay must stay within half a period of the first
        // measurement, or A and B would be swapped: warn before
        double period = run_settle + run_dwell;
        if (!calibrating && measured_delay >= 0 && std::fabs(run_delay - measured_delay) > 0.35 * period)
            ImGui::TextColored(ImVec4(1, 0.3, 0.3, 1), "Delay moved by %.2f s since the start:\nincrease Dwell, or Reset to measure it again", run_delay - measured_delay);
    }

    void FrequencySwitcher::compute_once()
    {
        // mtx must be held
        int frame = ImGui::GetFrameCount();
        if (frame != computed_frame)
        {
            compute();
            computed_frame = frame;
        }
    }

    void FrequencySwitcher::plot_folded(ImVec2 size, bool main_area)
    {
        // mtx must be held
        if (ImPlot::BeginPlot("##fswitch", size))
        {
            ImPlot::SetupAxes("Frequency (MHz)", use_db ? "(Psys+dP)/Psys (dB)" : "dP/Psys x1e-3", ImPlotAxisFlags_AutoFit, auto_y ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None);
            float lo = use_db ? y_min_db : y_min, hi = use_db ? y_max_db : y_max;
            if (!auto_y)
                ImPlot::SetupAxisLimits(ImAxis_Y1, std::min(lo, hi), std::max(lo, hi) + 1e-6, ImPlotCond_Always);
            if (plot_x.size() > 0)
                ImPlot::PlotLine("folded", plot_x.data(), plot_y.data(), plot_x.size());
            // In the main area, deselecting "folded" in the legend stops switching
            if (main_area)
            {
                ImPlotItem *item = ImPlot::GetItem("folded");
                if (item != nullptr && !item->Show)
                {
                    item->Show = true;
                    stop_requested = true;
                }
            }
            ImPlot::EndPlot();
        }
    }

    void FrequencySwitcher::draw_plot(ImVec2 size)
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (nbins == 0)
            return;
        compute_once();
        plot_folded(size, true);
    }

    void FrequencySwitcher::draw_window()
    {
        bool folded = show_window && !running;
        if (!folded && !show_unfolded_window)
            return;
        std::lock_guard<std::mutex> lk(mtx);
        if (nbins == 0) // never started
            return;
        compute_once();

        if (folded)
        {
            ImGui::SetNextWindowSize(ImVec2(700, 400), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Frequency switching", &show_window))
                plot_folded(ImGui::GetContentRegionAvail(), false);
            ImGui::End();
        }

        if (show_unfolded_window)
        {
            ImGui::SetNextWindowSize(ImVec2(700, 400), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Frequency switching - unfolded", &show_unfolded_window))
            {
                if (ImPlot::BeginPlot("##fswitch_unf", ImGui::GetContentRegionAvail()))
                {
                    // Receiver axis: bin k shown at f0 + its offset from the tuned frequency
                    ImPlot::SetupAxes("Frequency (MHz)", use_db ? "power ratio (dB)" : "ratio x1e-3", ImPlotAxisFlags_AutoFit, auto_y ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None);
                    float lo = use_db ? y_min_db : y_min, hi = use_db ? y_max_db : y_max;
                    if (!auto_y)
                        ImPlot::SetupAxisLimits(ImAxis_Y1, std::min(lo, hi), std::max(lo, hi) + 1e-6, ImPlotCond_Always);
                    if (unf_x.size() > 0)
                    {
                        ImPlot::PlotLine(use_db ? "o1-o2 : A/B" : "o1-o2 : (A-B)/B", unf_x.data(), unf_12.data(), unf_x.size());
                        ImPlot::PlotLine(use_db ? "o2-o1 : B/A" : "o2-o1 : (B-A)/A", unf_x.data(), unf_21.data(), unf_x.size());
                    }
                    ImPlot::EndPlot();
                }
            }
            ImGui::End();
        }
    }

    std::string FrequencySwitcher::save_csv(std::string directory)
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (n_a == 0 || n_b == 0 || nbins == 0)
            return "";

        time_t t = time(nullptr);
        struct tm tmv;
#ifdef _WIN32
        gmtime_s(&tmv, &t);
#else
        gmtime_r(&t, &tmv);
#endif
        char ts[32];
        strftime(ts, sizeof(ts), "%Y-%m-%d_%H-%M-%S", &tmv);
        std::string path = directory + "/fswitch_" + std::to_string((uint64_t)f0) + "Hz_" + ts + "Z.csv";

        std::ofstream out(path);
        if (!out)
            return "";

        out << "# SatDump frequency switching (experimental)\n";
        out << "# f0_hz=" << (uint64_t)f0 << " samplerate=" << samplerate << " fft_size=" << nbins << "\n";
        out << "# offset_hz=+/-" << delta_eff << " offset_bins=" << dk << " dwell_s=" << run_dwell << " settle_s=" << run_settle << " delay_s=" << run_delay << "\n";
        out << "# cycles=" << cycles << " frames_a=" << n_a << " frames_b=" << n_b << "\n";
        out << "# A = mean linear power tuned at f0+offset, B = tuned at f0-offset, indexed by FFT bin\n";
        out << "# q = (A-B)/B ; folded(j) = ((A-B)/B (j-d) + (B-A)/A (j+d))/2, sky frequency rf_hz, empty near the edges\n";
        out << "bin,if_offset_hz,rf_hz,a_lin,b_lin,q,folded\n";
        out.precision(10);
        for (int j = 0; j < nbins; j++)
        {
            double off = ((double)j / nbins - 0.5) * samplerate;
            double a = sum_a[j] / n_a, b = sum_b[j] / n_b;
            double q = b > 0 ? (a - b) / b : 0;
            out << j << "," << off << "," << (uint64_t)std::llround(f0 + off) << "," << a << "," << b << "," << q << ",";
            if (j >= dk && j < nbins - dk)
            {
                double aa = sum_a[j - dk] / n_a, ba = sum_b[j - dk] / n_b;
                double ab = sum_a[j + dk] / n_a, bb = sum_b[j + dk] / n_b;
                double q12 = ba > 0 ? (aa - ba) / ba : 0; // o1 - o2, reference B
                double q21 = ab > 0 ? (bb - ab) / ab : 0; // o2 - o1, reference A
                out << 0.5 * (q12 + q21);
            }
            out << "\n";
        }
        return path;
    }
} // namespace satdump
