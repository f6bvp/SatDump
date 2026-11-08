#include "recorder_fswitch.h"
#include "core/config.h"
#include "imgui/imgui.h"
#include "imgui/implot/implot.h"
#include "imgui/implot/implot_internal.h"
#include "logger.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <ctime>
#include <fstream>
#include <filesystem>

namespace satdump
{
    static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

    int FrequencySwitcher::dk_for(double off_khz, double sr, int fft_size)
    {
        // Use an offset that is a whole number of FFT bins, so that
        // the folding is exact
        int d = std::max<int>(1, std::lround(off_khz * 1e3 * fft_size / sr));
        if (d > fft_size / 2 - 1)
            d = fft_size / 2 - 1;
        return d;
    }

    bool FrequencySwitcher::compatible(double center_hz, double sr, int fft_size, int new_dk)
    {
        // mtx must be held
        return n_a > 0 && n_b > 0 && center_hz == f0 && sr == samplerate && fft_size == nbins && new_dk == dk && (int)sum_a.size() == nbins && (int)sum_b.size() == nbins;
    }

    bool FrequencySwitcher::can_resume(double center_hz, double sr, int fft_size, double off_khz)
    {
        if (running || sr <= 0 || fft_size <= 0)
            return false;
        std::lock_guard<std::mutex> lk(mtx);
        return compatible(center_hz, sr, fft_size, dk_for(off_khz, sr, fft_size));
    }

    bool FrequencySwitcher::init_run(double center_hz, double sr, int fft_size, const float *init_display, bool resume)
    {
        // mtx must be held. Common to local and server switching
        int new_dk = dk_for(offset_khz, sr, fft_size);
        bool keep = resume && compatible(center_hz, sr, fft_size, new_dk);
        if (resume && !keep && (n_a > 0 || n_b > 0))
            logger->warn("Frequency switching: frequency, samplerate, FFT size or offset changed, starting from zero");

        f0 = center_hz;
        samplerate = sr;
        nbins = fft_size;
        dk = new_dk;
        delta_eff = dk * samplerate / nbins;

        // Same limits as the panel
        run_dwell = std::clamp(dwell_s, DWELL_MIN_S, DWELL_MAX_S);
        run_settle = std::clamp(settle_s, 0.0f, SETTLE_MAX_S);

        events.clear();
        n_frames = 0;
        lag_s = 0;
        speed = 1;
        rate_hist.clear();
        speed_now = 1;
        if (!keep)
        {
            sum_a.assign(nbins, 0);
            sum_b.assign(nbins, 0);
            n_a = n_b = 0;
            cycles = 0;
            plot_x.clear();
            plot_y.clear();
            unf_x.clear();
            unf_12.clear();
            unf_21.clear();
        }
        n_a_run0 = n_a;
        n_b_run0 = n_b;
        cycles_run0 = cycles;
        resumed = keep;
        // A cycle cut by the stop is not written: start a new one
        std::fill(cyc_a.begin(), cyc_a.end(), 0);
        std::fill(cyc_b.begin(), cyc_b.end(), 0);
        cyc_na = cyc_nb = 0;
        cyc_last_ph = -1;
        disp.assign(nbins, -150.0f);
        // Keep showing the current spectrum until the first frames arrive
        if (init_display != nullptr)
            std::copy(init_display, init_display + nbins, disp.begin());
        disp_a.assign(nbins, 0.0f);
        disp_b.assign(nbins, 0.0f);

        close_cycles();
        if (save_cycles)
            open_cycles(keep && cyc_path.size() > 0);
        else if (!keep)
            cyc_path.clear();

        remote_status_seen = false;
        remote_iq_file.clear();
        remote_iq_error.clear();
        remote_iq_bytes = 0;
        remote_dwells = 0;
        return keep;
    }

    void FrequencySwitcher::start(double center_hz, double sr, int fft_size, double fdt, const float *init_display, bool resume)
    {
        stop();

        bool kept;
        {
            std::lock_guard<std::mutex> lk(mtx);
            remote = false;
            kept = init_run(center_hz, sr, fft_size, init_display, resume);
            frame_dt = fdt;
        }

        logger->info("Frequency switching%s: f0 %.0f Hz, offset +/- %.0f Hz (%d bins), dwell %.2f s, settle %.2f s", kept ? " resumed" : "", f0, delta_eff, dk, run_dwell, run_settle);

        running = true;
        thread = std::thread(&FrequencySwitcher::run, this);
    }

    bool FrequencySwitcher::start_remote(double center_hz, double sr, int fft_size, const float *init_display, bool resume)
    {
        stop();

        nlohmann::json st;
        uint64_t gen;
        bool kept;
        {
            std::lock_guard<std::mutex> lk(mtx);
            remote = true; // before init_run: written in the cycles file header
            kept = init_run(center_hz, sr, fft_size, init_display, resume);
            gen = ++remote_gen;
            // A = sky f0 + D, B = sky f0 - D, through the xconverter
            st["freq_a"] = sdr_freq_of(f0 + delta_eff);
            st["freq_b"] = sdr_freq_of(f0 - delta_eff);
            st["freq_back"] = sdr_freq_of(f0);
            st["sky_f0"] = f0;
            st["offset_hz"] = delta_eff;
            st["fft_size"] = nbins;
            st["dwell_s"] = run_dwell;
            st["settle_s"] = run_settle;
            st["save_iq"] = save_iq_on_server;
            st["iq_bits"] = iq_bits;
            remote = true;
            remote_t_start = now_s();
        }

        auto cb = std::make_shared<dsp::remote_fswitch::Callbacks>();
        std::weak_ptr<int> alive = alive_token;
        cb->on_status = [this, alive, gen](nlohmann::json j)
        {
            if (auto a = alive.lock())
            {
                if (gen == remote_gen)
                    on_remote_status(j);
            }
        };
        cb->on_dwell = [this, alive, gen](int phase, int frames, uint64_t index, const float *mean, int n)
        {
            if (auto a = alive.lock())
            {
                if (gen == remote_gen)
                    on_remote_dwell(phase, frames, index, mean, n);
            }
        };

        running = true;
        if (!remote_start(st, cb))
        {
            running = false;
            std::lock_guard<std::mutex> lk(mtx);
            remote = false;
            return false;
        }

        logger->info("Frequency switching on the server%s: f0 %.0f Hz, offset +/- %.0f Hz (%d bins), dwell %.2f s, settle %.2f s%s", kept ? " resumed" : "", f0, delta_eff, dk, run_dwell, run_settle,
                     save_iq_on_server ? ", IQ saved on the server" : "");
        return true;
    }

    void FrequencySwitcher::on_remote_status(nlohmann::json st)
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (!remote)
            return;
        remote_status_seen = true;
        if (st.contains("iq_file"))
            remote_iq_file = st["iq_file"].get<std::string>();
        if (st.contains("iq_bytes"))
            remote_iq_bytes = st["iq_bytes"].get<uint64_t>();
        if (st.contains("iq_error"))
            remote_iq_error = st["iq_error"].get<std::string>();
        if (st.contains("dwells"))
            remote_dwells = st["dwells"].get<uint64_t>();

        bool srv_running = st.value("running", false);
        std::string err = st.value("error", std::string());
        if (running && !srv_running)
        {
            // Refused at start (or stopped by the server): switch locally
            // if it never started, else just stop
            if (n_a == n_a_run0 && n_b == n_b_run0 && remote_dwells == 0)
            {
                fallback_why = err.size() ? err : "refused by the server";
                fallback_requested = true;
            }
            else
            {
                logger->warn("Frequency switching stopped by the server%s", err.size() ? (": " + err).c_str() : "");
                stop_requested = true;
            }
        }
        else if (srv_running && st.contains("samplerate") && st["samplerate"].get<double>() != samplerate)
        {
            fallback_why = "the server samplerate differs from the recorder";
            fallback_requested = true;
        }
    }

    void FrequencySwitcher::on_remote_dwell(int phase, int frames, uint64_t index, const float *mean, int n)
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (!remote || !running || n != nbins || frames <= 0 || (phase != 0 && phase != 1))
            return;
        remote_status_seen = true;
        remote_dwells = index + 1;

        // The dwell mean holds `frames` FFT frames: sum them as such
        std::vector<double> &s = phase == 0 ? sum_a : sum_b;
        for (int k = 0; k < nbins; k++)
            s[k] += (double)mean[k] * frames;
        if (phase == 0)
            n_a += frames;
        else
        {
            n_b += frames;
            cycles++;
        }

        cycle_add(phase, mean, nullptr, frames);

        // Each dwell is already an average of many frames: show the latest
        update_display_lin(mean, phase, 1.0f);
    }

    bool FrequencySwitcher::take_fallback(std::string &why)
    {
        {
            std::lock_guard<std::mutex> lk(mtx);
            // An older server does not know the request: no answer at all
            if (remote && running && !remote_status_seen && now_s() - remote_t_start > 3.0)
            {
                fallback_why = "no answer from the server (older SatDump server?)";
                fallback_requested = true;
            }
        }
        if (!fallback_requested.exchange(false))
            return false;
        std::lock_guard<std::mutex> lk(mtx);
        why = fallback_why;
        return true;
    }

    void FrequencySwitcher::stop()
    {
        if (remote && running)
        {
            running = false;
            remote_stop(); // the server tunes back to f0
            logger->info("Frequency switching on the server stopped after %llu cycles", (unsigned long long)cycles);
            std::lock_guard<std::mutex> lk(mtx);
            close_cycles();
            return;
        }
        if (!thread.joinable())
            return;
        {
            std::lock_guard<std::mutex> lk(cv_mtx);
            running = false;
        }
        cv.notify_all();
        thread.join();
        logger->info("Frequency switching stopped after %llu cycles", (unsigned long long)cycles);
        std::lock_guard<std::mutex> lk(mtx);
        close_cycles();
    }

    void FrequencySwitcher::reset()
    {
        std::lock_guard<std::mutex> lk(mtx);
        std::fill(sum_a.begin(), sum_a.end(), 0);
        std::fill(sum_b.begin(), sum_b.end(), 0);
        n_a = n_b = 0;
        cycles = 0;
        n_a_run0 = n_b_run0 = cycles_run0 = 0;
        resumed = false;
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

            auto until = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(run_settle + run_dwell));
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
        if (!running || remote)
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

        // Local SDR: its samples arrive within Settle after a retune
        int ph = phase_of(fstart, fend, 0.0);
        if (ph < 0)
            return;
        update_display(db, ph);

        // dB -> linear power: 10^(x/10) = exp(x * ln(10) / 10)
        std::vector<double> &s = ph == 0 ? sum_a : sum_b;
        lin_tmp.resize(nbins);
        for (int k = 0; k < nbins; k++)
        {
            lin_tmp[k] = std::exp((double)db[k] * 0.23025850929940458);
            s[k] += lin_tmp[k];
        }
        cycle_add(ph, nullptr, lin_tmp.data(), 1.0);
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

    void FrequencySwitcher::update_display_lin(const float *lin, int ph, float r)
    {
        // mtx must be held. As update_display(), from a linear spectrum
        std::vector<float> &buf = ph == 0 ? disp_a : disp_b;
        int shift = ph == 0 ? dk : -dk;
        for (int j = 0; j < nbins; j++)
        {
            int k = j - shift;
            if (k < 0 || k >= nbins)
                continue;
            float v = lin[k];
            buf[j] = buf[j] == 0 ? v : buf[j] * (1.0f - r) + v * r;
        }
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

        std::vector<double> a, b, q(nbins);
        mean_spectra(a, b);

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
        folded_raw = y;
        bl_curve.clear();
        bl_points = 0;
        bl_rms = 0;
        if (bl_remove)
            remove_baseline(plot_x, y);
        smooth_scale(y, plot_y);

        // Unfolded, on the receiver axis (bin k shown at f0 + its offset)
        unf_x.clear();
        for (int k = 0; k < nbins; k++)
            unf_x.push_back((f0 + ((double)k / nbins - 0.5) * samplerate) / 1e6);
        smooth_scale(q, unf_12);
        smooth_scale(q21, unf_21);
    }

    void FrequencySwitcher::rfi_detect(const std::vector<double> &m, std::vector<uint8_t> &flag, int &count)
    {
        // mtx must be held
        int n = m.size();
        flag.assign(n, 0);
        count = 0;
        int block = std::max(5, (int)std::lround(rfi_width_khz * 1e3 * n / samplerate));
        int nblocks = n / block;
        if (nblocks < 3)
            return;

        // Local level: median of each block, linearly interpolated between
        // block centres. A spur narrower than half a block does not move it
        std::vector<double> centre(nblocks), level(nblocks), tmp;
        for (int i = 0; i < nblocks; i++)
        {
            int lo = i * block, hi = (i == nblocks - 1) ? n : lo + block;
            tmp.assign(m.begin() + lo, m.begin() + hi);
            std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
            level[i] = tmp[tmp.size() / 2];
            centre[i] = 0.5 * (lo + hi - 1);
        }
        std::vector<double> r(n);
        int i = 0;
        for (int k = 0; k < n; k++)
        {
            while (i < nblocks - 2 && k > centre[i + 1])
                i++;
            double t = std::clamp((k - centre[i]) / (centre[i + 1] - centre[i]), 0.0, 1.0);
            double base = level[i] + t * (level[i + 1] - level[i]);
            r[k] = base > 0 ? m[k] / base - 1 : 0;
        }

        // Noise: robust spread of the residual over the whole band
        tmp = r;
        std::nth_element(tmp.begin(), tmp.begin() + n / 2, tmp.end());
        double med = tmp[n / 2];
        for (int k = 0; k < n; k++)
            tmp[k] = std::fabs(r[k] - med);
        std::nth_element(tmp.begin(), tmp.begin() + n / 2, tmp.end());
        double sigma = 1.4826 * tmp[n / 2];
        if (!(sigma > 0))
            return;

        // Spurs stand above the local level only
        double lim = med + rfi_sigma * sigma;
        for (int k = 0; k < n; k++)
        {
            if (r[k] > lim)
            {
                for (int j = std::max(0, k - rfi_grow); j <= std::min(n - 1, k + rfi_grow); j++)
                    flag[j] = 1;
            }
        }
        for (int k = 0; k < n; k++)
            count += flag[k];
    }

    void FrequencySwitcher::rfi_fill(std::vector<double> &m, const std::vector<uint8_t> &flag)
    {
        // Each run of flagged bins: straight line between its clean neighbours
        int n = m.size();
        if ((int)flag.size() != n)
            return;
        int k = 0;
        while (k < n)
        {
            if (!flag[k])
            {
                k++;
                continue;
            }
            int start = k;
            while (k < n && flag[k])
                k++;
            int left = start - 1, right = k; // clean neighbours, may be outside
            if (left < 0 && right >= n)
                return;
            double vl = left >= 0 ? m[left] : m[right];
            double vr = right < n ? m[right] : m[left];
            for (int j = start; j < k; j++)
                m[j] = vl + (vr - vl) * (j - left) / (double)(right - left);
        }
    }

    void FrequencySwitcher::mean_spectra(std::vector<double> &a, std::vector<double> &b)
    {
        // mtx must be held, n_a and n_b > 0
        a.resize(nbins);
        b.resize(nbins);
        for (int k = 0; k < nbins; k++)
        {
            a[k] = sum_a[k] / n_a;
            b[k] = sum_b[k] / n_b;
        }
        if (!rfi_remove)
        {
            rfi_a.clear();
            rfi_b.clear();
            rfi_count_a = rfi_count_b = 0;
            rfi_na = rfi_nb = 0;
            return;
        }
        // Detect again only when the averages or the settings changed
        if (n_a != rfi_na || n_b != rfi_nb || rfi_sigma != rfi_last_sigma || rfi_width_khz != rfi_last_width || (int)rfi_a.size() != nbins)
        {
            rfi_detect(a, rfi_a, rfi_count_a);
            rfi_detect(b, rfi_b, rfi_count_b);
            rfi_na = n_a;
            rfi_nb = n_b;
            rfi_last_sigma = rfi_sigma;
            rfi_last_width = rfi_width_khz;
        }
        rfi_fill(a, rfi_a);
        rfi_fill(b, rfi_b);
    }

    void FrequencySwitcher::remove_baseline(const std::vector<double> &x, std::vector<double> &y)
    {
        // Least squares with Legendre polynomials on x mapped to [-1, 1]
        // (well conditioned), points of the line window left out, then two
        // passes rejecting points beyond 3 sigma (residual spurs, noise peaks)
        int n = y.size(), m = std::clamp(bl_order, 0, 9) + 1;
        if (n < 2 * m)
            return;
        double x0 = x.front(), x1 = x.back();
        double flo = mhz_of_velocity(bl_vmax), fhi = mhz_of_velocity(bl_vmin); // higher v = lower f
        if (flo > fhi)
            std::swap(flo, fhi);
        std::vector<double> P(m);
        auto basis = [&](double xv)
        {
            double t = 2 * (xv - x0) / (x1 - x0) - 1;
            P[0] = 1;
            if (m > 1)
                P[1] = t;
            for (int k = 2; k < m; k++)
                P[k] = ((2 * k - 1) * t * P[k - 1] - (k - 1) * P[k - 2]) / k;
        };
        std::vector<char> use(n);
        for (int i = 0; i < n; i++)
            use[i] = !std::isnan(y[i]) && !(x[i] >= flo && x[i] <= fhi);
        std::vector<double> coef(m, 0.0);
        bool ok = false;
        for (int pass = 0; pass < 3; pass++)
        {
            std::vector<double> A(m * m, 0.0), r(m, 0.0);
            int cnt = 0;
            for (int i = 0; i < n; i++)
            {
                if (!use[i])
                    continue;
                basis(x[i]);
                for (int a = 0; a < m; a++)
                {
                    r[a] += P[a] * y[i];
                    for (int b = 0; b < m; b++)
                        A[a * m + b] += P[a] * P[b];
                }
                cnt++;
            }
            if (cnt < 2 * m)
                break;
            // Gaussian elimination with partial pivoting
            ok = true;
            for (int c = 0; c < m && ok; c++)
            {
                int pv = c;
                for (int k = c + 1; k < m; k++)
                    if (std::fabs(A[k * m + c]) > std::fabs(A[pv * m + c]))
                        pv = k;
                if (std::fabs(A[pv * m + c]) < 1e-300)
                {
                    ok = false;
                    break;
                }
                if (pv != c)
                {
                    for (int k = 0; k < m; k++)
                        std::swap(A[c * m + k], A[pv * m + k]);
                    std::swap(r[c], r[pv]);
                }
                for (int k = c + 1; k < m; k++)
                {
                    double f = A[k * m + c] / A[c * m + c];
                    for (int l = c; l < m; l++)
                        A[k * m + l] -= f * A[c * m + l];
                    r[k] -= f * r[c];
                }
            }
            if (!ok)
                break;
            for (int c = m - 1; c >= 0; c--)
            {
                double v = r[c];
                for (int k = c + 1; k < m; k++)
                    v -= A[c * m + k] * coef[k];
                coef[c] = v / A[c * m + c];
            }
            // Residual spread, then reject outliers for the next pass
            double ss = 0;
            for (int i = 0; i < n; i++)
            {
                if (!use[i])
                    continue;
                basis(x[i]);
                double fit = 0;
                for (int k = 0; k < m; k++)
                    fit += coef[k] * P[k];
                ss += (y[i] - fit) * (y[i] - fit);
            }
            bl_points = cnt;
            bl_rms = std::sqrt(ss / cnt);
            if (pass == 2)
                break;
            for (int i = 0; i < n; i++)
            {
                if (!use[i])
                    continue;
                basis(x[i]);
                double fit = 0;
                for (int k = 0; k < m; k++)
                    fit += coef[k] * P[k];
                if (std::fabs(y[i] - fit) > 3 * bl_rms)
                    use[i] = 0;
            }
        }
        if (!ok)
            return;
        bl_curve.assign(n, 0.0);
        for (int i = 0; i < n; i++)
        {
            basis(x[i]);
            double fit = 0;
            for (int k = 0; k < m; k++)
                fit += coef[k] * P[k];
            bl_curve[i] = fit;
            if (!std::isnan(y[i]))
                y[i] -= fit;
        }
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
        if (bl_remove && bl_points > 0)
            ImGui::Text("Baseline : order %d, rms %.3g x1e-3 outside the line (%d points)", bl_order, bl_rms * 1e3, bl_points);
        if (rfi_remove && (int)rfi_a.size() == nbins)
            ImGui::Text("Interference removed : A %d / B %d bins", rfi_count_a, rfi_count_b);
        if (resumed)
            ImGui::TextDisabled("Resumed: %llu cycles kept from before", (unsigned long long)cycles_run0);
        if (cyc_path.size())
            ImGui::TextWrapped("Cycles saved : %llu in %s", (unsigned long long)cyc_count, cyc_path.c_str());
        if (cyc_error.size())
            ImGui::TextColored(ImVec4(1, 0.3, 0.3, 1), "Cycles : %s", cyc_error.c_str());
        if (remote)
        {
            ImGui::Text("Switching done by the server");
            if (running && !remote_status_seen)
                ImGui::TextColored(ImVec4(1, 0.8, 0.2, 1), "Waiting for the server...");
            ImGui::Text("Dwells received : %llu", (unsigned long long)remote_dwells);
            if (remote_iq_file.size())
                ImGui::TextWrapped("IQ on the server : %s (%.2f GB)", remote_iq_file.c_str(), remote_iq_bytes / 1e9);
            if (remote_iq_error.size())
                ImGui::TextColored(ImVec4(1, 0.3, 0.3, 1), "IQ : %s", remote_iq_error.c_str());
            return;
        }
        ImGui::Text("Stream : %.2f MS/s (last 5 s)", speed_now * samplerate / 1e6);
        ImGui::Text("Samples lost : %.1f %% now, %.1f %% since start", std::max(0.0, (1 - speed_now) * 100), std::max(0.0, (1 - speed) * 100));
        if (running && n_frames > 20 && speed_now < 0.9)
            ImGui::TextColored(ImVec4(1, 0.6, 0.2, 1), "More than 10 %% of the samples are lost");
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

    void FrequencySwitcher::galactic_to_j2000(double l, double b, double &ra_deg, double &dec_deg)
    {
        // North galactic pole and longitude of the north celestial pole (J2000)
        const double D = M_PI / 180, ag = 192.85948 * D, dg = 27.12825 * D, lncp = 122.93192 * D;
        l *= D;
        b *= D;
        double sd = std::sin(b) * std::sin(dg) + std::cos(b) * std::cos(dg) * std::cos(lncp - l);
        double y = std::cos(b) * std::sin(lncp - l);
        double x = std::sin(b) * std::cos(dg) - std::cos(b) * std::sin(dg) * std::cos(lncp - l);
        ra_deg = std::fmod((ag + std::atan2(y, x)) / D + 360.0, 360.0);
        dec_deg = std::asin(sd) / D;
    }

    void FrequencySwitcher::lsr_correction(double jd, double lat, double lon, double alt_m, double ra, double dec, double &rot, double &orb, double &sun)
    {
        // Velocities projected on the source direction (km/s), to add to the
        // observed radial velocity. Low precision formulas (Meeus), J2000
        const double D = M_PI / 180;
        double s[3] = {std::cos(dec * D) * std::cos(ra * D), std::cos(dec * D) * std::sin(ra * D), std::sin(dec * D)};
        double T = (jd - 2451545.0) / 36525;

        // Earth orbit: from the geometric longitude of the Sun (Meeus ch. 25)
        double L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
        double M = (357.52911 + 35999.05029 * T - 0.0001537 * T * T) * D;
        double C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * std::sin(M) + (0.019993 - 0.000101 * T) * std::sin(2 * M) + 0.000289 * std::sin(3 * M);
        double e = 0.016708634 - 0.000042037 * T;
        double prec = 1.396971 * T + 0.0003086 * T * T; // equinox of date -> J2000
        double lsun = (L0 + C - prec) * D, peri = (102.93735 + 1.71946 * T - prec) * D;
        const double K = 29.7885, eps = 23.4392911 * D;
        double vx = K * (std::sin(lsun) - e * std::sin(peri));
        double vy = K * (-std::cos(lsun) + e * std::cos(peri));
        orb = vx * s[0] + vy * std::cos(eps) * s[1] + vy * std::sin(eps) * s[2];

        // Earth rotation, from the local sidereal time
        double d = jd - 2451545.0;
        double gmst = std::fmod(280.46061837 + 360.98564736629 * d + 0.000387933 * T * T, 360.0);
        double lst = (gmst + lon) * D;
        double vr = 7.2921159e-5 * (6378.137 + alt_m / 1000) * std::cos(lat * D);
        rot = vr * (-std::sin(lst) * s[0] + std::cos(lst) * s[1]);

        // Sun relative to the LSR
        const double ar = 270.959542 * D, dr = 30.004667 * D;
        sun = 20.0 * (std::cos(dr) * std::cos(ar) * s[0] + std::cos(dr) * std::sin(ar) * s[1] + std::sin(dr) * s[2]);
    }

    void FrequencySwitcher::update_lsr()
    {
        if (!lsr_qth_ok)
        {
            try
            {
                qth_lat = satdump_cfg.getValueFromSatDumpGeneral<double>("qth_lat");
                qth_lon = satdump_cfg.getValueFromSatDumpGeneral<double>("qth_lon");
                qth_alt = satdump_cfg.getValueFromSatDumpGeneral<double>("qth_alt");
                lsr_qth_ok = true;
            }
            catch (std::exception &)
            {
            }
        }
        if (!lsr)
        {
            lsr_corr = lsr_rot = lsr_orb = lsr_sun = 0;
            return;
        }
        double ra, dec;
        if (lsr_frame == 1)
            galactic_to_j2000(src_l, src_b, ra, dec);
        else
        {
            ra = src_ra_h * 15.0;
            dec = src_dec;
        }
        double jd = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count() / 86400.0 + 2440587.5;
        lsr_correction(jd, qth_lat, qth_lon, qth_alt, ra, dec, lsr_rot, lsr_orb, lsr_sun);
        lsr_corr = lsr_rot + lsr_orb + lsr_sun;
    }

    void FrequencySwitcher::velocity_ticks(double fmin, double fmax)
    {
        // Round velocity steps (1, 2 or 5 x 10^k km/s), about 8 ticks.
        // Frequency f (MHz) <-> v = c (f_rest - f) / f_rest
        const double c = 299792.458;
        vel_ticks.clear();
        vel_text.clear();
        vel_labels.clear();
        if (!(fmax > fmin) || !(rest_mhz > 0))
            return;
        // Observed velocities, plus the LSR correction (0 when off)
        double v1 = c * (rest_mhz - fmax) / rest_mhz + lsr_corr, v2 = c * (rest_mhz - fmin) / rest_mhz + lsr_corr;
        double raw = (v2 - v1) / 8, p = std::pow(10.0, std::floor(std::log10(raw))), step = p;
        for (double m : {1.0, 2.0, 5.0, 10.0})
            if (m * p >= raw)
            {
                step = m * p;
                break;
            }
        for (double v = std::ceil(v1 / step) * step; v <= v2 && vel_ticks.size() < 50; v += step)
        {
            double vv = std::fabs(v) < step * 1e-6 ? 0.0 : v;
            vel_ticks.push_back(rest_mhz * (1 - (vv - lsr_corr) / c));
            char t[32];
            snprintf(t, sizeof(t), step < 1 ? "%.1f" : "%.0f", vv);
            vel_text.push_back(t);
        }
        for (auto &t : vel_text)
            vel_labels.push_back(t.c_str());
    }

    void FrequencySwitcher::plot_folded(ImVec2 size, bool main_area)
    {
        // mtx must be held
        // Own mouse text: frequency once, velocity, value (ImPlot would show
        // the frequency twice, once per X axis)
        // Right click is ours (shows or hides the line window): no ImPlot
        // context menu nor box selection, the frequency axis fits by itself
        if (ImPlot::BeginPlot("##fswitch", size, ImPlotFlags_NoMouseText | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect))
        {
            ImPlot::SetupAxes("Frequency (MHz)", use_db ? "(Psys+dP)/Psys (dB)" : "dP/Psys x1e-3", ImPlotAxisFlags_AutoFit, auto_y ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None);
            float lo = use_db ? y_min_db : y_min, hi = use_db ? y_max_db : y_max;
            if (!auto_y)
                ImPlot::SetupAxisLimits(ImAxis_Y1, std::min(lo, hi), std::max(lo, hi) + 1e-6, ImPlotCond_Always);
            if (show_velocity && plot_x.size() > 1)
            {
                // Same frequency range as the bottom axis (last frame), labelled in km/s
                if (!(vel_lim[1] > vel_lim[0]))
                {
                    vel_lim[0] = plot_x.front();
                    vel_lim[1] = plot_x.back();
                }
                update_lsr();
                velocity_ticks(vel_lim[0], vel_lim[1]);
                ImPlot::SetupAxis(ImAxis_X2, lsr ? "Radial velocity, LSR (km/s)" : "Radial velocity, observer (km/s)", ImPlotAxisFlags_Opposite | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoMenus);
                ImPlot::SetupAxisLimits(ImAxis_X2, vel_lim[0], vel_lim[1], ImPlotCond_Always);
                if (vel_ticks.size() > 0)
                    ImPlot::SetupAxisTicks(ImAxis_X2, vel_ticks.data(), vel_ticks.size(), vel_labels.data());
            }
            if (show_velocity && plot_x.size() > 1)
            {
                ImPlotRect r = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                vel_lim[0] = r.X.Min;
                vel_lim[1] = r.X.Max;
            }
            if (plot_x.size() > 0)
                ImPlot::PlotLine("folded", plot_x.data(), plot_y.data(), plot_x.size());
            if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::GetIO().MouseDragMaxDistanceSqr[ImGuiMouseButton_Right] < 25)
                bl_show_lines = !bl_show_lines;
            if (bl_remove && bl_show_lines && plot_x.size() > 1 && rest_mhz > 0)
            {
                // Line window, draggable: the baseline is fitted outside it
                double fa = mhz_of_velocity(bl_vmin), fb = mhz_of_velocity(bl_vmax);
                ImVec4 col(1.0f, 0.75f, 0.2f, 0.9f);
                if (ImPlot::DragLineX(4201, &fa, col, 1.5f))
                    bl_vmin = rest_mhz > 0 ? (1 - fa / rest_mhz) * 299792.458 + lsr_corr : bl_vmin;
                if (ImPlot::DragLineX(4202, &fb, col, 1.5f))
                    bl_vmax = rest_mhz > 0 ? (1 - fb / rest_mhz) * 299792.458 + lsr_corr : bl_vmax;
                if (bl_vmin > bl_vmax)
                    std::swap(bl_vmin, bl_vmax);
            }
            if (ImPlot::IsPlotHovered())
            {
                ImPlotPoint mp = ImPlot::GetPlotMousePos(ImAxis_X1, ImAxis_Y1);
                // Snap to the nearest peak: start from the point of the curve
                // closest to the mouse, then climb, pixel column by pixel
                // column, to the top of that peak. A small peak between two
                // large ones is reached by pointing at it
                if (plot_x.size() > 1 && plot_x.size() == plot_y.size())
                {
                    const int R = 60;  // climbing range, in pixels
                    const int S = 15;  // start search range, in pixels
                    ImVec2 mpix = ImGui::GetMousePos();
                    // Highest point in each pixel column around the mouse
                    std::vector<long> col(2 * R + 1, -1);
                    for (int c = 0; c <= 2 * R; c++)
                    {
                        double xa = ImPlot::PixelsToPlot(ImVec2(mpix.x + c - R - 0.5f, mpix.y), ImAxis_X1, ImAxis_Y1).x;
                        double xb = ImPlot::PixelsToPlot(ImVec2(mpix.x + c - R + 0.5f, mpix.y), ImAxis_X1, ImAxis_Y1).x;
                        if (xa > xb)
                            std::swap(xa, xb);
                        size_t ia = std::lower_bound(plot_x.begin(), plot_x.end(), xa) - plot_x.begin();
                        size_t ib = std::upper_bound(plot_x.begin(), plot_x.end(), xb) - plot_x.begin();
                        for (size_t i = ia; i < ib; i++)
                            if (!std::isnan(plot_y[i]) && (col[c] < 0 || plot_y[i] > plot_y[col[c]]))
                                col[c] = i;
                    }
                    // Start: the column whose point is nearest to the mouse
                    int s = -1;
                    float dbest = 0;
                    for (int c = R - S; c <= R + S; c++)
                    {
                        if (col[c] < 0)
                            continue;
                        ImVec2 p = ImPlot::PlotToPixels(ImPlotPoint(plot_x[col[c]], plot_y[col[c]]), ImAxis_X1, ImAxis_Y1);
                        float dd = (p.x - mpix.x) * (p.x - mpix.x) + (p.y - mpix.y) * (p.y - mpix.y);
                        if (s < 0 || dd < dbest)
                            s = c, dbest = dd;
                    }
                    // Climb while a column within 2 pixels is higher
                    while (s >= 0)
                    {
                        int next = s;
                        for (int c = std::max(0, s - 2); c <= std::min(2 * R, s + 2); c++)
                            if (col[c] >= 0 && plot_y[col[c]] > plot_y[col[next]])
                                next = c;
                        if (next == s)
                            break;
                        s = next;
                    }
                    size_t best = s >= 0 ? (size_t)col[s] : plot_x.size();
                    if (best < plot_x.size())
                    {
                        mp = ImPlotPoint(plot_x[best], plot_y[best]);
                        ImVec2 p = ImPlot::PlotToPixels(mp, ImAxis_X1, ImAxis_Y1);
                        ImPlot::GetPlotDrawList()->AddCircle(p, 5.0f, ImGui::GetColorU32(ImGuiCol_Text), 0, 1.5f);
                    }
                }
                char txt[128];
                int n = snprintf(txt, sizeof(txt), "%.1f MHz", mp.x);
                if (show_velocity && rest_mhz > 0)
                    n += snprintf(txt + n, sizeof(txt) - n, "   %.1f km/s%s", 299792.458 * (rest_mhz - mp.x) / rest_mhz + lsr_corr, lsr ? " LSR" : "");
                snprintf(txt + n, sizeof(txt) - n, use_db ? "   %.3f dB" : "   %.2f", mp.y);
                ImVec2 pos = ImPlot::GetPlotPos(), sz = ImPlot::GetPlotSize(), ts = ImGui::CalcTextSize(txt);
                ImPlot::GetPlotDrawList()->AddText(ImVec2(pos.x + sz.x - ts.x - 8, pos.y + sz.y - ts.y - 6), ImGui::GetColorU32(ImGuiCol_Text), txt);
            }
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

    static std::string utc_now(bool for_file)
    {
        auto now = std::chrono::system_clock::now();
        time_t t = std::chrono::system_clock::to_time_t(now);
        int ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000;
        struct tm tmv;
#ifdef _WIN32
        gmtime_s(&tmv, &t);
#else
        gmtime_r(&t, &tmv);
#endif
        char ts[40];
        if (for_file)
            strftime(ts, sizeof(ts), "%Y-%m-%d_%H-%M-%S", &tmv);
        else
        {
            char b[32];
            strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tmv);
            snprintf(ts, sizeof(ts), "%s.%03dZ", b, ms);
        }
        return ts;
    }

    void FrequencySwitcher::open_cycles(bool append)
    {
        // mtx must be held, after init_run() set the run parameters
        cyc_a.assign(nbins, 0);
        cyc_b.assign(nbins, 0);
        cyc_na = cyc_nb = 0;
        cyc_last_ph = -1;
        cyc_error.clear();
        if (append)
        {
            // Resumed run: same columns, the cycle numbers go on
            cyc_file.open(cyc_path, std::ios::app);
            if (cyc_file)
                return;
            cyc_file.clear();
            logger->warn("Frequency switching: cannot reopen " + cyc_path + ", starting a new cycles file");
        }
        cyc_count = 0;
        cyc_path = cycles_dir + "/fswitch_cycles_" + std::to_string((uint64_t)std::llround(f0)) + "Hz_" + utc_now(true) + "Z.csv";
        cyc_file.open(cyc_path);
        if (!cyc_file)
        {
            cyc_error = "cannot write " + cyc_path;
            logger->error("Frequency switching: " + cyc_error);
            cyc_path.clear();
            return;
        }
        cyc_file << "# SatDump frequency switching, one folded spectrum per cycle (experimental)\n";
        cyc_file << "# f0_hz=" << (uint64_t)std::llround(f0) << " samplerate=" << samplerate << " fft_size=" << nbins << " offset_hz=+/-" << delta_eff
                 << " offset_bins=" << dk << " dwell_s=" << run_dwell << " settle_s=" << run_settle << " switching=" << (remote ? "server" : "client") << "\n";
        cyc_file << "# value = dP/Psys, linear, not smoothed: ((A-B)/B (j-d) + (B-A)/A (j+d)) / 2 with A, B the mean power of this cycle\n";
        cyc_file << "# utc = end of the cycle (client clock), frames = FFT frames in A and B. Columns after frames_b: sky frequency in Hz (next line)\n";
        cyc_file << "utc,unix_s,cycle,frames_a,frames_b";
        cyc_file.precision(12);
        for (int j = dk; j < nbins - dk; j++)
            cyc_file << "," << (uint64_t)std::llround(f0 + ((double)j / nbins - 0.5) * samplerate);
        cyc_file << "\n";
        cyc_file.precision(6);
    }

    void FrequencySwitcher::close_cycles()
    {
        // mtx must be held
        if (cyc_file.is_open())
        {
            cyc_file.close();
            // Nothing in it (stopped at once, or the server refused): remove it
            if (cyc_count == 0)
            {
                std::error_code ec;
                std::filesystem::remove(cyc_path, ec);
                cyc_path.clear();
                return;
            }
            logger->info("Frequency switching: %llu cycles saved to %s", (unsigned long long)cyc_count, cyc_path.c_str());
        }
    }

    void FrequencySwitcher::cycle_add(int ph, const float *lin_f, const double *lin_d, double weight)
    {
        // mtx must be held
        if (!cyc_file.is_open() || (int)cyc_a.size() != nbins)
            return;
        // A new A after a B: the previous cycle is complete
        if (ph == 0 && cyc_last_ph == 1)
            write_cycle();
        cyc_last_ph = ph;
        std::vector<double> &c = ph == 0 ? cyc_a : cyc_b;
        for (int k = 0; k < nbins; k++)
            c[k] += (lin_f ? (double)lin_f[k] : lin_d[k]) * weight;
        (ph == 0 ? cyc_na : cyc_nb) += (uint64_t)weight;
    }

    void FrequencySwitcher::write_cycle()
    {
        // mtx must be held
        if (cyc_na > 0 && cyc_nb > 0)
        {
            double unix_s = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
            cyc_file << utc_now(false) << ",";
            cyc_file.precision(13);
            cyc_file << unix_s;
            cyc_file.precision(6);
            cyc_file << "," << cyc_count << "," << cyc_na << "," << cyc_nb;
            // Means of this cycle, with the interference found in the whole averages removed
            std::vector<double> ca(nbins), cb(nbins);
            for (int k = 0; k < nbins; k++)
            {
                ca[k] = cyc_a[k] / cyc_na;
                cb[k] = cyc_b[k] / cyc_nb;
            }
            if (rfi_remove)
            {
                rfi_fill(ca, rfi_a);
                rfi_fill(cb, rfi_b);
            }
            for (int j = dk; j < nbins - dk; j++)
            {
                double aa = ca[j - dk], ba = cb[j - dk];
                double ab = ca[j + dk], bb = cb[j + dk];
                cyc_file << ",";
                if (ba > 0 && ab > 0)
                    cyc_file << 0.5 * ((aa - ba) / ba + (bb - ab) / ab);
            }
            cyc_file << "\n";
            cyc_file.flush();
            cyc_count++;
            if (!cyc_file)
            {
                cyc_error = "write error, cycles file closed";
                logger->error("Frequency switching: " + cyc_error);
                cyc_file.close();
            }
        }
        std::fill(cyc_a.begin(), cyc_a.end(), 0);
        std::fill(cyc_b.begin(), cyc_b.end(), 0);
        cyc_na = cyc_nb = 0;
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
        out << "# offset_hz=+/-" << delta_eff << " offset_bins=" << dk << " dwell_s=" << run_dwell << " settle_s=" << run_settle << "\n";
        out << "# switching=" << (remote ? "server" : "client");
        if (remote && remote_iq_file.size())
            out << " iq_on_server=" << remote_iq_file;
        out << "\n";
        out << "# cycles=" << cycles << " frames_a=" << n_a << " frames_b=" << n_b << "\n";
        out << "# A = mean linear power tuned at f0+offset, B = tuned at f0-offset, indexed by FFT bin\n";
        out << "# q = (A-B)/B ; folded(j) = ((A-B)/B (j-d) + (B-A)/A (j+d))/2, sky frequency rf_hz, empty near the edges\n";
        compute(); // current baseline (bl_curve), indexed by j - dk
        bool blc = bl_remove && (int)bl_curve.size() == nbins - 2 * dk;
        std::vector<double> ma, mb;
        mean_spectra(ma, mb);
        bool rfi = rfi_remove && (int)rfi_a.size() == nbins && (int)rfi_b.size() == nbins;
        if (rfi)
            out << "# interference removed: a_lin and b_lin interpolated where rfi_a / rfi_b = 1 (" << rfi_count_a << " / " << rfi_count_b
                << " bins, threshold " << rfi_sigma << " sigma, local median over " << rfi_width_khz << " kHz)\n";
        else
            out << "# interference removal off\n";
        update_lsr();
        out << "# v_kms = radial velocity from rf_hz, radio convention, rest frequency " << std::to_string(rest_mhz) << " MHz, observer frame\n";
        if (lsr)
            out << "# v_lsr_kms = v_kms + " << lsr_corr << " (LSR correction at save time: Earth rotation " << lsr_rot << ", Earth orbit " << lsr_orb << ", Sun " << lsr_sun
                << "), source " << (lsr_frame == 1 ? "l=" + std::to_string(src_l) + " b=" + std::to_string(src_b) : "RA=" + std::to_string(src_ra_h) + "h Dec=" + std::to_string(src_dec))
                << ", QTH " << qth_lat << " " << qth_lon << " " << qth_alt << " m\n";
        if (blc)
            out << "# folded_bl = folded minus a Legendre polynomial of order " << bl_order << " fitted outside " << bl_vmin << " .. " << bl_vmax << " km/s"
                << (lsr ? " (LSR)" : " (observer)") << ", residual rms " << bl_rms << " over " << bl_points << " points\n";
        out << "bin,if_offset_hz,rf_hz,a_lin,b_lin,q,folded,rfi_a,rfi_b,v_kms" << (lsr ? ",v_lsr_kms" : "") << (blc ? ",folded_bl" : "") << "\n";
        out.precision(10);
        for (int j = 0; j < nbins; j++)
        {
            double off = ((double)j / nbins - 0.5) * samplerate;
            double a = ma[j], b = mb[j];
            double q = b > 0 ? (a - b) / b : 0;
            out << j << "," << off << "," << (uint64_t)std::llround(f0 + off) << "," << a << "," << b << "," << q << ",";
            if (j >= dk && j < nbins - dk)
            {
                double aa = ma[j - dk], ba = mb[j - dk];
                double ab = ma[j + dk], bb = mb[j + dk];
                double q12 = ba > 0 ? (aa - ba) / ba : 0; // o1 - o2, reference B
                double q21 = ab > 0 ? (bb - ab) / ab : 0; // o2 - o1, reference A
                out << 0.5 * (q12 + q21);
            }
            out << "," << (rfi ? (int)rfi_a[j] : 0) << "," << (rfi ? (int)rfi_b[j] : 0) << "," << 299792.458 * (rest_mhz * 1e6 - (f0 + off)) / (rest_mhz * 1e6);
            if (lsr)
                out << "," << 299792.458 * (rest_mhz * 1e6 - (f0 + off)) / (rest_mhz * 1e6) + lsr_corr;
            if (blc)
            {
                out << ",";
                if (j >= dk && j < nbins - dk && !std::isnan(folded_raw[j - dk]))
                    out << folded_raw[j - dk] - bl_curve[j - dk];
            }
            out << "\n";
        }
        return path;
    }
} // namespace satdump
