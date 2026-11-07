#include "fswitch.h"
#include "fswitch_proto.h"
#include "main.h"
#include "remote.h"
#include "logger.h"
#include "common/dsp/window/window.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <chrono>

ServerFSwitch server_fswitch;

// Do not start an IQ file with less than this free, stop it below the second
static const uint64_t IQ_MIN_FREE_START = 2000000000ULL;
static const uint64_t IQ_MIN_FREE_RUN = 500000000ULL;

ServerFSwitch::~ServerFSwitch()
{
    std::lock_guard<std::mutex> lk(mtx);
    close_iq();
    free_fft();
}

bool ServerFSwitch::active()
{
    std::lock_guard<std::mutex> lk(mtx);
    return running;
}

void ServerFSwitch::free_fft()
{
    if (plan != nullptr)
        fftwf_destroy_plan(plan);
    if (fftw_in != nullptr)
        fftwf_free(fftw_in);
    if (fftw_out != nullptr)
        fftwf_free(fftw_out);
    plan = nullptr;
    fftw_in = fftw_out = nullptr;
}

nlohmann::json ServerFSwitch::status_locked(std::string error)
{
    nlohmann::json st;
    st["running"] = running;
    st["samplerate"] = samplerate;
    st["fft_size"] = fft_size;
    st["dwells"] = dwell_index;
    st["save_iq"] = save_iq;
    st["iq_file"] = iq_path;
    st["iq_bytes"] = iq_bytes;
    if (iq_error.size())
        st["iq_error"] = iq_error;
    if (error.size())
        st["error"] = error;
    return st;
}

nlohmann::json ServerFSwitch::start(const nlohmann::json &s, uint64_t sr, bool source_started)
{
    std::lock_guard<std::mutex> lk(mtx);
    if (running)
        return status_locked("already running");
    if (!source_started || sr == 0)
        return status_locked("the source is not started");

    try
    {
        freq_a = s["freq_a"].get<double>();
        freq_b = s["freq_b"].get<double>();
        freq_back = s["freq_back"].get<double>();
        sky_f0 = s.value("sky_f0", 0.0);
        offset_hz = s.value("offset_hz", 0.0);
        fft_size = s["fft_size"].get<int>();
        double dwell_s = s["dwell_s"].get<double>();
        double settle_s = s["settle_s"].get<double>();
        save_iq = s.value("save_iq", false);
        iq_bits = s.value("iq_bits", 16);

        // The TCP buffers hold 3 MB: 2^18 floats = 1 MB per dwell packet
        if (fft_size < 16 || fft_size > (1 << 18))
            return status_locked("bad FFT size");
        if (dwell_s < 0.1 || settle_s < 0)
            return status_locked("bad dwell or settle time");
        if (iq_bits != 8 && iq_bits != 16 && iq_bits != 32)
            iq_bits = 16;

        samplerate = sr;
        settle_n = (uint64_t)std::llround(settle_s * sr);
        dwell_n = (uint64_t)std::llround(dwell_s * sr);
        if (dwell_n < (uint64_t)fft_size)
            return status_locked("dwell shorter than one FFT");
    }
    catch (std::exception &e)
    {
        return status_locked(std::string("bad settings: ") + e.what());
    }

    // FFT, as in dsp::FFTPanBlock: Nuttall window with alternating sign
    // (spectrum shifted, DC in the middle)
    free_fft();
    taps.resize(fft_size);
    for (int i = 0; i < fft_size; i++)
        taps[i] = dsp::window::nuttall(i, fft_size - 1) * ((i % 2) ? 1.0f : -1.0f);
    fftw_in = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * fft_size);
    fftw_out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * fft_size);
    plan = fftwf_plan_dft_1d(fft_size, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);
    acc.assign(fft_size, 0.0);
    in_fft = 0;
    frames = 0;

    phase = 0;
    pos = 0;
    dwell_index = 0;
    total_n = 0;
    iq_bytes = 0;
    iq_path.clear();
    iq_error.clear();
    if (save_iq)
        open_iq();

    // First retune: to A
    state = WAIT_RETUNE;
    retune_freq = freq_a;
    retune_pending = true;
    running = true;
    retune_cv.notify_all();

    logger->info("Frequency switching: A %.0f Hz / B %.0f Hz, FFT %d, settle %llu dwell %llu samples%s", freq_a, freq_b, fft_size,
                 (unsigned long long)settle_n, (unsigned long long)dwell_n, save_iq ? (", IQ to " + iq_path).c_str() : "");
    return status_locked();
}

nlohmann::json ServerFSwitch::stop(std::string reason, bool tune_back)
{
    std::lock_guard<std::mutex> lk(mtx);
    if (!running)
        return status_locked();
    running = false;
    event("stop", -1, freq_back);
    close_iq();
    free_fft();
    retune_freq = freq_back;
    retune_pending = tune_back;
    retune_cv.notify_all();
    logger->info("Frequency switching stopped after %llu dwells%s", (unsigned long long)dwell_index, reason.size() ? (" (" + reason + ")").c_str() : "");
    return status_locked(reason);
}

void ServerFSwitch::wait_retune(int timeout_ms)
{
    std::unique_lock<std::mutex> lk(mtx);
    retune_cv.wait_for(lk, std::chrono::milliseconds(timeout_ms), [this] { return retune_pending.load(); });
}

bool ServerFSwitch::take_retune(double &freq)
{
    std::lock_guard<std::mutex> lk(mtx);
    if (!retune_pending)
        return false;
    retune_pending = false;
    freq = retune_freq;
    return true;
}

void ServerFSwitch::cancel_retune()
{
    std::lock_guard<std::mutex> lk(mtx);
    retune_pending = false;
}

void ServerFSwitch::retune_done()
{
    std::lock_guard<std::mutex> lk(mtx);
    if (!running || state != WAIT_RETUNE)
        return;
    // Every sample processed from now on was read after the retune
    // returned. Those still in the driver buffers are dropped by the
    // settle time.
    event("retune", phase, phase == 0 ? freq_a : freq_b);
    state = settle_n > 0 ? SETTLING : DWELLING;
    if (state == DWELLING)
        event("dwell_start", phase, phase == 0 ? freq_a : freq_b);
    pos = 0;
    in_fft = 0;
}

void ServerFSwitch::process(const complex_t *s, int n)
{
    std::lock_guard<std::mutex> lk(mtx);
    if (!running)
        return;

    if (save_iq)
        write_iq(s, n);

    int i = 0;
    while (i < n)
    {
        if (state == WAIT_RETUNE)
        {
            // Old frequency, or the retune is not done yet
            total_n += n - i;
            return;
        }

        if (state == SETTLING)
        {
            uint64_t left = settle_n - pos;
            int take = (int)std::min<uint64_t>(left, n - i);
            pos += take;
            i += take;
            total_n += take;
            if (pos >= settle_n)
            {
                state = DWELLING;
                event("dwell_start", phase, phase == 0 ? freq_a : freq_b);
            }
            continue;
        }

        // DWELLING: whole FFT frames only
        uint64_t left = settle_n + dwell_n - pos;
        int take = (int)std::min<uint64_t>(std::min<uint64_t>(left, n - i), fft_size - in_fft);
        memcpy(&fftw_in[in_fft], &s[i], take * sizeof(complex_t));
        in_fft += take;
        pos += take;
        i += take;
        total_n += take;
        if (in_fft == fft_size)
        {
            fft_frame();
            in_fft = 0;
        }
        if (pos >= settle_n + dwell_n)
        {
            event("dwell_end", phase, phase == 0 ? freq_a : freq_b);
            end_dwell();
        }
    }
}

void ServerFSwitch::fft_frame()
{
    // mtx held. fftw_in holds fft_size samples
    for (int k = 0; k < fft_size; k++)
    {
        fftw_in[k][0] *= taps[k];
        fftw_in[k][1] *= taps[k];
    }
    fftwf_execute(plan);
    // Same scale as volk_32fc_s32f_power_spectrum_32f(..., fft_size, ...)
    double norm = 1.0 / ((double)fft_size * fft_size);
    for (int k = 0; k < fft_size; k++)
    {
        double re = fftw_out[k][0], im = fftw_out[k][1];
        acc[k] += (re * re + im * im) * norm + 1e-20;
    }
    frames++;
}

void ServerFSwitch::end_dwell()
{
    // mtx held
    if (frames > 0)
    {
        remote_sdr::FSwitchDwell d;
        d.phase = phase;
        d.frames = frames;
        d.dwell_index = dwell_index;
        d.mean.resize(fft_size);
        for (int k = 0; k < fft_size; k++)
            d.mean[k] = acc[k] / frames;
        dsp::remote::sendPacketWithVector(tcp_server, dsp::remote::PKT_TYPE_FSWITCH_DWELL, remote_sdr::encode_fswitch_dwell(d));
    }
    dwell_index++;
    std::fill(acc.begin(), acc.end(), 0.0);
    frames = 0;
    in_fft = 0;

    // A status now and then, for the IQ file size
    if (save_iq && (dwell_index % 2) == 0)
        dsp::remote::sendPacketWithVector(tcp_server, dsp::remote::PKT_TYPE_FSWITCH_STATUS, nlohmann::json::to_cbor(status_locked()));

    phase ^= 1;
    state = WAIT_RETUNE;
    retune_freq = phase == 0 ? freq_a : freq_b;
    retune_pending = true;
    retune_cv.notify_all();
}

void ServerFSwitch::open_iq()
{
    // mtx held
    try
    {
        std::filesystem::create_directories(iq_dir);
        auto sp = std::filesystem::space(iq_dir);
        if (sp.available < IQ_MIN_FREE_START)
        {
            iq_error = "not enough free space in " + iq_dir;
            save_iq = false;
            return;
        }

        time_t t = time(nullptr);
        struct tm tmv;
#ifdef _WIN32
        gmtime_s(&tmv, &t);
#else
        gmtime_r(&t, &tmv);
#endif
        char ts[32];
        strftime(ts, sizeof(ts), "%Y-%m-%d_%H-%M-%S", &tmv);
        std::string base = iq_dir + "/fswitch_" + std::to_string((uint64_t)std::llround(sky_f0)) + "Hz_" + ts + "Z";
        std::string ext = iq_bits == 8 ? ".cs8" : iq_bits == 16 ? ".cs16" : ".cf32";
        iq_path = base + ext;

        iq_file.open(iq_path, std::ios::binary);
        ev_file.open(base + ".csv");
        if (!iq_file || !ev_file)
        {
            iq_error = "cannot write in " + iq_dir;
            close_iq();
            save_iq = false;
            return;
        }

        ev_file << "# SatDump frequency switching, IQ written by the server (experimental)\n";
        ev_file << "# iq_file=" << std::filesystem::path(iq_path).filename().string() << " format=" << ext.substr(1)
                << (iq_bits == 32 ? "" : " (interleaved I/Q integers, full scale = 1.0)") << "\n";
        ev_file << "# samplerate=" << samplerate << " sky_f0_hz=" << (uint64_t)std::llround(sky_f0) << " offset_hz=+/-" << offset_hz << "\n";
        ev_file << "# sdr_freq_a_hz=" << (uint64_t)std::llround(freq_a) << " sdr_freq_b_hz=" << (uint64_t)std::llround(freq_b)
                << " settle_samples=" << settle_n << " dwell_samples=" << dwell_n << " fft_size=" << fft_size << "\n";
        ev_file << "# start_utc=" << ts << "Z\n";
        ev_file << "# phase 0 = A (sky f0 + offset), 1 = B (sky f0 - offset). A dwell is [dwell_start, dwell_end).\n";
        ev_file << "# retune = first sample read after the SDR was retuned (samples just after it may still be at the old frequency)\n";
        ev_file << "sample_index,event,phase,sdr_freq_hz\n";
        iq_buf.resize(1 << 20);
    }
    catch (std::exception &e)
    {
        iq_error = std::string("IQ file: ") + e.what();
        close_iq();
        save_iq = false;
    }
}

void ServerFSwitch::close_iq()
{
    // mtx held
    if (iq_file.is_open())
        iq_file.close();
    if (ev_file.is_open())
    {
        ev_file << total_n << ",end,-1,0\n";
        ev_file.close();
    }
    std::vector<uint8_t>().swap(iq_buf);
}

void ServerFSwitch::event(const char *name, int ph, double freq)
{
    // mtx held
    if (ev_file.is_open())
        ev_file << total_n << "," << name << "," << ph << "," << (uint64_t)std::llround(freq) << "\n";
}

void ServerFSwitch::write_iq(const complex_t *s, int n)
{
    // mtx held. total_n is not advanced here: process() does it
    if (!iq_file.is_open())
        return;

    size_t bytes;
    if (iq_bits == 32)
    {
        bytes = n * sizeof(complex_t);
        iq_file.write((const char *)s, bytes);
    }
    else
    {
        int bps = iq_bits / 8;
        bytes = (size_t)n * 2 * bps;
        if (iq_buf.size() < bytes)
            iq_buf.resize(bytes);
        float scale = iq_bits == 8 ? 127.0f : 32767.0f;
        const float *f = (const float *)s;
        for (int i = 0; i < 2 * n; i++)
        {
            float v = std::max(-1.0f, std::min(1.0f, f[i])) * scale;
            if (bps == 1)
                ((int8_t *)iq_buf.data())[i] = (int8_t)std::lround(v);
            else
                ((int16_t *)iq_buf.data())[i] = (int16_t)std::lround(v);
        }
        iq_file.write((const char *)iq_buf.data(), bytes);
    }
    iq_bytes += bytes;

    if (!iq_file)
    {
        iq_error = "write error, IQ file closed";
        logger->error("Frequency switching: IQ " + iq_error);
        close_iq();
        save_iq = false;
        return;
    }

    // Check the free space now and then (about every 256 MB)
    if ((iq_bytes >> 28) != ((iq_bytes - bytes) >> 28))
    {
        std::error_code ec;
        auto sp = std::filesystem::space(iq_dir, ec);
        if (!ec && sp.available < IQ_MIN_FREE_RUN)
        {
            iq_error = "disk almost full, IQ file closed";
            logger->warn("Frequency switching: " + iq_error);
            close_iq();
            save_iq = false;
        }
    }
}
