#pragma once

/*
 * Frequency switching done by the server (experimental), see fswitch_proto.h.
 *
 * The stream thread feeds every sample to process(). Retunes are counted in
 * samples: after each retune, settle_s of samples are dropped, then dwell_s
 * of samples are turned into FFT frames and averaged. The mean spectrum of
 * the dwell is sent to the client, then the other frequency is requested.
 * The retune itself is done by the stream thread, outside the stream lock
 * (see streaming.cpp), so no sample is processed while it happens.
 *
 * Optionally the IQ is written to disk on the server with a CSV of events
 * (sample index of each retune and of each dwell), for post-processing.
 */

#include "common/dsp/complex.h"
#include "nlohmann/json.hpp"
#include <cstdint>
#include <fftw3.h>
#include <fstream>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

class ServerFSwitch
{
public:
    ~ServerFSwitch();

    // Directory for the IQ files on the server
    std::string iq_dir;

    // From the packet handler. Returns the status to send back
    nlohmann::json start(const nlohmann::json &settings, uint64_t samplerate, bool source_started);
    // From the packet handler or when the source stops. With tune_back,
    // asks for the frequency to be set back to freq_back
    nlohmann::json stop(std::string reason = "", bool tune_back = true);
    bool active();

    // From the stream thread, with the stream lock held
    void process(const complex_t *samples, int n);
    // From the retune thread (see streaming.cpp), with source_mtx held:
    // retune to do. retune_wanted() is a cheap check before locking.
    // The stream thread keeps reading meanwhile (some drivers, like
    // librtlsdr, cannot retune while their callback waits for the reader):
    // samples read until retune_done() are dropped.
    bool retune_wanted() { return retune_pending; }
    // Waits up to timeout_ms for a retune request
    void wait_retune(int timeout_ms);
    bool take_retune(double &freq);
    void retune_done();
    // With source_mtx held, when the client tunes by itself
    void cancel_retune();

private:
    std::mutex mtx;
    bool running = false;

    // Settings
    double freq_a = 0, freq_b = 0, freq_back = 0;
    double sky_f0 = 0, offset_hz = 0;
    uint64_t samplerate = 0;
    int fft_size = 0;
    uint64_t settle_n = 0, dwell_n = 0;

    // State
    enum State
    {
        WAIT_RETUNE,
        SETTLING,
        DWELLING,
    };
    State state = WAIT_RETUNE;
    int phase = 0; // 0 = A, 1 = B
    std::atomic<bool> retune_pending{false};
    std::condition_variable retune_cv;
    double retune_freq = 0;
    uint64_t pos = 0;         // samples since the retune
    uint64_t dwell_index = 0; // dwells sent
    uint64_t total_n = 0;     // samples seen since start

    // FFT, same as dsp::FFTPanBlock (Nuttall window, shifted, |X/N|^2)
    std::vector<float> taps;
    fftwf_complex *fftw_in = nullptr, *fftw_out = nullptr;
    fftwf_plan plan = nullptr;
    int in_fft = 0;
    std::vector<double> acc;
    int frames = 0;
    void free_fft();
    void fft_frame();
    void end_dwell();

    // IQ recording
    bool save_iq = false;
    int iq_bits = 16;
    std::ofstream iq_file, ev_file;
    std::string iq_path;
    uint64_t iq_bytes = 0;
    std::vector<uint8_t> iq_buf;
    std::string iq_error;
    void open_iq();
    void close_iq();
    void write_iq(const complex_t *s, int n);
    void event(const char *name, int ph, double freq);

    nlohmann::json status_locked(std::string error = "");
};

extern ServerFSwitch server_fswitch;
