#pragma once

/*
 * Experimental frequency switching for the recorder.
 *
 * Meant for very stable signals (radio astronomy, weak continuum or line
 * sources), not for satellites. The SDR is tuned alternately to f0 + D and
 * f0 - D. The FFT frames of each position are averaged in linear power
 * (A and B). The receiver bandpass (IF filter, ripple) stays at the same
 * place in both positions while the sky signal moves by 2D, so
 *
 *     o1 - o2 = (A(k) - B(k)) / B(k)     (B is the reference)
 *     o2 - o1 = (B(k) - A(k)) / A(k)     (A is the reference)
 *
 * remove the bandpass shape and give the signal in units of the system
 * power. In o1 - o2 the line is positive at k = j - d, in o2 - o1 at
 * k = j + d; both are shifted back onto the sky frequency axis and averaged
 * (in-band frequency switching, sqrt(2) gain in SNR):
 *
 *     T(j) = [ (A-B)/B (j - d) + (B-A)/A (j + d) ] / 2    with d = D in bins
 *
 * Each half uses its own reference, so a strong line is not underestimated.
 * The result keeps two negative ghosts of half amplitude at +/- 2D.
 */

#include "imgui/imgui.h"
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace satdump
{
    class FrequencySwitcher
    {
    public:
        // Settings (only applied at start)
        double offset_khz = 450;
        // Offset as a share of the largest one allowed (30 % of the samplerate):
        // the centre then stays inside the flat part of the IF filter
        float offset_pct = 60;
        static double max_offset_hz(double samplerate) { return 0.3 * samplerate; }
        float dwell_s = 2.0f;
        float settle_s = 0.5f;
        // Delay between a retune and its samples (server, network, buffers).
        // Measured during the first cycles when auto_delay is on.
        bool auto_delay = true;
        float delay_s = 0.0f;

        // Display
        int smoothing = 1;
        bool show_window = true;          // folded result
        bool show_unfolded_window = false; // o1 - o2 and o2 - o1, before folding

        // Vertical scale, in units of 1e-3 of the system power
        bool auto_y = false;
        float y_min = -5.0f, y_max = 5.0f;
        // dB display: 10 log10(1 + dP/Psys), own scale
        bool use_db = false;
        float y_min_db = -0.5f, y_max_db = 0.5f;
        // Bins more than this below the median level (band edges) are hidden
        float edge_db = 3.0f;

        // Tunes the SDR without changing the recorder frequency
        std::function<void(double)> set_sdr_frequency = [](double) {};

        // Main display while switching: each frame is moved back to the sky
        // frequency (by -/+ the offset) then averaged like the normal FFT
        float display_avg = 10;
        std::function<void(float *)> on_display = [](float *) {};
        float *display_buffer() { return disp.data(); }

    public:
        ~FrequencySwitcher() { stop(); }

        // frame_dt: duration of signal used by one FFT frame (seconds)
        // init_display: spectrum shown when switching starts (dB, fft_size values), may be null
        void start(double center_hz, double samplerate, int fft_size, double frame_dt, const float *init_display = nullptr);
        void stop();
        void reset();
        bool is_running() { return running; }

        // Called from the FFT thread with the power spectrum in dB
        void feed(const float *db, int size);

        // Draws the floating windows (call from the UI thread). The folded
        // window is only shown when not running: the recorder then draws the
        // folded result in place of its spectrum, with draw_plot()
        void draw_window();
        void draw_plot(ImVec2 size);
        // True once after the user deselected "folded" in the legend of the
        // main plot: the recorder then stops switching
        bool take_stop_request() { return stop_requested.exchange(false); }
        // Draws status lines (cycles, integration time)
        void draw_status();

        // Saves the result as CSV in the given directory, returns the file name
        std::string save_csv(std::string directory);

    private:
        void run();
        void compute();
        int phase_of(double fstart, double fend, double delay);
        // Best delay in [lo, hi), returns false if no data. contrast = best / mean score
        bool search_delay(double lo, double hi, double &best_d, double &contrast);
        void estimate_delay();
        void track_delay();
        void add_slot(const float *db, double fstart, double fend);

        std::mutex mtx;
        std::atomic<bool> running{false};
        std::thread thread;
        std::mutex cv_mtx;
        std::condition_variable cv;

        // Set at start
        double f0 = 0;
        double samplerate = 0;
        int nbins = 0;
        int dk = 0;
        double delta_eff = 0;
        double run_dwell = 0, run_settle = 0, run_delay = 0;

        // Retune events: time (s) and phase (0 = f0 + D (A), 1 = f0 - D (B), -1 = back to f0)
        std::deque<std::pair<double, int>> events;

        // Frames are dated by their arrival time: samples lost on the way
        // (server, network) do not shift the dates. The delay between a
        // retune and the arrival of its samples is measured, then tracked.
        double frame_dt = 0;
        double t_first = 0;
        uint64_t n_frames = 0;
        double lag_s = 0;
        double speed = 1;
        // Sliding measurement over the last seconds: (arrival time, frame count)
        std::deque<std::pair<double, uint64_t>> rate_hist;
        double speed_now = 1;

        // Delay measurement: short slots of decimated spectra
        struct Slot
        {
            double t_start, t_end;
            std::vector<double> sum;
            int count = 0;
        };
        static const int CAL_CHANNELS = 256;
        static const int CAL_CYCLES = 3;
        static const int TRACK_CYCLES = 4;
        bool calibrating = false;
        double slot_len = 0.05;
        std::vector<Slot> slots;
        double measured_delay = -1;
        double delay_score_ratio = 0;

        // Accumulators (linear power)
        std::vector<double> sum_a, sum_b;
        uint64_t n_a = 0, n_b = 0;
        uint64_t cycles = 0;

        // Computed for display
        std::vector<double> plot_x, plot_y;
        // Unfolded: o1 - o2 = (A-B)/B and o2 - o1 = (B-A)/A, receiver axis
        std::vector<double> unf_x, unf_12, unf_21;
        void smooth_scale(const std::vector<double> &y, std::vector<double> &out);
        int computed_frame = -1; // compute() once per UI frame
        std::atomic<bool> stop_requested{false};
        void compute_once();
        void plot_folded(ImVec2 size, bool main_area);

        // Main display spectrum (dB), on the sky frequency axis
        std::vector<float> disp;
        // Recentred averages of each position, linear power (0 = no data)
        std::vector<float> disp_a, disp_b;
        void update_display(const float *db, int ph);
    };
} // namespace satdump
