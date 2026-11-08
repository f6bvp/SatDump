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
#include "common/dsp_source_sink/remote_fswitch.h"
#include "nlohmann/json.hpp"
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <fstream>
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
        static constexpr float DWELL_MIN_S = 0.1f, DWELL_MAX_S = 10.0f, SETTLE_MAX_S = 5.0f;

        // One folded spectrum per cycle (A then B) written to a CSV file in
        // cycles_dir, for post-processing. Set before start
        bool save_cycles = false;
        std::string cycles_dir;

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

        // Radial velocity scale above the folded plot, radio convention:
        // v = c (f_rest - f) / f_rest, positive when the source recedes.
        // Relative to the observer (no correction to the LSR)
        bool show_velocity = true;
        double rest_mhz = 1420.405752; // HI; OH: 1665.402, 1667.359

        // Velocity relative to the local standard of rest (kinematic LSR,
        // 20 km/s toward 18h03m50.3s +30d00m17s J2000): observer velocity
        // projected on the source = Earth rotation + Earth orbit + Sun.
        // Computed for the current time; agrees with astropy to 0.02 km/s
        // Baseline: a polynomial (order bl_order) fitted to the folded
        // spectrum outside the line window [bl_vmin, bl_vmax] (km/s, same
        // frame as the velocity axis), then subtracted everywhere. Removes
        // the slow ripple left by the parts of the receiver response that
        // follow the RF frequency (antenna, filters, cables)
        bool bl_remove = false;
        int bl_order = 4;
        float bl_vmin = -150.0f, bl_vmax = 150.0f;
        bool bl_show_lines = true; // the two lines of the window (right click on the plot toggles)
        double bl_rms = 0;  // residual outside the window (dP/Psys)
        int bl_points = 0;
        double mhz_of_velocity(double v) { return rest_mhz * (1 - (v - lsr_corr) / 299792.458); }

        bool lsr = false;
        int lsr_frame = 0;               // 0 = RA/Dec J2000, 1 = galactic l/b
        double src_ra_h = 0.712, src_dec = 41.269; // M31
        double src_l = 120.0, src_b = 0.0;
        // Correction to add to the observed velocity (km/s), and its parts
        double lsr_corr = 0, lsr_rot = 0, lsr_orb = 0, lsr_sun = 0;
        bool lsr_qth_ok = false;
        double qth_lat = 0, qth_lon = 0, qth_alt = 0;
        // Updates lsr_* for the current time (call from the UI thread)
        void update_lsr();
        static void galactic_to_j2000(double l, double b, double &ra_deg, double &dec_deg);
        static void lsr_correction(double jd, double lat, double lon, double alt_m, double ra_deg, double dec_deg, double &rot, double &orb, double &sun);

        // Narrow interference (spurs, e.g. the Airspy clock harmonic at
        // 1420.000 MHz): in the mean spectra A and B, bins standing above the
        // local median by more than rfi_sigma times the noise are replaced
        // by an interpolation of their neighbours. The local median is taken
        // over blocks of rfi_width_khz, much narrower than the HI line but
        // much wider than a spur. The raw sums are kept: unticking restores.
        bool rfi_remove = true;
        float rfi_sigma = 6.0f;
        float rfi_width_khz = 10.0f;
        int rfi_grow = 2; // also flag this many bins on each side

        // Tunes the SDR without changing the recorder frequency
        std::function<void(double)> set_sdr_frequency = [](double) {};

        // With a remote SDR the server switches (start_remote): it knows
        // which samples come after each retune, and only one mean spectrum
        // per dwell crosses the network. A local SDR is switched here; its
        // samples arrive within Settle after a retune.
        bool save_iq_on_server = false;
        int iq_bits = 16; // 8, 16 or 32 (float)
        // Sky frequency -> SDR frequency (xconverter), without tuning
        std::function<double(double)> sdr_freq_of = [](double f) { return f; };
        // Asks the source to switch by itself; false if it cannot
        std::function<bool(nlohmann::json, std::shared_ptr<dsp::remote_fswitch::Callbacks>)> remote_start =
            [](nlohmann::json, std::shared_ptr<dsp::remote_fswitch::Callbacks>) { return false; };
        std::function<void()> remote_stop = []() {};

        // Main display while switching: each frame is moved back to the sky
        // frequency (by -/+ the offset) then averaged like the normal FFT
        float display_avg = 10;
        std::function<void(float *)> on_display = [](float *) {};
        float *display_buffer() { return disp.data(); }

    public:
        // In server mode, do not talk to the source here: it may be gone
        ~FrequencySwitcher()
        {
            if (remote)
                running = false;
            else
                stop();
        }

        // frame_dt: duration of signal used by one FFT frame (seconds)
        // init_display: spectrum shown when switching starts (dB, fft_size values), may be null
        // resume: keep the averages of the previous run (see can_resume)
        void start(double center_hz, double samplerate, int fft_size, double frame_dt, const float *init_display = nullptr, bool resume = false);
        // Same, done by the server. False if the source cannot do it
        bool start_remote(double center_hz, double samplerate, int fft_size, const float *init_display = nullptr, bool resume = false);
        // True when stopped with averages that a new run with these settings
        // can continue: same frequency, samplerate, FFT size and offset
        bool can_resume(double center_hz, double samplerate, int fft_size, double offset_khz);
        void stop();
        bool is_remote() { return remote; }
        // True once when the server refused or did not answer: the recorder
        // then switches locally. Call from the UI thread
        bool take_fallback(std::string &why);
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
        // File of the folded cycles (empty if none)
        std::string cycles_file()
        {
            std::lock_guard<std::mutex> lk(mtx);
            return cyc_path;
        }

    private:
        // Server mode
        bool remote = false;
        // Callbacks hold a weak copy: they do nothing once this object is
        // gone, and a late packet of an earlier run is ignored (generation)
        std::shared_ptr<int> alive_token = std::make_shared<int>(0);
        std::atomic<uint64_t> remote_gen{0};
        bool remote_status_seen = false;
        double remote_t_start = 0;
        std::atomic<bool> fallback_requested{false};
        std::string fallback_why;
        std::string remote_iq_file, remote_iq_error;
        uint64_t remote_iq_bytes = 0;
        uint64_t remote_dwells = 0;
        void on_remote_status(nlohmann::json st);
        void on_remote_dwell(int phase, int frames, uint64_t index, const float *mean, int n);
        // Returns true if the averages of the previous run were kept
        bool init_run(double center_hz, double sr, int fft_size, const float *init_display, bool resume);
        int dk_for(double offset_khz, double sr, int fft_size);
        bool compatible(double center_hz, double sr, int fft_size, int new_dk);

        void run();
        void compute();
        int phase_of(double fstart, double fend, double delay);

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
        double run_dwell = 0, run_settle = 0;

        // Retune events: time (s) and phase (0 = f0 + D (A), 1 = f0 - D (B), -1 = back to f0)
        std::deque<std::pair<double, int>> events;

        // Frames are dated by their arrival time: samples lost on the way
        // do not shift the dates.
        double frame_dt = 0;
        double t_first = 0;
        uint64_t n_frames = 0;
        double lag_s = 0;
        double speed = 1;
        // Sliding measurement over the last seconds: (arrival time, frame count)
        std::deque<std::pair<double, uint64_t>> rate_hist;
        double speed_now = 1;

        // Accumulators (linear power)
        std::vector<double> sum_a, sum_b;
        std::vector<double> lin_tmp;
        uint64_t n_a = 0, n_b = 0;
        uint64_t cycles = 0;
        // Their values when this run started (0, or those of a resumed run)
        uint64_t n_a_run0 = 0, n_b_run0 = 0, cycles_run0 = 0;
        bool resumed = false;

        // Per cycle folding, saved to cyc_file
        std::ofstream cyc_file;
        std::string cyc_path, cyc_error;
        std::vector<double> cyc_a, cyc_b;
        uint64_t cyc_na = 0, cyc_nb = 0, cyc_count = 0;
        int cyc_last_ph = -1;
        // append: continue the file of the resumed run
        void open_cycles(bool append = false);
        void close_cycles();
        // Adds `weight` frames of linear power (mean) to the current cycle;
        // a new A after a B closes the previous cycle
        void cycle_add(int ph, const float *lin_f, const double *lin_d, double weight);
        void write_cycle();

        // Computed for display
        std::vector<double> plot_x, plot_y;
        // Unfolded: o1 - o2 = (A-B)/B and o2 - o1 = (B-A)/A, receiver axis
        std::vector<double> unf_x, unf_12, unf_21;
        void smooth_scale(const std::vector<double> &y, std::vector<double> &out);
        // Fits and subtracts the baseline from y (x in MHz, NaN = hidden);
        // the subtracted curve is kept in bl_curve
        void remove_baseline(const std::vector<double> &x, std::vector<double> &y);
        std::vector<double> bl_curve, folded_raw;

        // Interference flags of A and B (receiver axis), from the whole
        // averages, used for the result and the cycles file too
        std::vector<uint8_t> rfi_a, rfi_b;
        int rfi_count_a = 0, rfi_count_b = 0;
        uint64_t rfi_na = 0, rfi_nb = 0; // n_a, n_b when last detected
        float rfi_last_sigma = 0, rfi_last_width = 0;
        void rfi_detect(const std::vector<double> &m, std::vector<uint8_t> &flag, int &count);
        static void rfi_fill(std::vector<double> &m, const std::vector<uint8_t> &flag);
        // Mean A and B, with the interference removed if enabled
        void mean_spectra(std::vector<double> &a, std::vector<double> &b);
        int computed_frame = -1; // compute() once per UI frame
        std::atomic<bool> stop_requested{false};
        void compute_once();
        void plot_folded(ImVec2 size, bool main_area);
        // Velocity axis: frequency range shown last frame, ticks and labels
        double vel_lim[2] = {0, 0};
        std::vector<double> vel_ticks;
        std::vector<std::string> vel_text;
        std::vector<const char *> vel_labels;
        void velocity_ticks(double fmin_mhz, double fmax_mhz);

        // Main display spectrum (dB), on the sky frequency axis
        std::vector<float> disp;
        // Recentred averages of each position, linear power (0 = no data)
        std::vector<float> disp_a, disp_b;
        void update_display(const float *db, int ph);
        // Same with a linear spectrum, mixed with weight r
        void update_display_lin(const float *lin, int ph, float r);
    };
} // namespace satdump
