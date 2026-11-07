#pragma once

/*
 * Frequency switching done by a remote SDR server (experimental).
 *
 * The recorder fires these events with its source. A source that can
 * switch by itself (the remote source plugin) sets *accepted and then
 * calls the callbacks from its network thread. Events are used rather than
 * new virtual functions, so DSPSampleSource and the plugins built against
 * it stay unchanged.
 */

#include "dsp_sample_source.h"
#include "nlohmann/json.hpp"
#include <functional>
#include <memory>

namespace dsp
{
    namespace remote_fswitch
    {
        struct Callbacks
        {
            // Status from the server (CBOR decoded): running, error, iq_file...
            std::function<void(nlohmann::json)> on_status = [](nlohmann::json) {};
            // Mean linear power spectrum of one dwell (FFT shifted)
            std::function<void(int phase, int frames, uint64_t dwell_index, const float *mean, int nbins)> on_dwell =
                [](int, int, uint64_t, const float *, int) {};
        };

        // settings: freq_a, freq_b, freq_back (SDR Hz), sky_f0, offset_hz,
        // fft_size, dwell_s, settle_s, save_iq, iq_bits
        struct StartEvent
        {
            DSPSampleSource *source;
            nlohmann::json settings;
            std::shared_ptr<Callbacks> callbacks;
            bool *accepted;
        };

        struct StopEvent
        {
            DSPSampleSource *source;
        };
    } // namespace remote_fswitch
} // namespace dsp
