#pragma once

/*
 * Frequency switching done by the SDR server (experimental).
 *
 * The server knows exactly which samples were taken before and after each
 * retune, so it does the switching itself: no delay has to be measured on
 * the client, and the IQ stream does not cross the network while switching.
 *
 * For each dwell the server sends the mean power spectrum (linear, the same
 * scale as the FFT of the recorder) and the number of FFT frames it holds.
 * The client sums them into A and B and does the folding.
 *
 * FSWITCH_DWELL payload, little endian:
 *   int32  phase       0 = tuned for f0 + D (A), 1 = tuned for f0 - D (B)
 *   int32  frames      FFT frames averaged in this dwell
 *   uint64 dwell_index counted from the start, both phases
 *   int32  nbins
 *   float  mean[nbins] mean linear power per bin, FFT shifted (DC in the middle)
 */

#include <cstdint>
#include <cstring>
#include <vector>

namespace remote_sdr
{
    struct FSwitchDwell
    {
        int32_t phase = 0;
        int32_t frames = 0;
        uint64_t dwell_index = 0;
        std::vector<float> mean;
    };

    inline std::vector<uint8_t> encode_fswitch_dwell(const FSwitchDwell &d)
    {
        int32_t nbins = d.mean.size();
        std::vector<uint8_t> out(4 + 4 + 8 + 4 + nbins * sizeof(float));
        uint8_t *p = out.data();
        memcpy(p, &d.phase, 4);
        memcpy(p + 4, &d.frames, 4);
        memcpy(p + 8, &d.dwell_index, 8);
        memcpy(p + 16, &nbins, 4);
        if (nbins > 0)
            memcpy(p + 20, d.mean.data(), nbins * sizeof(float));
        return out;
    }

    // Returns false if the payload is malformed
    inline bool decode_fswitch_dwell(const uint8_t *in, int len, FSwitchDwell &d)
    {
        if (len < 20)
            return false;
        int32_t nbins;
        memcpy(&d.phase, in, 4);
        memcpy(&d.frames, in + 4, 4);
        memcpy(&d.dwell_index, in + 8, 8);
        memcpy(&nbins, in + 16, 4);
        if (nbins < 0 || len < 20 + (int64_t)nbins * (int)sizeof(float))
            return false;
        d.mean.resize(nbins);
        if (nbins > 0)
            memcpy(d.mean.data(), in + 20, nbins * sizeof(float));
        return true;
    }
} // namespace remote_sdr
