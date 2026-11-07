#include "streaming.h"
#include "main.h"
#include "remote.h"
#include "iq_pkt.h"
#include "fswitch.h"

std::mutex source_stream_mtx;
bool source_should_stream = false;
std::atomic<int> streaming_bit_depth = 32;

void sourceStreamThread()
{
    uint8_t *buffer_tx = new uint8_t[dsp::STREAM_BUFFER_SIZE * sizeof(complex_t)];
    float *mag_buffer = new float[dsp::STREAM_BUFFER_SIZE];
    while (1)
    {
        int swrite_ret = 1;
        source_stream_mtx.lock();
        if (source_should_stream)
        {
            int nsamples = current_sample_source->output_stream->read();

            if (nsamples <= 0)
            {
                source_stream_mtx.unlock();
                continue;
            }

            if (server_fswitch.active())
            {
                // Frequency switching: the samples stay on the server, only
                // the mean spectrum of each dwell is sent
                server_fswitch.process(current_sample_source->output_stream->readBuf, nsamples);
            }
            else
            {
                int pktlen = 1;
                buffer_tx[0] = dsp::remote::PKT_TYPE_IQ;

                pktlen += remote_sdr::encode_iq_pkt(&buffer_tx[1], current_sample_source->output_stream->readBuf, mag_buffer, nsamples, streaming_bit_depth);
                swrite_ret = tcp_server->swrite(buffer_tx, pktlen);
            }

            // logger->trace(nsamples);
            current_sample_source->output_stream->flush();
        }
        source_stream_mtx.unlock();

        if(swrite_ret <= 0)
            tcp_server->closeconn();

        if (!source_should_stream)
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    delete[] buffer_tx;
    delete[] mag_buffer;
}

// Retunes asked by the frequency switching, in their own thread: the stream
// thread must keep reading, or drivers whose callback blocks on the stream
// (librtlsdr) cannot complete the retune. The switcher drops the samples read
// until retune_done(). The retune is taken with source_mtx held, so a
// frequency set by the client (which cancels it under the same lock) wins.
void sourceRetuneThread()
{
    while (1)
    {
        server_fswitch.wait_retune(100);
        if (!server_fswitch.retune_wanted())
            continue;
        double freq;
        bool done = false;
        {
            std::lock_guard<std::mutex> lock(source_mtx);
            if (server_fswitch.take_retune(freq))
            {
                if (source_is_open && current_sample_source)
                    current_sample_source->set_frequency(freq);
                done = true;
            }
        }
        if (done)
            server_fswitch.retune_done();
    }
}
