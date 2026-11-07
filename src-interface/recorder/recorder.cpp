#include "recorder.h"
#include "common/utils.h"
#include "common/widgets/frequency_input.h"
#include "common/widgets/stepped_slider.h"
#include "core/config.h"
#include "core/plugin.h"
#include "core/style.h"
#include "handlers/handler.h"
#include "imgui/imgui.h"
#include "imgui/imgui_stdlib.h"
#include "logger.h"
#include <math.h>

#include "main_ui.h"

#include "imgui/imgui_internal.h"

#include "core/resources.h"
#include "pipeline/pipeline.h"

namespace satdump
{
    RecorderApplication::RecorderApplication(nlohmann::json cli_set) : handlers::Handler(), pipeline_selector(true)
    {
        automated_live_output_dir = satdump_cfg.main_cfg["satdump_directories"]["live_processing_autogen"]["value"].get<bool>();
        processing_modules_floating_windows = satdump_cfg.main_cfg["user_interface"]["recorder_floating_windows"]["value"].get<bool>();
        remaining_disk_space_time = satdump_cfg.main_cfg["user_interface"]["remaining_disk_space_time"]["value"].get<int>();

        load_rec_path_data();
        dsp::registerAllSources();
        sources = dsp::getAllAvailableSources();

        for (dsp::SourceDescriptor src : sources)
        {
            logger->debug("Device " + src.name);
            sdr_select_string += src.name + '\0';
        }

        // First, check to see if a source was passed via command line
        if (satdump::satdump_cfg.main_cfg.contains("cli") && satdump::satdump_cfg.main_cfg["cli"].contains("source"))
        {
            auto &cli_settings = satdump::satdump_cfg.main_cfg["cli"];
            std::string source = cli_settings["source"];
            for (int i = 0; i < (int)sources.size(); i++)
            {
                if (sources[i].source_type == source)
                {
                    if (cli_settings.contains("source_id") && cli_settings["source_id"].get<std::string>() != sources[i].unique_id)
                        continue;

                    try
                    {
                        source_ptr = dsp::getSourceFromDescriptor(sources[i]);
                        source_ptr->open();
                        sdr_select_id = i;
                        break;
                    }
                    catch (std::runtime_error &e)
                    {
                        logger->error(e.what());
                    }
                }
            }
        }

        // If not, or it failed to open, revert to last used SDR
        auto dcfg = db->get_user_json("recorder_sdr_settings");
        if (sdr_select_id == -1 && dcfg.contains("last_used_sdr"))
        {
            std::string last_used_sdr = dcfg["last_used_sdr"].get<std::string>();
            for (int i = 0; i < (int)sources.size(); i++)
            {
                if (sources[i].name != last_used_sdr)
                    continue;
                try
                {
                    source_ptr = dsp::getSourceFromDescriptor(sources[i]);
                    source_ptr->open();
                    sdr_select_id = i;
                    break;
                }
                catch (std::runtime_error &e)
                {
                    logger->error(e.what());
                    break;
                }
            }
        }

        // If no source has been opened yet, open the first possible one
        for (int i = 0; sdr_select_id == -1 && i < (int)sources.size(); i++)
        {
            try
            {
                source_ptr = dsp::getSourceFromDescriptor(sources[i]);
                source_ptr->open();
                sdr_select_id = i;
            }
            catch (std::runtime_error &e)
            {
                logger->error(e.what());
            }
        }

        // Load waterfall palettes
        {
            waterfall_palettes_str = "";
            std::filesystem::recursive_directory_iterator paletteIterator(resources::getResourcePath("waterfall"));
            std::error_code iteratorError;
            while (paletteIterator != std::filesystem::recursive_directory_iterator())
            {
                if (!std::filesystem::is_directory(paletteIterator->path()))
                {
                    auto map = colormaps::loadMap(paletteIterator->path().string());
                    waterfall_palettes.push_back(map);
                    waterfall_palettes_str += map.name + " [" + map.author + "]" + '\0';
                }

                paletteIterator.increment(iteratorError);
                if (iteratorError)
                    logger->critical(iteratorError.message());
            }
        }

        set_frequency(frequency_hz);
        try_load_sdr_settings();

        splitter = std::make_shared<dsp::SplitterBlock>(source_ptr->output_stream);
        splitter->add_output("record");
        splitter->add_output("live");
        splitter->add_output("fswitch");

        fft = std::make_shared<dsp::FFTPanBlock>(splitter->output_stream);
        fft->set_fft_settings(fft_size, get_samplerate(), fft_rate);
        fft->avg_num = 10;
        fft->start();

        file_sink = std::make_shared<dsp::FileSinkBlock>(splitter->get_output("record"));
        file_sink->start();

        fft_plot = std::make_shared<widgets::FFTPlot>(fft->output_stream->writeBuf, fft_size, -150, 150, 10);
        fft_plot->frequency = frequency_hz;
        waterfall_plot = std::make_shared<widgets::WaterfallPlot>(fft_sizes_lut[0], 500);

        fft->on_fft = [this](float *v)
        {
            // While switching, the waterfall is fed with recentred frames
            if (!fswitch.is_running())
                waterfall_plot->push_fft(v);
        };
        fswitch.on_display = [this](float *v) { waterfall_plot->push_fft(v); };
        fswitch.set_sdr_frequency = [this](double f) { set_sdr_frequency(f); };
        fswitch.sdr_freq_of = [this](double f) { return sdr_freq_of(f); };
        fswitch.remote_start = [this](nlohmann::json settings, std::shared_ptr<dsp::remote_fswitch::Callbacks> cb)
        {
            bool accepted = false;
            if (source_ptr)
                eventBus->fire_event<dsp::remote_fswitch::StartEvent>({source_ptr.get(), settings, cb, &accepted});
            return accepted;
        };
        fswitch.remote_stop = [this]()
        {
            if (source_ptr)
                eventBus->fire_event<dsp::remote_fswitch::StopEvent>({source_ptr.get()});
        };

        // Load config
        auto cfg = db->get_user_json("recorder_state");
        if (!cfg.is_null())
            deserialize_config(cfg);

        file_sink->set_output_sample_type(baseband_format);
        fft_plot->set_size(fft_size);
        waterfall_plot->set_size(fft_size);
        waterfall_plot->set_rate(fft_rate, waterfall_rate);

        // Attempt to apply provided CLI settings
        if (!cli_set.is_null())
        {
            auto &cli_settings = cli_set;
            if (source_ptr)
            {
                if (cli_settings.contains("samplerate"))
                    source_ptr->set_samplerate(cli_settings["samplerate"].get<uint64_t>());
                source_ptr->set_settings(cli_settings);
            }

            if (cli_settings.contains("start_recorder_device") && cli_settings["start_recorder_device"].get<bool>())
            {
                logger->warn("Recorder was asked to autostart!");
                start();
            }

            if (cli_settings.contains("engage_autotrack") && cli_settings["engage_autotrack"].get<bool>())
            {
                tracking_started_cli = true;
                try_init_tracking_widget(cli_set);
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        eventBus->register_handler<RecorderSetFrequencyEvent>([this](const RecorderSetFrequencyEvent &evt) { set_frequency(evt.frequency); });

        eventBus->register_handler<RecorderStartDeviceEvent>([this](const RecorderStartDeviceEvent &) { start(); });
        eventBus->register_handler<RecorderStopDeviceEvent>([this](const RecorderStopDeviceEvent &) { stop(); });
        eventBus->register_handler<RecorderSetDeviceSamplerateEvent>([this](const RecorderSetDeviceSamplerateEvent &evt) { source_ptr->set_samplerate(evt.samplerate); });
        eventBus->register_handler<RecorderSetDeviceDecimationEvent>([this](const RecorderSetDeviceDecimationEvent &evt) { current_decimation = evt.decimation; });
        eventBus->register_handler<RecorderSetDeviceLoOffsetEvent>([this](const RecorderSetDeviceLoOffsetEvent &evt) { xconverter_frequency = evt.offset; });

        eventBus->register_handler<RecorderStartProcessingEvent>(
            [this](const RecorderStartProcessingEvent &evt)
            {
                pipeline_selector.select_pipeline(evt.pipeline_id);
                start_processing();
            });
        eventBus->register_handler<RecorderStopProcessingEvent>([this](const RecorderStopProcessingEvent &) { stop_processing(); });

        eventBus->register_handler<RecorderSetFFTSettingsEvent>(
            [this](const RecorderSetFFTSettingsEvent &evt)
            {
                if (evt.fft_min != -1)
                    fft_plot->scale_min = waterfall_plot->scale_min = evt.fft_min;
                if (evt.fft_max != -1)
                    fft_plot->scale_max = waterfall_plot->scale_max = evt.fft_max;
                if (evt.fft_avgn != -1)
                    fft->avg_num = evt.fft_avgn;
                if (evt.fft_size != -1 || evt.fft_rate != -1 || evt.waterfall_rate != -1)
                {
                    if (evt.fft_size != -1)
                        fft_size = evt.fft_size;
                    if (evt.fft_rate != -1)
                        fft_rate = evt.fft_rate;
                    if (evt.waterfall_rate != -1)
                        waterfall_rate = evt.waterfall_rate;
                    fft->set_fft_settings(fft_size, get_samplerate(), fft_rate);
                    waterfall_plot->set_rate(fft_rate, waterfall_rate);
                    waterfall_plot->set_size(fft_size);
                    fft_plot->set_size(fft_size);
                    for (int i = 0; i < (int)fft_sizes_lut.size(); i++)
                        if (fft_sizes_lut[i] == fft_size)
                            selected_fft_size = i;
                }
            });
    }

    RecorderApplication::~RecorderApplication()
    {
        save_settings();

    retry_vfo:
        for (auto &vfo : vfo_list)
        {
            del_vfo(vfo.id);
            goto retry_vfo;
        }

        stop_processing();
        if (is_started)
            stop();
        source_ptr->close();
        splitter->input_stream = std::make_shared<dsp::stream<complex_t>>();
        splitter->stop();
        fft->stop();
        file_sink->stop();

        if (tracking_widget != nullptr)
            delete tracking_widget;
        if (constellation_debug != nullptr)
            delete constellation_debug;
    }

    void RecorderApplication::drawMenu()
    {
        bool assume_started = is_started;
        if (ImGui::CollapsingHeader(_("Source"), tracking_started_cli ? ImGuiTreeNodeFlags_None : ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Spacing();
            if (assume_started)
                style::beginDisabled();
            if (ImGui::Combo("##Source", &sdr_select_id, sdr_select_string.c_str()))
            {
                // Try to open a device, if it doesn't work, we re-open a device we can
                try
                {
                    source_ptr.reset();
                    source_ptr = getSourceFromDescriptor(sources[sdr_select_id]);
                    source_ptr->open();
                }
                catch (std::runtime_error &e)
                {
                    logger->error("Could not open device! %s", e.what());

                    for (int i = 0; i < (int)sources.size(); i++)
                    {
                        try
                        {
                            source_ptr.reset();
                            source_ptr = dsp::getSourceFromDescriptor(sources[i]);
                            source_ptr->open();
                            sdr_select_id = i;
                            break;
                        }
                        catch (std::runtime_error &e)
                        {
                            logger->error(e.what());
                        }
                    }
                }

                set_frequency(frequency_hz);
                try_load_sdr_settings();
            }
            ImGui::SameLine();
            if (ImGui::Button(u8" \uead2 "))
            {
                sources = dsp::getAllAvailableSources();

                sdr_select_string.clear();
                for (dsp::SourceDescriptor src : sources)
                {
                    logger->debug("Device " + src.name);
                    sdr_select_string += src.name + '\0';
                }

                while (sdr_select_id >= (int)sources.size())
                    sdr_select_id--;

                source_ptr = getSourceFromDescriptor(sources[sdr_select_id]);
                source_ptr->open();
                set_frequency(frequency_hz);
                try_load_sdr_settings();
            }
            /*
            if (ImGui::IsItemHovered())
                {
                    ImGui::BeginTooltip();
                    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
                    ImGui::TextUnformatted("Refresh source list");
                    ImGui::PopTextWrapPos();
                    ImGui::EndTooltip();
                }
            */

            ImGui::InputInt(_("Decimation##recorderdecimation"), &current_decimation);
            if (current_decimation < 1)
                current_decimation = 1;

            bool disableLO = satdump_cfg.main_cfg["user_interface"]["lock_xconverter_input"]["value"];
            if (assume_started && !disableLO)
                style::endDisabled();

            bool pushed_color_xconv = xconverter_frequency != 0;
            if (pushed_color_xconv)
                ImGui::PushStyleColor(ImGuiCol_Text, style::theme.green.Value);

            if (ImGui::InputDouble(_("MHz (LO offset)##downupconverter"), &xconverter_frequency))
                set_frequency(frequency_hz);

            if (pushed_color_xconv)
                ImGui::PopStyleColor();

            if (assume_started && disableLO)
                style::endDisabled();

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            if (widgets::FrequencyInput("Hz##mainfreq", &frequency_hz))
                set_frequency(frequency_hz);

            ImGui::Spacing();
            source_ptr->drawControlUI();

            if (!assume_started)
            {
                if (ImGui::Button(_("Start")))
                    start();
            }
            else
            {
                if (ImGui::Button(_("Stop")))
                    stop();
            }

            sdr_error.draw();
        }

        if (ImGui::CollapsingHeader(_("FFT"), tracking_started_cli ? ImGuiTreeNodeFlags_None : ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Combo(_("FFT Size"), &selected_fft_size,
                             "131072\0"
                             "65536\0"
                             "32768\0"
                             "16384\0"
                             "8192\0"
                             "4096\0"
                             "2048\0"
                             "1024\0"))
            {
                if (fswitch.is_running())
                    stop_fswitch();
                fft_size = fft_sizes_lut[selected_fft_size];

                fft->set_fft_settings(fft_size, get_samplerate(), fft_rate);
                fft_plot->set_size(fft_size);
                waterfall_plot->set_size(fft_size);

                logger->info("Set FFT size to %d", fft_size);
            }
            int old_rate = fft_rate;
            if (ImGui::InputInt(_("FFT Rate"), &fft_rate))
            {
                if (fft_rate <= 0)
                    fft_rate = old_rate;
                else
                {
                    fft_size = fft_sizes_lut[selected_fft_size];
                    fft->set_fft_settings(fft_size, get_samplerate(), fft_rate);
                    waterfall_plot->set_rate(fft_rate, waterfall_rate);
                    logger->info("Set FFT rate to %d", fft_rate);
                }
            }
            if (ImGui::IsItemDeactivatedAfterEdit() && fft_rate < waterfall_rate)
            {
                waterfall_rate = fft_rate;
                waterfall_plot->set_rate(fft_rate, waterfall_rate);
                logger->info("Set Waterfall rate to %d", waterfall_rate);
            }
            old_rate = waterfall_rate;
            if (ImGui::InputInt(_("Waterfall Rate"), &waterfall_rate))
            {
                if (waterfall_rate <= 0)
                    waterfall_rate = old_rate;
                else
                {
                    waterfall_plot->set_rate(fft_rate, waterfall_rate);
                    logger->info("Set Waterfall rate to %d", waterfall_rate);
                }
            }
            if (ImGui::IsItemDeactivatedAfterEdit() && waterfall_rate > fft_rate)
            {
                fft_rate = waterfall_rate;
                fft->set_fft_settings(fft_size, get_samplerate(), fft_rate);
                waterfall_plot->set_rate(fft_rate, waterfall_rate);
                logger->info("Set FFT rate to %d", fft_rate);
            }
            widgets::SteppedSliderFloat(_("FFT Max"), &fft_plot->scale_max, -160, 150);
            widgets::SteppedSliderFloat(_("FFT Min"), &fft_plot->scale_min, -160, 150);
            widgets::SteppedSliderFloat(_("Avg Num"), &fft->avg_num, 1, 500, 1);
            fswitch.display_avg = fft->avg_num;
            if (ImGui::Combo(_("Palette"), &selected_waterfall_palette, waterfall_palettes_str.c_str()))
                waterfall_plot->set_palette(waterfall_palettes[selected_waterfall_palette]);
            ImGui::Checkbox(_("Show Waterfall"), &show_waterfall);
            ImGui::Checkbox(_("Frequency Scale"), &fft_plot->enable_freq_scale);
        }

        if (fft_plot->scale_max < fft_plot->scale_min)
        {
            fft_plot->scale_min = waterfall_plot->scale_min;
            fft_plot->scale_max = waterfall_plot->scale_max;
        }
        else if (fft_plot->scale_min > fft_plot->scale_max)
        {
            fft_plot->scale_min = waterfall_plot->scale_min;
            fft_plot->scale_max = waterfall_plot->scale_max;
        }
        else
        {
            waterfall_plot->scale_min = fft_plot->scale_min;
            waterfall_plot->scale_max = fft_plot->scale_max;
        }

        bool fswitch_open = ImGui::CollapsingHeader(_("Frequency switching (experimental)"));
        if (fswitch_open)
        {
            ImGui::TextWrapped("%s", _("For very stable signals only (radio astronomy), not for satellites. Recording, processing and tracking must be off."));
            bool running = fswitch.is_running();
            if (running)
                ImGui::BeginDisabled();
            {
                // Offset in % of the largest one allowed, shown in kHz too
                double max_khz = FrequencySwitcher::max_offset_hz(get_samplerate()) / 1e3;
                char fmt[64];
                snprintf(fmt, sizeof(fmt), "%%.0f %%%% (%.0f kHz)", fswitch.offset_pct / 100.0 * max_khz);
                widgets::SteppedSliderFloat(_("Offset (% of max)"), &fswitch.offset_pct, 5, 100, 1, fmt);
                ImGui::TextDisabled("%s %.0f kHz", _("Max offset at this samplerate:"), max_khz);
            }
            ImGui::InputFloat(_("Dwell (s)"), &fswitch.dwell_s, 0.1f, 1.0f, "%.2f");
            ImGui::InputFloat(_("Settle (s)"), &fswitch.settle_s, 0.05f, 0.5f, "%.2f");
            // Remote SDR: the server switches by itself (exact timing, no IQ on
            // the network)
            if (source_is_remote())
            {
                ImGui::Indent();
                ImGui::Checkbox(_("Save IQ on the server"), &fswitch.save_iq_on_server);
                if (fswitch.save_iq_on_server)
                {
                    int bits_sel = fswitch.iq_bits == 8 ? 0 : fswitch.iq_bits == 32 ? 2 : 1;
                    if (ImGui::Combo(_("IQ format##fswitch"), &bits_sel, "8 bits (cs8)\0" "16 bits (cs16)\0" "32 bits float (cf32)\0"))
                        fswitch.iq_bits = bits_sel == 0 ? 8 : bits_sel == 2 ? 32 : 16;
                    double gb_h = get_samplerate() * 2.0 * (fswitch.iq_bits / 8) * 3600 / 1e9;
                    ImGui::TextDisabled("%s %.0f GB/h", _("About"), gb_h);
                }
                ImGui::Unindent();
            }
            ImGui::Checkbox(_("Save each folded cycle (CSV, recording folder)"), &fswitch.save_cycles);
            if (running)
                ImGui::EndDisabled();

            if (!running)
            {
                bool can_start = is_started && !is_recording && !is_processing && vfo_list.size() == 0;
                if (!can_start)
                    ImGui::BeginDisabled();
                if (ImGui::Button(_("Start switching")))
                    start_fswitch();
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                    ImGui::SetTooltip("%s", _("New run: the averages start from zero"));
                // Continue the averages where Stop left them
                uint64_t sr_now = get_samplerate();
                double off_khz = fswitch.offset_pct / 100.0 * FrequencySwitcher::max_offset_hz(sr_now) / 1e3;
                bool can_resume = fswitch.can_resume(frequency_hz, sr_now, fft_size, off_khz);
                ImGui::SameLine();
                if (!can_resume)
                    ImGui::BeginDisabled();
                if (ImGui::Button(_("Resume switching")))
                    start_fswitch(true);
                if (!can_resume)
                    ImGui::EndDisabled();
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                    ImGui::SetTooltip("%s", can_resume ? _("Continue the averages where Stop left them")
                                                       : _("Nothing to continue, or frequency, samplerate, FFT size or offset changed"));
                if (!can_start)
                    ImGui::EndDisabled();
            }
            else if (ImGui::Button(_("Stop switching")))
                stop_fswitch();
            ImGui::SameLine();
            if (ImGui::Button(_("Reset##fswitch")))
            {
                // Full restart: the delay is measured again from a clean pattern
                if (fswitch.is_running())
                {
                    stop_fswitch();
                    start_fswitch();
                }
                else
                    fswitch.reset();
            }
            ImGui::SameLine();
            if (ImGui::Button(_("Save CSV")))
            {
                load_rec_path_data();
                std::string file = fswitch.save_csv(recording_path);
                if (file.size())
                    logger->info("Frequency switching result saved to " + file);
                else
                    logger->error("Frequency switching: nothing to save, or could not write in " + recording_path);
            }

            ImGui::SliderInt(_("Smoothing (bins)"), &fswitch.smoothing, 1, 256);
            if (ImGui::RadioButton(_("Linear##fswitch"), !fswitch.use_db))
                fswitch.use_db = false;
            ImGui::SameLine();
            if (ImGui::RadioButton(_("dB##fswitch"), fswitch.use_db))
                fswitch.use_db = true;
            ImGui::Checkbox(_("Auto Y scale"), &fswitch.auto_y);
            if (fswitch.auto_y)
                ImGui::BeginDisabled();
            if (fswitch.use_db)
            {
                widgets::SteppedSliderFloat(_("Y Max (dB)"), &fswitch.y_max_db, -20, 20, 0.001, "%.3f", ImGuiSliderFlags_Logarithmic);
                widgets::SteppedSliderFloat(_("Y Min (dB)"), &fswitch.y_min_db, -20, 20, 0.001, "%.3f", ImGuiSliderFlags_Logarithmic);
            }
            else
            {
                widgets::SteppedSliderFloat(_("Y Max (x1e-3)"), &fswitch.y_max, -1000, 1000, 0.01, "%.2f", ImGuiSliderFlags_Logarithmic);
                widgets::SteppedSliderFloat(_("Y Min (x1e-3)"), &fswitch.y_min, -1000, 1000, 0.01, "%.2f", ImGuiSliderFlags_Logarithmic);
            }
            if (fswitch.auto_y)
                ImGui::EndDisabled();
            widgets::SteppedSliderFloat(_("Edge mask (dB)"), &fswitch.edge_db, 0.5, 20, 0.5);
            ImGui::Checkbox(_("Velocity scale (km/s)"), &fswitch.show_velocity);
            if (fswitch.show_velocity)
            {
                ImGui::Indent();
                ImGui::InputDouble(_("Rest frequency (MHz)"), &fswitch.rest_mhz, 0, 0, "%.6f");
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s", _("HI 1420.405752, OH 1665.402 / 1667.359. v = c (f_rest - f) / f_rest, positive = receding"));
                ImGui::Checkbox(_("LSR correction"), &fswitch.lsr);
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s", _("Velocity relative to the local standard of rest: removes the Earth rotation, the Earth orbit and the Sun motion.\nNeeds the source position and the station position (SatDump settings)."));
                if (fswitch.lsr)
                {
                    ImGui::Combo(_("Coordinates##fswitchlsr"), &fswitch.lsr_frame, "RA/Dec J2000\0Galactic l/b\0");
                    if (fswitch.lsr_frame == 0)
                    {
                        ImGui::InputDouble(_("RA (hours)"), &fswitch.src_ra_h, 0, 0, "%.4f");
                        ImGui::InputDouble(_("Dec (degrees)"), &fswitch.src_dec, 0, 0, "%.4f");
                    }
                    else
                    {
                        ImGui::InputDouble(_("l (degrees)"), &fswitch.src_l, 0, 0, "%.3f");
                        ImGui::InputDouble(_("b (degrees)"), &fswitch.src_b, 0, 0, "%.3f");
                    }
                    fswitch.update_lsr();
                    if (fswitch.lsr_qth_ok)
                        ImGui::TextDisabled("Correction %+.2f km/s (rotation %+.2f, orbit %+.2f, Sun %+.2f)", fswitch.lsr_corr, fswitch.lsr_rot, fswitch.lsr_orb, fswitch.lsr_sun);
                    else
                        ImGui::TextColored(ImVec4(1, 0.3, 0.3, 1), "%s", _("Station position not set (SatDump settings)"));
                }
                ImGui::Unindent();
            }
            ImGui::Checkbox(_("Remove baseline"), &fswitch.bl_remove);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", _("Fits a polynomial outside the line window (orange lines, draggable) and subtracts it,\nto flatten the slow ripple left by the receiver."));
            if (fswitch.bl_remove)
            {
                ImGui::Indent();
                ImGui::Checkbox(_("Show line window"), &fswitch.bl_show_lines);
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s", _("The two orange lines; a right click on the plot also shows or hides them"));
                ImGui::SliderInt(_("Order"), &fswitch.bl_order, 0, 9);
                ImGui::InputFloat(_("Line from (km/s)"), &fswitch.bl_vmin, 10, 50, "%.0f");
                ImGui::InputFloat(_("Line to (km/s)"), &fswitch.bl_vmax, 10, 50, "%.0f");
                ImGui::Unindent();
            }
            ImGui::Checkbox(_("Remove narrow interference"), &fswitch.rfi_remove);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", _("Spurs much narrower than the HI line (e.g. 1420.000 MHz) are removed from A and B before folding.\nThe raw averages are kept: untick to see them again."));
            if (fswitch.rfi_remove)
            {
                ImGui::Indent();
                widgets::SteppedSliderFloat(_("Threshold (x noise)"), &fswitch.rfi_sigma, 3, 50, 0.5, "%.1f");
                widgets::SteppedSliderFloat(_("Local median width (kHz)"), &fswitch.rfi_width_khz, 1, 100, 1, "%.0f");
                ImGui::Unindent();
            }
            ImGui::Checkbox(_("Show folded window after stop"), &fswitch.show_window);
            ImGui::Checkbox(_("Show unfolded window (o1-o2, o2-o1)"), &fswitch.show_unfolded_window);
            if (fswitch_error.size())
                ImGui::TextColored(ImVec4(1, 0.3, 0.3, 1), "%s", fswitch_error.c_str());
            fswitch.draw_status();
        }

        // With Frequency switching open (or running), only Source, FFT and
        // Frequency switching are useful: hide the rest to keep the panel short
        if (fswitch_open || fswitch.is_running())
            show_tracking = false;
        else
        {
            if (ImGui::CollapsingHeader(_("Processing")))
            {
                // Settings & Selection menu
                bool assume_processing = is_processing;
                if (assume_processing)
                    style::beginDisabled();
                pipeline_selector.renderSelectionBox(ImGui::GetContentRegionAvail().x);
                if (!automated_live_output_dir)
                    pipeline_selector.drawMainparamsLive();
                pipeline_selector.renderParamTable();
                if (assume_processing)
                    style::endDisabled();

                // Preset Menu
                if (pipeline_selector.selected_pipeline.preset.frequencies.size() > 0)
                {
                    if (ImGui::BeginCombo(_("Freq###presetscombo"), pipeline_selector.selected_pipeline.preset.frequencies[pipeline_preset_id].second == frequency_hz
                                                                        ? pipeline_selector.selected_pipeline.preset.frequencies[pipeline_preset_id].first.c_str()
                                                                        : ""))
                    {
                        for (int n = 0; n < (int)pipeline_selector.selected_pipeline.preset.frequencies.size(); n++)
                        {
                            const bool is_selected = (pipeline_preset_id == n);
                            if (ImGui::Selectable(pipeline_selector.selected_pipeline.preset.frequencies[n].first.c_str(), is_selected))
                            {
                                pipeline_preset_id = n;

                                if (pipeline_selector.selected_pipeline.preset.frequencies[pipeline_preset_id].second != 0)
                                {
                                    frequency_hz = pipeline_selector.selected_pipeline.preset.frequencies[pipeline_preset_id].second;
                                    set_frequency(frequency_hz);
                                }
                            }

                            if (is_selected)
                                ImGui::SetItemDefaultFocus();
                        }
                        ImGui::EndCombo();
                    }
                }

                if (!assume_started)
                    style::beginDisabled();

                bool assume_stopping_processing = is_stopping_processing;
                if (!assume_processing)
                {
                    if (ImGui::Button(_("Start###startprocessing")))
                        start_processing();
                }
                else if (assume_stopping_processing)
                {
                    style::beginDisabled();
                    ImGui::Button(_("Stopping...##stoppingprocessing"));
                    style::endDisabled();
                }
                else
                {
                    if (ImGui::Button(_("Stop##stopprocessing")))
                        ui_thread_pool.push([=](int) { stop_processing(); });
                }

                error.draw();

                if (!assume_started)
                    style::endDisabled();
            }

            if (ImGui::CollapsingHeader(_("Recording")))
            {
                bool assume_recording = is_recording;
                if (assume_recording)
                    style::beginDisabled();

                if (baseband_format.draw_record_combo(_("Format##basebandrecordformat")))
                    file_sink->set_output_sample_type(baseband_format);

                if (assume_recording)
                    style::endDisabled();

                uint64_t file_written = file_sink->get_written();
                uint64_t estimated_available = 0;
                if (file_written <= disk_available)
                    estimated_available = disk_available - file_written;

                if (file_written < 1e9)
                    ImGui::Text(_("Size : %.2f MB"), file_written / 1e6);
                else
                    ImGui::Text(_("Size : %.2f GB"), file_written / 1e9);

                ImGui::Text(_("Free Space: %.2f GB"), estimated_available / pow(1024, 3));

                int timeleft;
                switch (baseband_format)
                {
                case dsp::CF_32:
                case dsp::CS_32:
                    timeleft = estimated_available / (8 * get_samplerate());
                    break;
                case dsp::CS_16:
                case dsp::WAV_16:
                    timeleft = estimated_available / (4 * get_samplerate());
                    break;
                case dsp::CS_8:
                case dsp::CU_8:
                    timeleft = estimated_available / (2 * get_samplerate());
                    break;
                default:
                    // Silence GCC warns
                    timeleft = 0;
                    break;
                }
    #ifdef BUILD_ZIQ
                if (baseband_format != dsp::ZIQ)
    #endif
                {
    #ifdef BUILD_ZIQ2
                    if (baseband_format != dsp::ZIQ2)
    #endif
                    {
                        if (is_recording && remaining_disk_space_time > timeleft && !been_warned)
                        {
                            logger->warn(_("!!!!WARNING - LOW AMOUNT OF FREE DISK SPACE!!!!"));
                            been_warned = true;
                        }

                        int day = timeleft / (24 * 3600);

                        timeleft = timeleft % (24 * 3600);
                        int hour = timeleft / 3600;

                        timeleft %= 3600;
                        int minutes = timeleft / 60;

                        timeleft %= 60;
                        int seconds = timeleft;
                        ImGui::Text(_("Time left: %02d:%02d:%02d:%02d"), day, hour, minutes, seconds);
                    }
                }

    #ifdef BUILD_ZIQ
                if (baseband_format == dsp::ZIQ)
                {
                    if (file_sink->get_written_raw() < 1e9)
                        ImGui::Text(_("Size (raw) : %.2f MB"), file_sink->get_written_raw() / 1e6);
                    else
                        ImGui::Text(_("Size (raw) : %.2f GB"), file_sink->get_written_raw() / 1e9);
                }
    #endif

                ImGui::Text(_("File : %s"), recorder_filename.c_str());

                ImGui::Spacing();

                if (!assume_recording)
                    ImGui::TextColored(style::theme.red, _("IDLE"));
                else
                    ImGui::TextColored(style::theme.green, _("RECORDING"));

                ImGui::Spacing();

                if (!assume_recording)
                {
                    if (ImGui::Button(_("Start###startrecording")))
                        start_recording();
                }
                else
                {
                    if (ImGui::Button(_("Stop##stoprecording")))
                        stop_recording();
                }
            }

            show_tracking = ImGui::CollapsingHeader(_("Tracking"), tracking_started_cli ? ImGuiTreeNodeFlags_DefaultOpen : ImGuiTreeNodeFlags_None);
            if (show_tracking)
            {
                try_init_tracking_widget();
                tracking_widget->render();
            }

            if (ImGui::CollapsingHeader(_("VFOs")))
            {
                vfos_mtx.lock();
                std::string to_delete = "";
                if (vfo_list.size() == 0)
                {
                    const char *no_vfo_text = _("No Active VFOs");
                    ImGui::SetCursorPosX(ImGui::GetContentRegionAvail().x / 2 - ImGui::CalcTextSize(no_vfo_text).x / 2);
                    ImGui::TextDisabled("%s", no_vfo_text);
                }
                else
                {
                    for (auto &vfo : vfo_list)
                    {
                        ImGui::PushStyleColor(ImGuiCol_Text, style::theme.green.Value);
                        ImGui::SeparatorText(vfo.name.c_str());
                        ImGui::PopStyleColor();
                        ImGui::BulletText(_("Frequency: %s"), format_notated(vfo.freq, "Hz").c_str());
                        if (vfo.selected_pipeline.id != "")
                        {
                            ImGui::BulletText(_("Pipeline: %s"), vfo.selected_pipeline.id.c_str());
                            ImGui::BulletText(_("Directory: %s"), vfo.output_dir.c_str());
                        }
                        else if (vfo.file_sink)
                        {
                            if (vfo.file_sink->get_written() < 1e9)
                                ImGui::BulletText(_("IQ Size: %.2f MB"), vfo.file_sink->get_written() / 1e6);
                            else
                                ImGui::BulletText(_("IQ Size: %.2f GB"), vfo.file_sink->get_written() / 1e9);
                        }
                        ImGui::Spacing();
                        if (ImGui::Button(std::string(_("Stop##") + vfo.id).c_str(), ImVec2(ImGui::GetContentRegionAvail().x, 0)))
                            to_delete = vfo.id;
                    }
                }
                vfos_mtx.unlock();
                if (to_delete != "")
                    del_vfo(to_delete);
    #if 0
                if (ImGui::Button("Add Test"))
                {
                    int idp = 0;
                    for (int n = 0; n < (int)pipeline::pipelines.size(); n++)
                    {
                        if ("meteor_m2-x_lrpt" == pipeline::pipelines[n].id)
                            idp = n;
                    }

                    add_vfo_live("meteor_test", "Meteor Test", 137.9e6, pipeline::pipelines[idp], {});
                }

                if (ImGui::Button("Del Test"))
                {
                    del_vfo("meteor_test");
                }
    #endif
            }

            if (ImGui::CollapsingHeader(_("Debug")))
            {
                if (constellation_debug == nullptr)
                    constellation_debug = new widgets::ConstellationViewer();
                if (is_started)
                    constellation_debug->pushComplex(source_ptr->output_stream->readBuf, 256);
                constellation_debug->draw();
            }
        }

        eventBus->fire_event<RecorderDrawPanelEvent>({});
    }

    void RecorderApplication::drawContents(ImVec2 win_size)
    {
        fswitch.draw_window();

        float right_width = win_size.x;
        float wf_size = win_size.y;
        float left_width = ImGui::GetCursorPosX();

        float recorder_size_x = right_width;

        float _live_height = 250;

        float wf_size_offset = 0;
        if (is_processing && !processing_modules_floating_windows)
            wf_size_offset = _live_height * ui_scale;
        if (vfo_list.size() > 0)
        {
            wf_size_offset = 270 * ui_scale;
            _live_height = 270;
        }

        wf_size = win_size.y - wf_size_offset; // + 13 * ui_scale;

        // float wf_size_offset = 0;
        // if (is_processing && !processing_modules_floating_windows)
        //     wf_size_offset = 250 * ui_scale;
        // if (vfo_list.size() > 0)
        //     wf_size_offset = 270 * ui_scale;

        // ImVec2 recorder_size = ImGui::GetContentRegionAvail();
        // float wf_size_offset = 0;
        // if (is_processing && !processing_modules_floating_windows)
        //     wf_size_offset = 250 * ui_scale;
        // if (vfo_list.size() > 0)
        //     wf_size_offset = 270 * ui_scale;
        // float wf_size = recorder_size.y - wf_size_offset; // + 13 * ui_scale;

        {
            ImGui::BeginChild("RecorderFFT", {right_width, wf_size}, false, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
            if (fswitch.is_running())
            {
                // While switching, the folded result replaces spectrum and waterfall
                ImGui::SetCursorScreenPos(ImVec2(left_width + 8, 25 * ui_scale));
                fswitch.draw_plot(ImVec2(right_width - 9 * ui_scale, wf_size - 8 * ui_scale));
                if (fswitch.take_stop_request())
                    stop_fswitch(); // back to the normal spectrum and acquisition
                std::string why;
                if (fswitch.take_fallback(why))
                {
                    stop_fswitch();
                    fswitch_error = "Switching on the server not possible: " + why;
                    logger->error("Frequency switching: " + fswitch_error);
                }
            }
            else
            {
                // Fix for "could not allocate waterfall memory"
                if (waterfall_ratio > 0.9)
                    waterfall_ratio = 0.9;

                float fft_height = wf_size * (show_waterfall ? waterfall_ratio : 1.0);
                float wf_height = wf_size * (1 - waterfall_ratio) + 15 * ui_scale;
                float wfft_widht = right_width - 9 * ui_scale;
#ifdef __ANDROID__
                int offset = 8;
#else
                int offset = 30;
#endif
                // ImGui::SetNextWindowSizeConstraints(ImVec2((right_width + offset * ui_scale), 50), ImVec2((right_width + offset * ui_scale), wf_size));
                //                   ImGui::SetNextWindowSize(ImVec2((right_width + offset * ui_scale), show_waterfall ? waterfall_ratio * wf_size : wf_size));
                ImGui::SetCursorScreenPos(ImVec2(left_width + 8, 25 * ui_scale));
                // The #fft child is resizable (ResizeY), so ImGui only uses its size
                // argument on first use. Force the height while the waterfall is hidden
                // (FFT and its frequency scale fill the area), and once when it is shown
                // again (back to the height set by waterfall_ratio).
                if (!show_waterfall || !last_show_waterfall)
                    ImGui::SetNextWindowSize(ImVec2((right_width + offset * ui_scale), show_waterfall ? waterfall_ratio * wf_size : wf_size));
                last_show_waterfall = show_waterfall;
                if (ImGui::BeginChild("#fft", ImVec2((right_width + offset * ui_scale), (show_waterfall ? waterfall_ratio * wf_size : wf_size)), ImGuiChildFlags_ResizeY, //
                                      ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_ChildWindow | ImGuiWindowFlags_NoScrollbar |
                                          ImGuiWindowFlags_NoScrollWithMouse))
                {
                    fft_plot->draw({float(wfft_widht), fft_height});
                    if (show_waterfall) // && ImGui::IsMouseDragging(ImGuiMouseButton_Left))
                        waterfall_ratio = ImGui::GetWindowHeight() / wf_size;
                    if (ImGui::IsWindowHovered())
                    {
                        ImVec2 mouse_pos = ImGui::GetMousePos();
                        float ratio = (mouse_pos.x - left_width - 16 * ui_scale) / (right_width - 9 * ui_scale) - 0.5;
                        ImGui::SetTooltip(
                            "%s",
                            ((ratio >= 0 ? "" : "- ") + format_notated(abs(ratio) * get_samplerate(), "Hz\n") + format_notated(source_ptr->get_frequency() + ratio * get_samplerate(), "Hz")).c_str());
                    }

                    ImGui::EndChild();
                }

                if (show_waterfall)
                {
                    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 8 * ui_scale);
                    ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 4 * ui_scale);
                    waterfall_plot->draw({wfft_widht, wf_height}, is_started);
                }
            }
            ImGui::EndChild();
        }

        if (vfo_list.size() > 0)
        {
            if (ImGui::BeginTabBar("VFO TabBar", ImGuiTabBarFlags_None))
            {
                if (is_processing)
                {
                    if (ImGui::BeginTabItem(_("Live Processing")))
                    {
                        float y_pos = ImGui::GetCursorPosY(); //+ 35 * ui_scale;
                        float live_width = recorder_size_x + 16 * ui_scale;
                        float live_height = _live_height * ui_scale;
                        float winwidth = live_pipeline->modules.size() > 0 ? live_width / live_pipeline->modules.size() : live_width;
                        float currentPos = ImGui::GetCursorPosX();
                        ImGui::PushStyleColor(ImGuiCol_TitleBgActive, ImGui::GetStyleColorVec4(ImGuiCol_TitleBg));
                        for (std::shared_ptr<pipeline::ProcessingModule> &module : live_pipeline->modules)
                        {
                            ImGui::SetNextWindowPos({currentPos, y_pos});
                            ImGui::SetNextWindowSize({(float)winwidth, (float)live_height});
                            module->drawUI(false);
                            currentPos += winwidth;

                            // if (ImGui::GetCurrentContext()->last_window != NULL)
                            //     currentPos += ImGui::GetCurrentContext()->last_window->Size.x;
                            //  logger->info(ImGui::GetCurrentContext()->last_window->Name);
                        }
                        ImGui::PopStyleColor();
                        ImGui::EndTabItem();
                    }
                }

                vfos_mtx.lock();
                for (auto &vfo : vfo_list)
                {
                    if (ImGui::BeginTabItem(vfo.name.c_str()))
                    {
                        if (vfo.selected_pipeline.id != "")
                        {
                            float y_pos = ImGui::GetCursorPosY(); //+ 35 * ui_scale;
                            float live_width = recorder_size_x + 16 * ui_scale;
                            float live_height = _live_height * ui_scale;
                            float winwidth = vfo.live_pipeline->modules.size() > 0 ? live_width / vfo.live_pipeline->modules.size() : live_width;
                            float currentPos = ImGui::GetCursorPosX();
                            ;
                            ImGui::PushStyleColor(ImGuiCol_TitleBgActive, ImGui::GetStyleColorVec4(ImGuiCol_TitleBg));
                            for (std::shared_ptr<pipeline::ProcessingModule> &module : vfo.live_pipeline->modules)
                            {
                                ImGui::SetNextWindowPos({currentPos, y_pos});
                                ImGui::SetNextWindowSize({(float)winwidth, (float)live_height});
                                module->drawUI(false);
                                currentPos += winwidth;

                                // if (ImGui::GetCurrentContext()->last_window != NULL)
                                //     currentPos += ImGui::GetCurrentContext()->last_window->Size.x;
                                //  logger->info(ImGui::GetCurrentContext()->last_window->Name);
                            }
                            ImGui::PopStyleColor();
                        }
                        ImGui::EndTabItem();
                    }
                }
                vfos_mtx.unlock();

                ImGui::EndTabBar();
            }
        }
        else if (is_processing)
        {
            if (processing_modules_floating_windows)
            {
                for (std::shared_ptr<pipeline::ProcessingModule> &module : live_pipeline->modules)
                    module->drawUI(true);
            }
            else
            {
                float y_pos = ImGui::GetCursorPosY(); //+ 35 * ui_scale;
                float live_width = recorder_size_x + 16 * ui_scale;
                float live_height = _live_height * ui_scale;
                float winwidth = live_pipeline->modules.size() > 0 ? live_width / live_pipeline->modules.size() : live_width;
                float currentPos = ImGui::GetCursorPosX();
                ImGui::PushStyleColor(ImGuiCol_TitleBgActive, ImGui::GetStyleColorVec4(ImGuiCol_TitleBg));
                for (std::shared_ptr<pipeline::ProcessingModule> &module : live_pipeline->modules)
                {
                    ImGui::SetNextWindowPos({currentPos, y_pos});
                    ImGui::SetNextWindowSize({(float)winwidth, (float)live_height});
                    module->drawUI(false);
                    currentPos += winwidth;

                    // if (ImGui::GetCurrentContext()->last_window != NULL)
                    //     currentPos += ImGui::GetCurrentContext()->last_window->Size.x;
                    //  logger->info(ImGui::GetCurrentContext()->last_window->Name);
                }
                ImGui::PopStyleColor();
            }
        }

        // Keyboard shortcuts
        {
            // FFT
            if (ImGui::IsKeyDown(ImGuiKey_ReservedForModShift) && ImGui::IsKeyDown(ImGuiKey_X))
            {
                if (ImGui::IsKeyDown(ImGuiKey_DownArrow))
                    fft_plot->scale_max -= 0.5;
                else if (ImGui::IsKeyDown(ImGuiKey_UpArrow))
                    fft_plot->scale_max += 0.5;
                else if (ImGui::IsKeyDown(ImGuiKey_LeftArrow))
                    fft_plot->scale_min -= 0.5;
                else if (ImGui::IsKeyDown(ImGuiKey_RightArrow))
                    fft_plot->scale_min += 0.5;
#if 0
                else if (ImGui::IsKeyDown(ImGuiKey_PageUp))
                    fft->avg_rate += 0.001;
                else if (ImGui::IsKeyDown(ImGuiKey_PageDown))
                    fft->avg_rate -= 0.001;

                if (fft->avg_rate > 0.99)
                    fft->avg_rate = 0.99;
                else if (fft->avg_rate < 0.01)
                    fft->avg_rate = 0.01;
#endif
            }
        }
    }
}; // namespace satdump
