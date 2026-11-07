#pragma once

#include "dll_export.h"
#include <string>

namespace satdump
{
    SATDUMP_DLL extern const std::string RESOURCES_PATH;
    SATDUMP_DLL extern const std::string LIBRARIES_PATH;
    SATDUMP_DLL extern const std::string SATDUMP_VERSION;
    SATDUMP_DLL extern const std::string SATDUMP_VERSION_TAG;

    SATDUMP_DLL extern std::string RESPATH;
    SATDUMP_DLL extern std::string LIBPATH;

    std::string getSatDumpVersionName();

    // Name and logo (resource file) shown on the loading and welcome
    // screens. A derived application (AstroDump) sets them at start
    SATDUMP_DLL extern std::string APP_NAME;
    SATDUMP_DLL extern std::string APP_ICON;
    // Path of APP_ICON, or of the SatDump logo if that file is missing
    SATDUMP_DLL std::string getAppIconPath();
} // namespace satdump