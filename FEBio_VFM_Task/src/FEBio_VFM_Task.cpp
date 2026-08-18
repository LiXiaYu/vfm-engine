// FEBio_VFM_Task.cpp: 定义应用程序的入口点。
//
#include "FEBio_VFM_Task.h"

FECORE_EXPORT unsigned int GetSDKVersion()
{
    //return 198656;
    return FE_SDK_VERSION;
}

FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
    FECoreKernel::SetInstance(&febio);

    REGISTER_FECORE_CLASS(VFMTask, "VFM");
}

FECORE_EXPORT void PluginCleanup()
{
    // Clean up plugin resources 
}

FECORE_EXPORT void GetPluginVersion(int& major, int& minor, int& patch)
{
    major = 1;
    minor = 1; // mmap
    patch = 2;
}
