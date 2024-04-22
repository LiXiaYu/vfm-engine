#pragma once

#include "common_FEBio.h"

#include <pybind11/embed.h>

class VFMTask_configure
{
public:
    double bpm = 0.0;
    bool isRead_FEMresult_fromsavefile = false;
    bool isReadfromsaveOptimfunc = false;
    int displacementdataNumber = 0;
    int stressdataNumber = 0;
    int pressurevalueNumber = 0;
    std::vector<std::tuple<std::string, std::string>> solution;
    std::vector<std::tuple<std::string, std::string>> pressure_load;
    std::vector<std::tuple<std::string, std::string>> fixed;
};

class VFMTask : public FECoreTask
{
public:
    VFMTask(FEModel* pfem);

    ~VFMTask();

    bool Init(const char* szfile);

    bool Run();

    ::std::string szfile;

    ::std::string outlogfile;
    ::std::string outSavefile;
    ::std::string dumpfile;
public:
    ::std::vector<vec3d> nodes;
    ::std::vector<::std::vector<double>> initialCoordinate;


    ::std::vector<double> timestep;
    ::std::vector<::std::vector<::std::vector<double>>> timedisplacement;
    ::std::vector<::std::vector<::std::vector<double>>> timestress;

    ::std::vector<int> fixednode;
    
    ::std::vector<::std::vector<int>> forceFacet;

    int nVirtualFields;

    ::std::vector<int> solution_elementsID;

public:
    pybind11::scoped_interpreter guard{};
// input parameters
public:
    VFMTask_configure configure;

    ::std::vector<double> objects_values;

    ::std::vector<int> readfromsaveOptimfunc_index;
};


std::vector<std::vector<double>> createGrid(const std::vector<double>& xl, const std::vector<double>& xu, const std::vector<double>& xd, int now_index = -1, ::std::vector<double>* now_xs = nullptr);


bool everytimestep_withinited_savedata(FEModel* fem, unsigned int when, void* pd);