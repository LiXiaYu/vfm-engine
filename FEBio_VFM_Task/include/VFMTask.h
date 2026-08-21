#pragma once

#include "common_FEBio.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

class VFMTask_configure
{
public:
    double bpm = 0.0;
    int start_index = 0;
    int end_index = 0;
    bool run_febio_solve = true;
    bool reuse_saved_result_buffer = false;
    bool isReadfromsaveOptimfunc = false;
    bool isSetDisplacmentAndPressure = false;
    bool allow_failed_solve_postprocessing = false;
    ::std::string FEBio_dump_path;
    ::std::string Optim_dump_path;

    bool isReadfromfeblogfile = false;
    ::std::string feblog_uxuyuz_path;
    ::std::string feblog_sxsyszsxysyzsxz_path;

    bool isLaplaceVFM = false;

    int InterpolatedDisplacementTimes = 0;

    int displacementdataNumber = 0;
    int stressdataNumber = 0;
    int pressurevalueNumber = 0;


    std::vector<std::tuple<std::string, std::string>> solution;
    std::vector<std::tuple<std::string, std::string>> pressure_load;
    std::vector<std::tuple<std::string, std::string>> constraint_load;
    std::vector<std::tuple<std::string, std::string>> fixed;
    ::std::vector<int> fixednode;

    ::std::vector<::std::vector<double>> initialCoordinate;

    ::std::vector<::std::function<::std::vector<double>(const ::std::vector<double>&, int)>> vf_u_functions;
    
    // python: get_timestep, set_timestep
    span<double> timestep;
    ::std::optional<::std::vector<double>> pending_timestep;

    ::std::function<::std::vector<int>(const int, const ::std::vector<int>&, const ::std::vector<::std::vector<double>>&)> select_solution_element_function;
    ::std::function<::std::vector<double>(const ::std::vector<double>&, int, bool)> setdisplacement_function;
    ::std::function<::std::vector<double>(const ::std::vector<::std::vector<double>>&, int)> setstress_function;

    int optim_function_iter_max = 50;

    ::std::function<double(const ::std::vector<double>&)> optim_function;
    ::std::function<double(const ::std::vector<double>&)> optim_function_T;
    ::std::function<::std::vector<::std::vector<double>>(double gamma, double tau)> optim_function_E__gamma_tau;
    ::std::function<::std::vector<::std::vector<double>>(double gamma, double tau)> optim_function_gamma__E_tau;
    ::std::function<::std::vector<::std::vector<double>>(double gamma, double tau)> optim_function_tau__E_gamma;
    ::std::function<double(const ::std::vector<double>&)> optim_function_E_gamma_tau;
    ::std::function<double(const ::std::vector<double>&)> optim_function_E_gamma_tau_invF;
    ::std::function<double(const ::std::vector<double>&)> optim_function_Et;
    ::std::function<double(const ::std::vector<double>&)> optim_function_elastic_E;

    // One-time diagnostic data prepared while the optimization callbacks are
    // initialized. Per-evaluation data is controlled separately below.
    bool optim_function_output_debug_info = false;
    bool optim_function_output_iteration_stress = false;
    bool optim_function_output_iteration_virtual_work = false;

    // OpenMP execution used inside the C++ objective functions. A thread count
    // of zero selects the OpenMP runtime default; one forces serial execution.
    bool optim_function_parallel = true;
    int optim_function_num_threads = 0;

    bool isSetNewLaplaceVFMs = false;
    ::std::vector<::std::complex<double>> LaplaceVFM_s;
    bool isTalbotLaplaceVFM_s = false;
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

    span3d<double> timedisplacement;
    span3d<double> timestress;
    span3d<double> timenodalforce;
    span2d<double> timeconstraintpressure;
    span2d<uint8_t> timeconstraintactivate; // bool
    
    ::std::vector<::std::vector<int>> forceFacet;

    int nVirtualFields;

    ::std::vector<int> solution_elementsID;

// python interperter
public:
    pybind11::scoped_interpreter guard{};
    pybind11::module_ pyfile_module;
// input parameters
public:
    VFMTask_configure configure;

    ::std::vector<double> objects_values;

    ::std::vector<int> readfromsaveOptimfunc_index;

public:
    // mmap
    mio::mmap_sink mmap;

    size_t total_steps = 0;
    size_t recorded_steps = 0;
    size_t ndisplacement = 0;
    size_t nstress = 0;
    size_t nnodalforce = 0;
    size_t nconstraint = 0; // pressure and activate

    void reinit_mmap_and_spans();
};


std::vector<std::vector<double>> createGrid(const std::vector<double>& xl, const std::vector<double>& xu, const std::vector<double>& xd, int now_index = -1, ::std::vector<double>* now_xs = nullptr);

