#pragma once

#include "FEBio_refunction.h"
#include "common_FEBio.h"
#include "VFMTask.h"
#include "optim.h"

bool read_inited_information(FEModel* fem, unsigned int when, void* pd);

bool read_stepsolved_information(FEModel* fem, unsigned int when, void* pd);

bool everytimestep_withinited_savedata(FEModel* fem, unsigned int when, void* pd);

bool read_solved_information(FEModel* fem, unsigned int when, void* pd);

void L_externalVirtualWork_LaplcaeTransform(std::vector<std::vector<double>>& externalVirtualWork, std::vector<std::vector<Eigen::dcomplex>>& externalVirtualWork_laplace, std::vector<double>& timeArray, VFMTask* task);

void L_StressPK2_withJc_LaplaceTransform(std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<Eigen::dcomplex>>>>& Sepk2_dotdot_vStrain, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

::std::vector<::std::vector<::std::complex<double>>> L_fs(std::vector<std::vector<double>>& externalVirtualWork, std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

::std::vector<::std::vector<double>> inv_ft(::std::vector<::std::vector<::std::complex<double>>> fs, ::std::vector<double> timeArray);

::std::tuple<::std::vector<::std::vector<double>>, ::std::vector<::std::vector<double>>, ::std::vector<int>, ::std::vector<::std::vector<::std::vector<mat3ds>>>, ::std::vector<::std::vector<::std::vector<::std::vector<double>>>>> cal_internal_normal_linerE_vw(std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

::std::pair<int, int> find_last_two_cycles(const ::std::span<double>& ts, double bpm);

static bool write_tecplot_nodal_position_and_displacement(const std::string& filepath, const std::span<double>& time, const std::vector<std::vector<double>>& initialCoordinate, const span3d<double>& timeDisplacement, const std::vector<int>& selected_node_list, FEMesh& mesh, const std::vector<int>& selected_element_list, int start_index, int end_index_inclusive);