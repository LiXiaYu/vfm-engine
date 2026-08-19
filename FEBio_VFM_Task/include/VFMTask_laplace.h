#pragma once

#include "VFMTask.h"

void L_externalVirtualWork_LaplcaeTransform(
    std::vector<std::vector<double>>& externalVirtualWork,
    std::vector<std::vector<Eigen::dcomplex>>& externalVirtualWork_laplace,
    std::vector<double>& timeArray,
    VFMTask* task);

void L_StressPK2_withJc_LaplaceTransform(
    std::vector<double>& timeArray,
    VFMTask* task,
    FEModel* fem,
    std::vector<int>& solution_elementsDomainID,
    std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray,
    std::vector<std::vector<std::vector<double>>>& trueJArray,
    std::vector<std::vector<std::vector<std::vector<Eigen::dcomplex>>>>& Sepk2_dotdot_vStrain,
    std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

std::vector<std::vector<std::complex<double>>> L_fs(
    std::vector<std::vector<double>>& externalVirtualWork,
    std::vector<double>& timeArray,
    VFMTask* task,
    FEModel* fem,
    std::vector<int>& solution_elementsDomainID,
    std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray,
    std::vector<std::vector<std::vector<double>>>& trueJArray,
    std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

std::vector<std::vector<double>> inv_ft(
    std::vector<std::vector<std::complex<double>>> fs,
    std::vector<double> timeArray);

std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<int>,
    std::vector<std::vector<std::vector<mat3ds>>>,
    std::vector<std::vector<std::vector<std::vector<double>>>>>
cal_internal_normal_linerE_vw(
    std::vector<double>& timeArray,
    VFMTask* task,
    FEModel* fem,
    std::vector<int>& solution_elementsDomainID,
    std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray,
    std::vector<std::vector<std::vector<double>>>& trueJArray,
    std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV);

std::pair<int, int> find_last_two_cycles(const std::span<double>& ts, double bpm);
