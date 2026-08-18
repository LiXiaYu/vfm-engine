#pragma once

#include "VFMTask.h"

bool write_tecplot_nodal_position_and_displacement(
    const std::string& filepath,
    const std::span<double>& time,
    const std::vector<std::vector<double>>& initialCoordinate,
    const span3d<double>& timeDisplacement,
    const std::vector<int>& selected_node_list,
    FEMesh& mesh,
    const std::vector<int>& selected_element_list,
    int start_index,
    int end_index_inclusive);
