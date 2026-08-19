#include "VFMTask_debug_output.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
bool write_tecplot_nodal_position_and_displacement(
	const std::string& filepath,
	const std::span<double>& time,
	const std::vector<std::vector<double>>& initialCoordinate,
	const span3d<double>& timeDisplacement,
	const std::vector<int>& selected_node_list,
    FEMesh& mesh,
	const std::vector<int>& selected_element_list,
	int start_index,
	int end_index_inclusive)
{
	std::ofstream out(filepath, std::ios::out | std::ios::trunc);
	if (!out.is_open()) return false;

    const int nNodes = static_cast<int>(initialCoordinate.size());

	out << "TITLE = \"NodalPositionAndDisplacement\"\n";
	out << "VARIABLES = \"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"Time\",\"NodeId\"\n";

   out << std::setprecision(16);

	// determine which nodes to output (global node indices)
	std::vector<int> output_node_ids;
	output_node_ids.reserve(!selected_node_list.empty() ? selected_node_list.size() : static_cast<size_t>(nNodes));
	if (!selected_node_list.empty())
	{
		output_node_ids.assign(selected_node_list.begin(), selected_node_list.end());
	}
	else
	{
		for (int nid = 0; nid < nNodes; ++nid) output_node_ids.push_back(nid);
	}

	const int output_nodes = static_cast<int>(output_node_ids.size());

	// build global->local node index map (local index is 0-based)
	std::vector<int> global_to_local(nNodes, -1);
	for (int i = 0; i < output_nodes; ++i)
	{
		const int gid = output_node_ids[i];
		if (gid >= 0 && gid < nNodes) global_to_local[gid] = i;
	}

	// build element connectivity using local node indices (Tecplot uses 1-based indices)
	std::vector<std::vector<int>> conn_hex8;
	std::vector<std::vector<int>> conn_penta6;
	std::vector<std::vector<int>> conn_tet4;
	std::vector<std::vector<int>> conn_tri3;
	std::vector<std::vector<int>> conn_seg2;

	const bool use_element_selection = !selected_element_list.empty();
	if (use_element_selection)
	{
		for (int ei = 0; ei < (int)selected_element_list.size(); ++ei)
		{
			const int element_index = selected_element_list[ei];
			FEElement* pe = mesh.Element(element_index);
			if (pe == nullptr) continue;

			const int nn = pe->Nodes();
			std::vector<int> conn(nn);
			bool ok = true;
			for (int k = 0; k < nn; ++k)
			{
				const int gid = pe->m_node[k];
				if (gid < 0 || gid >= nNodes) { ok = false; break; }
				const int lid = global_to_local[gid];
				if (lid < 0) { ok = false; break; }
				conn[k] = lid + 1;
			}
			if (!ok) continue;

			switch (nn)
			{
			case 8: conn_hex8.push_back(std::move(conn)); break;
			case 6: conn_penta6.push_back(std::move(conn)); break;
			case 4: conn_tet4.push_back(std::move(conn)); break;
			case 3: conn_tri3.push_back(std::move(conn)); break;
			case 2: conn_seg2.push_back(std::move(conn)); break;
			default:
				// unsupported element type for Tecplot FE zone in this exporter
				break;
			}
		}
	}

    auto write_zone = [&](int it_step, const double t, const char* zone_suffix, const char* zone_type, const std::vector<std::vector<int>>& conn, int nodes_per_elem)
	{
		if (conn.empty()) return;

		out << "ZONE T=\"t=" << t << " " << zone_suffix << "\""
			<< ", N=" << output_nodes
			<< ", E=" << static_cast<int>(conn.size())
			<< ", ZONETYPE=" << zone_type
			<< ", DATAPACKING=POINT\n";

		for (int local_i = 0; local_i < output_nodes; ++local_i)
		{
			const int gid = output_node_ids[local_i];
            const double u = timeDisplacement[it_step][gid][0];
			const double v = timeDisplacement[it_step][gid][1];
			const double w = timeDisplacement[it_step][gid][2];

			const double x = initialCoordinate[gid][0] + u;
			const double y = initialCoordinate[gid][1] + v;
			const double z = initialCoordinate[gid][2] + w;

			out << x << " " << y << " " << z << " "
				<< u << " " << v << " " << w << " "
				<< t << " " << (gid + 1) << "\n";
		}

		for (const auto& econn : conn)
		{
			if ((int)econn.size() != nodes_per_elem) continue;
			for (int k = 0; k < nodes_per_elem; ++k)
			{
				out << econn[k];
				out << ((k == nodes_per_elem - 1) ? '\n' : ' ');
			}
		}
	};

	for (int it = start_index; it <= end_index_inclusive; ++it)
	{
		const double t = time[it];

		// If no element selection provided (or none matched), still output a node-only zone
		if (!use_element_selection || (conn_hex8.empty() && conn_penta6.empty() && conn_tet4.empty() && conn_tri3.empty() && conn_seg2.empty()))
		{
			out << "ZONE T=\"t=" << t << " scatter\""
				<< ", I=" << output_nodes
				<< ", J=1, K=1"
				<< ", DATAPACKING=POINT\n";

			for (int local_i = 0; local_i < output_nodes; ++local_i)
			{
				const int gid = output_node_ids[local_i];
				const double u = timeDisplacement[it][gid][0];
				const double v = timeDisplacement[it][gid][1];
				const double w = timeDisplacement[it][gid][2];

				const double x = initialCoordinate[gid][0] + u;
				const double y = initialCoordinate[gid][1] + v;
				const double z = initialCoordinate[gid][2] + w;

				out << x << " " << y << " " << z << " "
					<< u << " " << v << " " << w << " "
					<< t << " " << (gid + 1) << "\n";
			}
			continue;
		}

     write_zone(it, t, "HEX8", "FEBRICK", conn_hex8, 8);
		write_zone(it, t, "PENTA6", "FEPRISM", conn_penta6, 6);
		write_zone(it, t, "TET4", "FETETRAHEDRON", conn_tet4, 4);
		write_zone(it, t, "TRI3", "FETRIANGLE", conn_tri3, 3);
		write_zone(it, t, "SEG2", "FELINESEG", conn_seg2, 2);
	}

	return true;
}
