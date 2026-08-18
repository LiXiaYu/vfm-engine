"""VFM configuration for the trabecular-meshwork DIC validation case."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import List

import numpy as np
from scipy.optimize import differential_evolution

import VFMTask_pybind


CASE_DIR = Path(__file__).resolve().parents[1]
POSITIONS_PATH = Path(
    os.environ.get(
        "VFM_DIC_POSITIONS",
        CASE_DIR / "input" / "data" / "dic_positions.npy",
    )
).resolve()
METRICS_PATH = Path(
    os.environ.get("VFM_VALIDATION_METRICS", Path.cwd() / "validation_metrics.json")
).resolve()
VALIDATION_SEED = int(os.environ.get("VFM_VALIDATION_SEED", "20260818"))


def InitVFMTask(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    configure = VFMTask_pybind.VFMTask_configure()

    configure.bpm = 15  # Retained from the legacy demo; currently ineffective.
    configure.isRead_FEMresult_fromsavefile = False
    configure.isReadfromsaveOptimfunc = False
    configure.isSetDisplacmentAndPressure = True
    configure.isLaplaceVFM = False
    configure.isSetNewLaplaceVFMs = False
    configure.isTalbotLaplaceVFM_s = False

    configure.displacementdataNumber = 0
    configure.stressdataNumber = 1
    configure.solution = [("elements", "solid")]
    configure.pressure_load = [("FEPressureLoad", "Pressure1")]
    configure.fixed = [("surface", "NormalDisplacement1")]
    configure.optim_function_output_debug_info = True
    return configure


def SetPeriodTimeIndex(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    configure.start_index = 0
    configure.end_index = 1
    return configure


def SetSelectSolutionElementFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    def select_solution_element_function(element_index, nodeids, node_xyzs):
        return []

    configure.select_solution_element_function = select_solution_element_function
    return configure


def SetDisplacmentAndPressureFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    if not POSITIONS_PATH.is_file():
        raise FileNotFoundError(f"DIC positions were not found: {POSITIONS_PATH}")

    position_array = np.load(POSITIONS_PATH, allow_pickle=False)
    if not (
        isinstance(position_array, np.ndarray)
        and position_array.ndim == 3
        and position_array.shape[2] == 2
    ):
        raise RuntimeError("DIC positions must have shape (T+1, N, 2)")

    time_array = [0.0, 0.02 * 9]  # Two selected frames from a 50 fps sequence.
    pressure_array = [0.0, 0.004]  # MPa.

    def setstress(coordinates, index_timestep):
        pressure = pressure_array[index_timestep] / 3.0
        return [pressure, pressure, pressure, 0.0, 0.0, 0.0]

    def setdisplacement(
        coordinate: List[float], index_timestep: int, debugtag: bool
    ) -> List[float]:
        point = np.asarray(coordinate[:2], dtype=float)
        frame_count = position_array.shape[0]
        point_count = position_array.shape[1]
        timestep = min(max(int(index_timestep), 0), frame_count - 1)

        initial_positions = position_array[0]
        current_positions = position_array[timestep]
        distances = np.linalg.norm(initial_positions - point.reshape((1, 2)), axis=1)

        exact_matches = np.where(distances < 1e-8)[0]
        if exact_matches.size:
            index = int(exact_matches[0])
            displacement = current_positions[index] - initial_positions[index]
            result = [float(displacement[0]), float(displacement[1]), 0.0]
            if debugtag:
                print("setdisplacement exact match:", index, result)
            return result

        neighbour_count = min(8, point_count)
        neighbour_indices = np.argsort(distances)[:neighbour_count]
        neighbour_distances = distances[neighbour_indices].copy()
        neighbour_distances[neighbour_distances < 1e-12] = 1e-12
        weights = 1.0 / neighbour_distances**2
        weights /= np.sum(weights)

        displacement_vectors = (
            current_positions[neighbour_indices] - initial_positions[neighbour_indices]
        )
        ux = float(np.sum(weights * displacement_vectors[:, 0]))
        uy = float(np.sum(weights * displacement_vectors[:, 1]))
        result = [ux, uy, 0.0]
        if debugtag:
            print("setdisplacement IDW:", neighbour_count, result)
        return result

    configure.set_timestep(time_array)
    configure.setdisplacement_function = setdisplacement
    configure.setstress_function = setstress
    configure.InterpolatedDisplacementTimes = 2
    return configure


X_SCALE = 1e-8
Y_SCALE = 1e-8
Z_SCALE = 1e-8


def _deterministic_component(index_coordinate: int, component: int) -> float:
    """Return a value independent of callback order and thread scheduling."""
    local_seed = VALIDATION_SEED + 1_000_003 * component + int(index_coordinate)
    return float(np.random.default_rng(local_seed).uniform(0.0, 1.0))


def SetVirtualDisplacementFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    fixed_nodes = set(configure.fixednode)

    def vf_u_random_x(coordinate, index_coordinate):
        if index_coordinate in fixed_nodes:
            return [0.0, 0.0, 0.0]
        return [
            _deterministic_component(index_coordinate, 0) * X_SCALE,
            0.0,
            0.0,
        ]

    def vf_u_random_y(coordinate, index_coordinate):
        if index_coordinate in fixed_nodes:
            return [0.0, 0.0, 0.0]
        return [
            0.0,
            _deterministic_component(index_coordinate, 1) * Y_SCALE,
            0.0,
        ]

    def vf_u_random_z(coordinate, index_coordinate):
        if index_coordinate in fixed_nodes:
            return [0.0, 0.0, 0.0]
        return [
            0.0,
            0.0,
            _deterministic_component(index_coordinate, 2) * Z_SCALE,
        ]

    configure.vf_u_functions = [vf_u_random_x, vf_u_random_y, vf_u_random_z]
    return configure


def BeforeOptim(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    def objective_function(params):
        return configure.optim_function_elastic_E([float(params[0])])

    result = differential_evolution(
        objective_function,
        [(0.00001, 9.0)],
        disp=True,
        seed=VALIDATION_SEED,
        workers=1,
        updating="immediate",
    )

    metrics = {
        "schema_version": 1,
        "case": "trabecular_meshwork_dic",
        "seed": VALIDATION_SEED,
        "elastic_E": {
            "x": float(result.x[0]),
            "fun": float(result.fun),
            "success": bool(result.success),
            "message": str(result.message),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
        },
    }
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"Validation metrics written to {METRICS_PATH}")
    return configure
