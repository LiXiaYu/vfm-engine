"""Analytical Gaussian-integration checks executed inside the FEBio VFM task."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

import VFMTask_pybind


ELEMENT_TYPE = os.environ["VFM_GAUSS_ELEMENT_TYPE"]
METRICS_PATH = Path(os.environ["VFM_GAUSS_METRICS"]).resolve()
VIRTUAL_SCALE = 1.0e-4
STRESS = 2.0


def InitVFMTask(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    configure = VFMTask_pybind.VFMTask_configure()
    configure.run_febio_solve = True
    configure.reuse_saved_result_buffer = False
    configure.isReadfromsaveOptimfunc = False
    configure.isSetDisplacmentAndPressure = True
    configure.allow_failed_solve_postprocessing = False
    configure.isLaplaceVFM = False
    configure.isSetNewLaplaceVFMs = False
    configure.isTalbotLaplaceVFM_s = False
    configure.displacementdataNumber = 0
    configure.stressdataNumber = 1
    configure.pressurevalueNumber = 0
    configure.solution = [("elements", "Solid")]
    configure.pressure_load = [("FEPressureLoad", "Pressure1")]
    configure.fixed = []
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
    configure.select_solution_element_function = (
        lambda element_index, nodeids, node_xyzs: []
    )
    return configure


def SetDisplacmentAndPressureFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    configure.set_timestep([0.0, 1.0])

    def setdisplacement(coordinate, index_timestep, debugtag):
        return [0.0, 0.0, 0.0]

    def setstress(coordinates, index_timestep):
        value = 0.0 if int(index_timestep) == 0 else STRESS
        return [value, value, value, 0.0, 0.0, 0.0]

    configure.setdisplacement_function = setdisplacement
    configure.setstress_function = setstress
    configure.InterpolatedDisplacementTimes = 2
    return configure


def SetVirtualDisplacementFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    def linear_x(coordinate, index_coordinate):
        return [VIRTUAL_SCALE * float(coordinate[0]), 0.0, 0.0]

    configure.vf_u_functions = [linear_x]
    return configure


def _matrix(path: Path) -> np.ndarray:
    values = np.loadtxt(path, delimiter=",")
    if values.ndim == 0:
        values = values.reshape((1, 1))
    elif values.ndim == 1:
        values = values.reshape((-1, 1))
    return np.asarray(values, dtype=float)


def BeforeOptim(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    result_dir = Path.cwd() / "temp" / "debug" / "result"
    integration_measures = _matrix(result_dir / "trueJArray.csv")[-1]
    internal = np.loadtxt(
        result_dir / "true_internalVirtualWork.csv", delimiter=","
    ).reshape(-1)[-1]
    external = np.loadtxt(
        result_dir / "true_externalVirtualWork.csv", delimiter=","
    ).reshape(-1)[-1]

    metrics = {
        "schema_version": 1,
        "case": "gauss_integration_solid_elements",
        "element_type": ELEMENT_TYPE,
        "virtual_scale": VIRTUAL_SCALE,
        "stress": STRESS,
        "gauss_point_count": int(integration_measures.size),
        "gauss_point_measures": [float(value) for value in integration_measures],
        "integrated_volume": float(np.sum(integration_measures)),
        "internal_virtual_work": float(internal),
        "external_virtual_work": float(external),
    }
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return configure
