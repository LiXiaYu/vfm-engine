"""Final-cycle periodic viscoelastic identification for the MA05P case."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

import VFMTask_pybind


CASE_DIR = Path(__file__).resolve().parents[1]
METRICS_PATH = Path(
    os.environ.get("VFM_VALIDATION_METRICS", Path.cwd() / "validation_metrics.json")
).resolve()
VALIDATION_SEED = int(os.environ.get("VFM_VALIDATION_SEED", "20260820"))

PERIOD = 4.0
TRUE_PARAMETERS = np.asarray([0.3, 8.0, 1.0], dtype=float)  # E, gamma, tau
PARAMETER_BOUNDS = [(0.1, 0.9), (1.0, 9.0), (1.0, 9.0)]
MAX_PARAMETER_RELATIVE_ERROR = 0.05
MAX_VIRTUAL_WORK_RELATIVE_L1 = 1.0e-3

_selected_window: dict[str, float | int] = {}


def InitVFMTask(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    configure = VFMTask_pybind.VFMTask_configure()

    # This validation intentionally performs the FEBio forward solve first and
    # consumes the newly captured result buffer in the same process.
    configure.run_febio_solve = True
    configure.reuse_saved_result_buffer = False
    configure.allow_failed_solve_postprocessing = False
    configure.isReadfromsaveOptimfunc = False
    configure.isSetDisplacmentAndPressure = False

    configure.isLaplaceVFM = False
    configure.isSetNewLaplaceVFMs = False
    configure.isTalbotLaplaceVFM_s = False

    configure.displacementdataNumber = 0
    configure.stressdataNumber = 1
    configure.pressurevalueNumber = 0
    configure.solution = [
        ("elements", "Part5"),
        ("elements", "Part5a"),
    ]
    configure.pressure_load = [("FEPressureLoad", "Pressure1")]
    configure.fixed = [("surface", "ZeroDisplacement1")]
    # Required here because the current plugin emits the direct FEBio-stress
    # internal/external virtual-work matrices under this diagnostic switch.
    configure.optim_function_output_debug_info = True

    print("MA05P validation: FEBio forward solve is required before VFM inversion.")
    return configure


def SetPeriodTimeIndex(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    times = np.asarray(configure.get_timestep, dtype=float)
    if times.ndim != 1 or times.size < 3:
        raise RuntimeError("FEBio did not provide enough solved time states")

    end_index = int(times.size - 1)
    target_start = float(times[end_index] - PERIOD)
    start_index = int(np.argmin(np.abs(times - target_start)))
    selected_span = float(times[end_index] - times[start_index])
    time_tolerance = max(1.0e-8, 0.51 * float(np.max(np.diff(times))))
    if abs(selected_span - PERIOD) > time_tolerance:
        raise RuntimeError(
            f"Could not select the final {PERIOD:g} s cycle: span={selected_span!r}"
        )

    configure.start_index = start_index
    configure.end_index = end_index
    _selected_window.update(
        {
            "full_state_count": int(times.size),
            "start_index": start_index,
            "end_index": end_index,
            "start_time": float(times[start_index]),
            "end_time": float(times[end_index]),
            "span": selected_span,
            "selected_state_count": int(end_index - start_index + 1),
            "solved_cycle_count": int(
                round(float(times[end_index] - times[0]) / PERIOD)
            ),
            "selected_cycle_count": 1,
            "cycle_selection": "final_complete_cycle",
        }
    )
    return configure


def SetSelectSolutionElementFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    def select_solution_element_function(element_index, nodeids, node_xyzs):
        return []

    configure.select_solution_element_function = select_solution_element_function
    return configure


def SetVirtualDisplacementFunction(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    fixed_nodes = set(configure.fixednode)
    scale = 1.0e-5

    def _fixed(index_coordinate: int) -> bool:
        return int(index_coordinate) in fixed_nodes

    def vf_rational_x(coordinate, index_coordinate):
        if _fixed(index_coordinate):
            return [0.0, 0.0, 0.0]
        x = float(coordinate[0])
        return [scale * x / (x + 1.0), 0.0, 0.0]

    def vf_linear_x(coordinate, index_coordinate):
        if _fixed(index_coordinate):
            return [0.0, 0.0, 0.0]
        x = float(coordinate[0])
        return [scale * x / 0.5, 0.0, 0.0]

    def vf_constant_x(coordinate, index_coordinate):
        if _fixed(index_coordinate):
            return [0.0, 0.0, 0.0]
        return [1.0e-3, 0.0, 0.0]

    # Preserve the three virtual fields used by the original periodic
    # constitutive-identification example.
    configure.vf_u_functions = [vf_rational_x, vf_linear_x, vf_constant_x]
    return configure


def _load_matrix(path: Path) -> np.ndarray:
    matrix = np.loadtxt(path, delimiter=",")
    if matrix.ndim == 1:
        matrix = matrix.reshape((-1, 1))
    return np.asarray(matrix, dtype=float)


def _virtual_work_metrics() -> dict[str, float | int | bool]:
    result_dir = Path.cwd() / "temp" / "debug" / "result"
    external = _load_matrix(result_dir / "true_externalVirtualWork.csv")
    internal = _load_matrix(result_dir / "true_internalVirtualWork.csv")
    if external.shape != internal.shape:
        raise RuntimeError(
            f"Virtual-work shape mismatch: external={external.shape}, internal={internal.shape}"
        )

    residual = internal - external
    denominator = max(float(np.sum(np.abs(external))), np.finfo(float).tiny)
    relative_l1 = float(np.sum(np.abs(residual)) / denominator)
    max_abs = float(np.max(np.abs(residual)))
    rms = float(np.sqrt(np.mean(np.square(residual))))
    return {
        "state_count": int(external.shape[0]),
        "virtual_field_count": int(external.shape[1]),
        "relative_l1": relative_l1,
        "max_abs": max_abs,
        "rms": rms,
        "threshold": MAX_VIRTUAL_WORK_RELATIVE_L1,
        "passed": bool(relative_l1 <= MAX_VIRTUAL_WORK_RELATIVE_L1),
    }


def BeforeOptim(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    def objective_function(params) -> float:
        return float(configure.optim_function_T([float(value) for value in params]))

    truth_objective = objective_function(TRUE_PARAMETERS)
    result = differential_evolution(
        objective_function,
        PARAMETER_BOUNDS,
        seed=VALIDATION_SEED,
        workers=1,
        updating="immediate",
        popsize=8,
        maxiter=200,
        tol=1.0e-7,
        atol=1.0e-14,
        polish=True,
        disp=True,
    )

    recovered = np.asarray(result.x, dtype=float)
    relative_error = np.abs(recovered - TRUE_PARAMETERS) / np.abs(TRUE_PARAMETERS)
    virtual_work = _virtual_work_metrics()

    metrics = {
        "schema_version": 1,
        "case": "ma05p_2element",
        "seed": VALIDATION_SEED,
        "forward_solve": {
            "required": True,
            "reuse_saved_result_buffer": False,
            **_selected_window,
        },
        "truth": {
            "E": float(TRUE_PARAMETERS[0]),
            "gamma": float(TRUE_PARAMETERS[1]),
            "tau": float(TRUE_PARAMETERS[2]),
        },
        "virtual_work_balance": virtual_work,
        "inverse": {
            "E": float(recovered[0]),
            "gamma": float(recovered[1]),
            "tau": float(recovered[2]),
            "fun": float(result.fun),
            "truth_fun": float(truth_objective),
            "constitutive_update": "VFM periodic S(t) final-cycle formulation",
            "virtual_field_set": "legacy_ma05p_three_fields",
            "success": bool(result.success),
            "message": str(result.message),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
            "relative_error": {
                "E": float(relative_error[0]),
                "gamma": float(relative_error[1]),
                "tau": float(relative_error[2]),
                "maximum": float(np.max(relative_error)),
            },
            "recovery_threshold": MAX_PARAMETER_RELATIVE_ERROR,
            "recovered": bool(np.all(relative_error <= MAX_PARAMETER_RELATIVE_ERROR)),
        },
    }

    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"MA05P validation metrics written to {METRICS_PATH}")
    return configure
