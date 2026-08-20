"""Forward-solve and inverse-identification setup for the MA05P two-element case."""

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
PARAMETER_BOUNDS = [(0.1, 0.7), (2.0, 12.0), (0.25, 3.0)]
MAX_PARAMETER_RELATIVE_ERROR = 0.05
MAX_VIRTUAL_WORK_RELATIVE_L1 = 1.0e-3

_selected_window: dict[str, float | int] = {}
_selected_times = np.empty(0, dtype=float)


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
    global _selected_times

    times = np.asarray(configure.get_timestep, dtype=float)
    if times.ndim != 1 or times.size < 3:
        raise RuntimeError("FEBio did not provide enough solved time states")

    # Keep the original model configuration as the kinematic reference.  The
    # current plugin rebases displacements at start_index, so selecting only a
    # late cycle would discard the finite-strain prestress/history needed by
    # the FEBio constitutive law.
    start_index = 0
    end_index = int(times.size - 1)
    selected_span = float(times[end_index] - times[start_index])

    configure.start_index = start_index
    configure.end_index = end_index
    _selected_times = times[start_index : end_index + 1].copy()
    _selected_window.update(
        {
            "full_state_count": int(times.size),
            "start_index": start_index,
            "end_index": end_index,
            "start_time": float(times[start_index]),
            "end_time": float(times[end_index]),
            "span": selected_span,
            "selected_state_count": int(end_index - start_index + 1),
            "loading_cycle_count": float(selected_span / PERIOD),
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
    scale = 1.0e-4

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
        return [scale * (x + 0.5), 0.0, 0.0]

    def vf_target_element_x(coordinate, index_coordinate):
        if _fixed(index_coordinate):
            return [0.0, 0.0, 0.0]
        # Equal values at the interface and loaded face produce zero virtual
        # strain in the known elastic support and isolate the target element.
        return [scale, 0.0, 0.0]

    configure.vf_u_functions = [
        vf_rational_x,
        vf_linear_x,
        vf_target_element_x,
    ]
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


def _build_transient_objective():
    """Build the one-term Prony objective with FEBio 4.13's recurrence."""
    result_dir = Path.cwd() / "temp" / "debug" / "result"
    unit_pk2_by_element = _load_matrix(result_dir / "S_e_0.csv")
    jacobians = _load_matrix(result_dir / "trueJArray.csv")
    virtual_strain = _load_matrix(result_dir / "virtualstrainArrayV.csv")
    target_external = _load_matrix(result_dir / "T_exter_nEvw.csv")[:, 0]
    visco_mask = np.loadtxt(
        result_dir / "ElementInSolution_visco_mask.csv", delimiter=",", dtype=int
    ).reshape(-1)

    state_count = int(_selected_times.size)
    element_count = int(unit_pk2_by_element.shape[0])
    if state_count < 2 or element_count < 1:
        raise RuntimeError("The transient constitutive history is empty")
    if unit_pk2_by_element.shape[1] != state_count * 6:
        raise RuntimeError("Unexpected S_e_0.csv shape")
    if visco_mask.size != element_count or int(visco_mask[0]) != 1:
        raise RuntimeError("The first solution element must be the viscoelastic target")
    if jacobians.shape[0] != state_count:
        raise RuntimeError("Unexpected trueJArray.csv state count")
    if jacobians.shape[1] % element_count != 0:
        raise RuntimeError("Cannot infer the number of Gauss points")
    if virtual_strain.shape != (state_count, element_count * 6):
        raise RuntimeError("Unexpected virtualstrainArrayV.csv shape")
    if target_external.size != state_count:
        raise RuntimeError("Unexpected T_exter_nEvw.csv state count")

    gauss_points = int(jacobians.shape[1] // element_count)
    unit_pk2 = unit_pk2_by_element[0].reshape((state_count, 6))
    target_jacobian_sum = np.sum(jacobians[:, :gauss_points], axis=1)
    target_virtual_strain = virtual_strain[:, :6]
    symmetric_dotdot_weights = np.asarray([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
    time_steps = np.diff(_selected_times)
    if np.any(time_steps <= 0.0):
        raise RuntimeError("The selected FEBio time states are not strictly increasing")

    def objective_function(params) -> float:
        elastic_modulus, gamma, tau = [float(value) for value in params]
        if tau <= 0.0:
            return float("inf")

        history = np.zeros(6, dtype=float)
        pk2 = np.empty_like(unit_pk2)
        pk2[0] = elastic_modulus * unit_pk2[0]
        for state_index, dt in enumerate(time_steps, start=1):
            decay = float(np.exp(-dt / tau))
            ramp = float((1.0 - decay) / (dt / tau))
            history = (
                decay * history
                + ramp * (unit_pk2[state_index] - unit_pk2[state_index - 1])
            )
            pk2[state_index] = elastic_modulus * (
                unit_pk2[state_index] + gamma * history
            )

        target_internal = target_jacobian_sum * np.sum(
            pk2 * target_virtual_strain * symmetric_dotdot_weights, axis=1
        )
        return float(np.sum(np.abs(target_internal - target_external)))

    diagnostics = {
        "constitutive_update": "FEBio 4.13 classic viscoelastic discrete recurrence",
        "history_initialization": "zero at model time zero",
        "virtual_field_index": 0,
        "target_element_index": 0,
        "gauss_point_count": gauss_points,
    }
    return objective_function, diagnostics


def BeforeOptim(
    configure: VFMTask_pybind.VFMTask_configure,
) -> VFMTask_pybind.VFMTask_configure:
    objective_function, objective_diagnostics = _build_transient_objective()

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
            **objective_diagnostics,
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
