"""Validate MA05P physical recovery and compare it with a reviewed baseline."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _check_close(name: str, actual: float, specification: dict) -> None:
    expected = float(specification["value"])
    absolute_tolerance = float(specification.get("atol", 0.0))
    relative_tolerance = float(specification.get("rtol", 0.0))
    if not math.isclose(
        actual, expected, abs_tol=absolute_tolerance, rel_tol=relative_tolerance
    ):
        raise AssertionError(
            f"{name}: actual={actual!r}, expected={expected!r}, "
            f"atol={absolute_tolerance!r}, rtol={relative_tolerance!r}"
        )


def _check_physics(actual: dict) -> None:
    if actual.get("case") != "ma05p_2element":
        raise AssertionError(f"Unexpected case name: {actual.get('case')!r}")

    forward = actual["forward_solve"]
    if not forward.get("required", False):
        raise AssertionError("The validation did not require a FEBio forward solve")
    if forward.get("reuse_saved_result_buffer", True):
        raise AssertionError("The validation reused an old result buffer")

    balance = actual["virtual_work_balance"]
    if not balance.get("passed", False):
        raise AssertionError(
            "FEBio/VFM virtual-work balance failed: "
            f"relative_l1={balance['relative_l1']!r}, "
            f"threshold={balance['threshold']!r}"
        )

    inverse = actual["inverse"]
    if not inverse.get("recovered", False):
        raise AssertionError(
            "VFM did not recover the known material parameters: "
            f"relative_error={inverse['relative_error']!r}, "
            f"threshold={inverse['recovery_threshold']!r}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    arguments = parser.parse_args()

    if not arguments.actual.is_file():
        raise FileNotFoundError(f"Generated metrics not found: {arguments.actual}")

    actual = _load_json(arguments.actual)
    _check_physics(actual)

    if not arguments.expected.is_file():
        print(
            "UNBASELINED: FEBio/VFM consistency and parameter recovery passed, "
            f"but no reviewed baseline exists at {arguments.expected}"
        )
        return 0

    expected = _load_json(arguments.expected)
    _check_close(
        "virtual_work_balance.relative_l1",
        float(actual["virtual_work_balance"]["relative_l1"]),
        expected["virtual_work_balance"]["relative_l1"],
    )
    for name in ("E", "gamma", "tau", "fun", "truth_fun"):
        _check_close(
            f"inverse.{name}",
            float(actual["inverse"][name]),
            expected["inverse"][name],
        )

    print("PASS: forward solve, VFM recovery, and reviewed baseline all agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
