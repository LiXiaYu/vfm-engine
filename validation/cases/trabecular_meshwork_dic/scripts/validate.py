"""Check generated validation metrics against a reviewed baseline."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _check_value(name: str, actual: float, specification: dict) -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    arguments = parser.parse_args()

    if not arguments.actual.is_file():
        raise FileNotFoundError(f"Generated metrics not found: {arguments.actual}")

    actual = _load_json(arguments.actual)
    optimization = actual["elastic_E"]
    if not optimization.get("success", False):
        raise AssertionError(
            "Elastic modulus optimization did not report success: "
            f"{optimization.get('message', 'no message')}"
        )

    if not arguments.expected.is_file():
        print(
            "UNBASELINED: execution produced valid metrics, but no reviewed "
            f"baseline exists at {arguments.expected}"
        )
        return 0

    expected = _load_json(arguments.expected)
    _check_value("elastic_E.x", float(optimization["x"]), expected["elastic_E"]["x"])
    _check_value(
        "elastic_E.fun", float(optimization["fun"]), expected["elastic_E"]["fun"]
    )
    print("PASS: generated metrics are within the reviewed tolerances")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
