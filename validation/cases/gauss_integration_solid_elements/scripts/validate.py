"""Validate volume and virtual-work integration against analytical values."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


RTOL = 1.0e-9
ATOL = 1.0e-12


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _close(name: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=RTOL, abs_tol=ATOL):
        raise AssertionError(
            f"{name}: actual={actual:.17g}, expected={expected:.17g}, "
            f"rtol={RTOL:g}, atol={ATOL:g}"
        )


def _surface_integral(surface_nodes: list[list[float]], stress: float, scale: float) -> float:
    points = np.asarray(surface_nodes, dtype=float)
    triangles = [(0, 1, 2)] if len(points) == 3 else [(0, 1, 2), (0, 2, 3)]
    value = 0.0
    for i, j, k in triangles:
        a, b, c = points[i], points[j], points[k]
        area_vector = 0.5 * np.cross(b - a, c - a)
        area = float(np.linalg.norm(area_vector))
        normal = area_vector / area
        centroid_x = float((a[0] + b[0] + c[0]) / 3.0)
        value += stress * float(normal[0]) * scale * centroid_x * area
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    arguments = parser.parse_args()

    actual = _load(arguments.actual)
    expected_cases = _load(arguments.expected)
    element_type = actual["element_type"]
    expected = expected_cases[element_type]

    gauss_points = int(actual["gauss_point_count"])
    if gauss_points != int(expected["gauss_points"]):
        raise AssertionError(
            f"{element_type}.gauss_points: actual={gauss_points}, "
            f"expected={expected['gauss_points']}"
        )

    volume = float(expected["volume"])
    stress = float(actual["stress"])
    scale = float(actual["virtual_scale"])
    green_strain_xx = scale + 0.5 * scale * scale
    expected_internal = stress * green_strain_xx * volume
    expected_external = _surface_integral(
        expected["surface_nodes"], stress=stress, scale=scale
    )

    _close(f"{element_type}.integrated_volume", float(actual["integrated_volume"]), volume)
    _close(
        f"{element_type}.internal_virtual_work",
        float(actual["internal_virtual_work"]),
        expected_internal,
    )
    _close(
        f"{element_type}.external_virtual_work",
        float(actual["external_virtual_work"]),
        expected_external,
    )
    print(f"PASS: {element_type} volume and virtual-work integration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
