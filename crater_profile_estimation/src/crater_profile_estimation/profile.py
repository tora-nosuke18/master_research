"""Fresh crater profile from the synthterrain FTmod implementation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

INNER_COEFFICIENTS = (-0.228809953, 0.227533882, 0.083116795, -0.039499407)
OUTER_COEFFICIENTS = (0.188253307, -0.187050452, 0.01844746, 0.01505647)
RIM_HEIGHT_RATIO = 0.036822095


def stopar_fresh_dd(diameter_m: float) -> float:
    """Return the fresh d/D step value used by synthterrain."""
    if diameter_m < 0:
        raise ValueError("diameter_m must be non-negative")
    for lower, dd in zip(
        (400.0, 200.0, 100.0, 40.0, 10.0, 0.0),
        (0.21, 0.17, 0.15, 0.13, 0.11, 0.10),
    ):
        if diameter_m >= lower:
            return dd
    raise AssertionError("unreachable")


def _polynomial(values: NDArray[np.float64], coefficients: tuple[float, ...]):
    result = np.zeros_like(values, dtype=float)
    for coefficient in reversed(coefficients):
        result = result * values + coefficient
    return result


def fresh_profile(
    radius_fraction: ArrayLike,
    diameter_m: float,
    depth_ratio: float | None = None,
) -> NDArray[np.float64]:
    """Evaluate the modified Fassett-Thomson fresh crater profile.

    ``radius_fraction`` is distance from the center divided by crater radius.
    Elevations are returned in meters relative to the pre-existing surface.
    """
    if diameter_m <= 0:
        raise ValueError("diameter_m must be positive")
    r = np.asarray(radius_fraction, dtype=float)
    if np.any(r < 0):
        raise ValueError("radius_fraction must be non-negative")

    dd = stopar_fresh_dd(diameter_m) if depth_ratio is None else depth_ratio
    if dd <= 0:
        raise ValueError("depth_ratio must be positive")

    profile = np.zeros_like(r)
    inner = r <= 0.98
    rim = (r > 0.98) & (r <= 1.02)
    outer = (r > 1.02) & (r <= 1.5)
    profile[inner] = _polynomial(r[inner], INNER_COEFFICIENTS)
    profile[rim] = RIM_HEIGHT_RATIO
    profile[outer] = _polynomial(r[outer], OUTER_COEFFICIENTS)

    floor = RIM_HEIGHT_RATIO - dd
    profile = np.maximum(profile, floor)
    return profile * diameter_m

