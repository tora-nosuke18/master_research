"""Initial and current d/D distributions used by the estimator."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import ndtr

from .profile import stopar_fresh_dd

# Quantiles digitized from Fassett et al. (2022), Figure 3 at 0.05 CDF
# intervals. The original figure has no tabulated numerical data.
FIGURE3_QUANTILES = np.arange(0.05, 1.0, 0.05)
FIGURE3_DD = {
    "1-10": np.array(
        [
            0.03946, 0.04589, 0.05015, 0.05341, 0.05632, 0.05898,
            0.06152, 0.06429, 0.06702, 0.06959, 0.07221, 0.07476,
            0.07764, 0.08083, 0.08471, 0.08865, 0.09301, 0.09862,
            0.10736,
        ]
    ),
    "10-39.8": np.array(
        [
            0.04413, 0.04973, 0.05387, 0.05704, 0.06011, 0.06285,
            0.06552, 0.06837, 0.07112, 0.07376, 0.07649, 0.07933,
            0.08286, 0.08647, 0.08989, 0.09381, 0.09918, 0.10526,
            0.11518,
        ]
    ),
    "39.8-177": np.array(
        [
            0.04984, 0.05475, 0.05820, 0.06096, 0.06396, 0.06679,
            0.06947, 0.07265, 0.07604, 0.07955, 0.08308, 0.08726,
            0.09168, 0.09611, 0.10057, 0.10573, 0.11179, 0.11970,
            0.13070,
        ]
    ),
}


def figure3_bin(diameter_m: float) -> str:
    if not 1.0 <= diameter_m <= 177.0:
        raise ValueError("Fassett Figure 3 supports diameters from 1 to 177 m")
    if diameter_m < 10.0:
        return "1-10"
    if diameter_m < 39.8:
        return "10-39.8"
    return "39.8-177"


@dataclass(frozen=True)
class Figure3Distribution:
    diameter_bin: str
    bandwidth: float = 0.003

    @classmethod
    def for_diameter(cls, diameter_m: float, bandwidth: float = 0.003):
        return cls(figure3_bin(diameter_m), bandwidth)

    @property
    def samples(self) -> NDArray[np.float64]:
        # Equal-CDF quantiles form an equal-weight empirical approximation.
        return FIGURE3_DD[self.diameter_bin]

    def pdf(self, values: ArrayLike) -> NDArray[np.float64]:
        x = np.asarray(values, dtype=float)[..., None]
        z = (x - self.samples) / self.bandwidth
        kernels = np.exp(-0.5 * z * z) / (self.bandwidth * math.sqrt(2.0 * math.pi))
        return kernels.mean(axis=-1)

    def cdf(self, values: ArrayLike) -> NDArray[np.float64]:
        x = np.asarray(values, dtype=float)[..., None]
        return ndtr((x - self.samples) / self.bandwidth).mean(axis=-1)


def initial_depth_grid(
    diameter_m: float, standard_deviation: float, points: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if standard_deviation <= 0:
        raise ValueError("standard_deviation must be positive")
    if points < 5:
        raise ValueError("points must be at least 5")
    mean = stopar_fresh_dd(diameter_m)
    lower = max(0.04, mean - 4.0 * standard_deviation)
    upper = min(0.30, mean + 4.0 * standard_deviation)
    values = np.linspace(lower, upper, points)
    weights = np.exp(-0.5 * ((values - mean) / standard_deviation) ** 2)
    weights /= np.trapezoid(weights, values)
    return values, weights


def initial_depth_density(
    diameter_m: float,
    values: ArrayLike,
    standard_deviation: float,
) -> NDArray[np.float64]:
    """Evaluate and normalize p0 on an existing lookup-table x0 grid."""
    if standard_deviation <= 0:
        raise ValueError("standard_deviation must be positive")
    x0 = np.asarray(values, dtype=float)
    if x0.ndim != 1 or len(x0) < 2 or np.any(np.diff(x0) <= 0):
        raise ValueError("values must be a strictly increasing 1-D grid")
    mean = stopar_fresh_dd(diameter_m)
    density = np.exp(-0.5 * ((x0 - mean) / standard_deviation) ** 2)
    integral = np.trapezoid(density, x0)
    if integral <= 0:
        raise ValueError("initial d/D distribution has no support on lookup grid")
    return density / integral


def gaussian_pdf(value: float, means: ArrayLike, sigma: float) -> NDArray[np.float64]:
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    z = (value - np.asarray(means, dtype=float)) / sigma
    return np.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi))
