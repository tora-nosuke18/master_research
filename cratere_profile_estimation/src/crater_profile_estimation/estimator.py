"""Onboard maximum-likelihood estimation using an offline lookup table."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.special import logsumexp

from .distributions import (
    Figure3Distribution,
    initial_depth_density,
)
from .lookup import DEFAULT_LOOKUP_PATH, LookupTable

SLOPE_DD_COEFFICIENT = 151.377
LOG_FLOOR = np.log(np.finfo(float).tiny)


@dataclass(frozen=True)
class EstimatorConfig:
    initial_dd_std: float = 0.02
    figure3_bandwidth: float = 0.003
    low_likelihood_threshold: float = 1.0e-6
    multimodal_relative_height: float = 0.25


@dataclass(frozen=True)
class EstimateResult:
    diameter_m: float
    lookup_diameter_m: float
    slope_deg: float | None
    sigma_slope_deg: float | None
    lambda_hat: float
    diffusion_amount_m2: float
    diffusion_q10_m2: float
    diffusion_q50_m2: float
    diffusion_q90_m2: float
    current_dd_median: float
    figure3_bin: str
    radius_m: NDArray[np.float64]
    height_q10_m: NDArray[np.float64]
    height_q50_m: NDArray[np.float64]
    height_q90_m: NDArray[np.float64]
    diffusion_values_m2: NDArray[np.float64]
    log_likelihood: NDArray[np.float64]
    likelihood: NDArray[np.float64]
    quality_flags: dict[str, bool]


def _weighted_quantile_1d(
    values: NDArray[np.float64], weights: NDArray[np.float64], q: float
) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    cumulative /= cumulative[-1]
    return float(sorted_values[np.searchsorted(cumulative, q, side="left")])


def _weighted_profile_quantile(
    values: NDArray[np.float64], weights: NDArray[np.float64], q: float
) -> NDArray[np.float64]:
    order = np.argsort(values, axis=0)
    sorted_values = np.take_along_axis(values, order, axis=0)
    expanded_weights = np.broadcast_to(weights[:, None], values.shape)
    sorted_weights = np.take_along_axis(expanded_weights, order, axis=0)
    cumulative = np.cumsum(sorted_weights, axis=0)
    cumulative /= cumulative[-1]
    indices = np.argmax(cumulative >= q, axis=0)
    return np.take_along_axis(sorted_values, indices[None, :], axis=0)[0]


def _multimodal(likelihood: NDArray[np.float64], relative_height: float) -> bool:
    if len(likelihood) < 3 or likelihood.max() <= 0:
        return False
    peaks = np.flatnonzero(
        (likelihood[1:-1] > likelihood[:-2])
        & (likelihood[1:-1] >= likelihood[2:])
        & (likelihood[1:-1] >= relative_height * likelihood.max())
    )
    return len(peaks) > 1


def estimate_profile(
    diameter_m: float,
    slope_deg: float | None = None,
    sigma_slope_deg: float | None = None,
    config: EstimatorConfig | None = None,
    lookup: LookupTable | Path | str = DEFAULT_LOOKUP_PATH,
    occluded_inner_wall: bool = False,
) -> EstimateResult:
    """Estimate a current profile without running diffusion online."""
    if diameter_m <= 0:
        raise ValueError("diameter_m must be positive")
    if (slope_deg is None) != (sigma_slope_deg is None):
        raise ValueError("slope_deg and sigma_slope_deg must be provided together")
    if sigma_slope_deg is not None and sigma_slope_deg <= 0:
        raise ValueError("sigma_slope_deg must be positive")
    cfg = config or EstimatorConfig()
    table = lookup if isinstance(lookup, LookupTable) else LookupTable.load(lookup)
    diameter_index = table.nearest_diameter_index(diameter_m)
    lookup_diameter = float(table.diameter_grid_m[diameter_index])
    x0 = table.x0_grid
    s_values = table.s_grid_m2[diameter_index]
    predicted_dd = table.current_dd[diameter_index]

    p0 = initial_depth_density(lookup_diameter, x0, cfg.initial_dd_std)
    current = Figure3Distribution.for_diameter(
        lookup_diameter, cfg.figure3_bandwidth
    )
    log_terms = np.log(np.maximum(current.pdf(predicted_dd), np.finfo(float).tiny))
    log_terms += np.log(np.maximum(p0[:, None], np.finfo(float).tiny))

    if slope_deg is not None and sigma_slope_deg is not None:
        predicted_slope = SLOPE_DD_COEFFICIENT * predicted_dd
        z = (slope_deg - predicted_slope) / sigma_slope_deg
        log_terms += -0.5 * z * z - np.log(
            sigma_slope_deg * np.sqrt(2.0 * np.pi)
        )

    # Linear diffusion may not increase d/D.
    physical = predicted_dd <= x0[:, None] + 1.0e-7
    log_terms = np.where(physical, log_terms, -np.inf)
    dx0 = np.gradient(x0)
    log_terms += np.log(dx0[:, None])
    log_likelihood = logsumexp(log_terms, axis=0)
    if not np.any(np.isfinite(log_likelihood)):
        raise RuntimeError("all lookup combinations have zero likelihood")
    best_index = int(np.nanargmax(log_likelihood))
    if best_index == len(s_values) - 1:
        raise RuntimeError(
            "maximum likelihood is at the lookup diffusion boundary; rebuild "
            "the table with a larger lambda_max"
        )

    likelihood = np.exp(log_likelihood - np.nanmax(log_likelihood))
    likelihood /= likelihood.sum()
    log_weights = log_terms[:, best_index]
    weights = np.exp(log_weights - logsumexp(log_weights))

    profiles = table.profiles_m[diameter_index, :, best_index, :]
    # Nearest-D lookup followed by scale interpolation preserves the observed D.
    scale = diameter_m / lookup_diameter
    profiles = profiles * scale
    radius_m = table.radius_fraction * diameter_m / 2.0
    dd_best = predicted_dd[:, best_index]
    s_hat = float(s_values[best_index] * scale * scale)
    s_values_observed = s_values * scale * scale

    slope_outlier = False
    if slope_deg is not None and sigma_slope_deg is not None:
        possible = SLOPE_DD_COEFFICIENT * predicted_dd[physical]
        slope_outlier = bool(
            slope_deg < possible.min() - 3.0 * sigma_slope_deg
            or slope_deg > possible.max() + 3.0 * sigma_slope_deg
        )
    quality_flags = {
        "flag_no_slope_observation": slope_deg is None,
        "flag_low_likelihood": bool(
            np.nanmax(log_likelihood) < np.log(cfg.low_likelihood_threshold)
        ),
        "flag_multimodal_likelihood": _multimodal(
            likelihood, cfg.multimodal_relative_height
        ),
        "flag_out_of_distribution_D": False,
        "flag_slope_outlier": slope_outlier,
        "flag_occluded_inner_wall": bool(occluded_inner_wall),
    }

    return EstimateResult(
        diameter_m=diameter_m,
        lookup_diameter_m=lookup_diameter,
        slope_deg=slope_deg,
        sigma_slope_deg=sigma_slope_deg,
        lambda_hat=s_hat / (diameter_m * diameter_m),
        diffusion_amount_m2=s_hat,
        diffusion_q10_m2=_weighted_quantile_1d(s_values_observed, likelihood, 0.10),
        diffusion_q50_m2=_weighted_quantile_1d(s_values_observed, likelihood, 0.50),
        diffusion_q90_m2=_weighted_quantile_1d(s_values_observed, likelihood, 0.90),
        current_dd_median=_weighted_quantile_1d(dd_best, weights, 0.50),
        figure3_bin=current.diameter_bin,
        radius_m=radius_m,
        height_q10_m=_weighted_profile_quantile(profiles, weights, 0.10),
        height_q50_m=_weighted_profile_quantile(profiles, weights, 0.50),
        height_q90_m=_weighted_profile_quantile(profiles, weights, 0.90),
        diffusion_values_m2=s_values_observed,
        log_likelihood=log_likelihood,
        likelihood=likelihood,
        quality_flags=quality_flags,
    )
