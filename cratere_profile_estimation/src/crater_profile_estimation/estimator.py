"""Maximum-likelihood estimation of crater degradation and profile."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .diffusion import diffuse_snapshots
from .distributions import (
    Figure3Distribution,
    gaussian_pdf,
    initial_depth_grid,
)

SLOPE_DD_COEFFICIENT = 151.377


@dataclass(frozen=True)
class EstimatorConfig:
    initial_dd_std: float = 0.02
    initial_dd_points: int = 25
    lambda_max: float = 0.08
    lambda_points: int = 81
    domain_size: int = 81
    figure3_bandwidth: float = 0.003
    profile_radius_points: int = 151


@dataclass(frozen=True)
class EstimateResult:
    diameter_m: float
    slope_deg: float | None
    sigma_slope_deg: float | None
    lambda_hat: float
    diffusion_amount_m2: float
    current_dd_median: float
    figure3_bin: str
    radius_m: NDArray[np.float64]
    height_q10_m: NDArray[np.float64]
    height_q50_m: NDArray[np.float64]
    height_q90_m: NDArray[np.float64]
    lambda_values: NDArray[np.float64]
    likelihood: NDArray[np.float64]


def _weighted_quantile(values: NDArray[np.float64], weights: NDArray[np.float64], q: float):
    order = np.argsort(values, axis=0)
    sorted_values = np.take_along_axis(values, order, axis=0)
    expanded_weights = np.broadcast_to(weights[:, None], values.shape)
    sorted_weights = np.take_along_axis(expanded_weights, order, axis=0)
    cumulative = np.cumsum(sorted_weights, axis=0)
    cumulative /= cumulative[-1]
    indices = np.argmax(cumulative >= q, axis=0)
    return np.take_along_axis(sorted_values, indices[None, :], axis=0)[0]


def estimate_profile(
    diameter_m: float,
    slope_deg: float | None = None,
    sigma_slope_deg: float | None = None,
    config: EstimatorConfig | None = None,
) -> EstimateResult:
    """Estimate degradation amount and radial profile quantiles."""
    if diameter_m <= 0:
        raise ValueError("diameter_m must be positive")
    if (slope_deg is None) != (sigma_slope_deg is None):
        raise ValueError("slope_deg and sigma_slope_deg must be provided together")
    cfg = config or EstimatorConfig()
    x0, prior_density = initial_depth_grid(
        diameter_m, cfg.initial_dd_std, cfg.initial_dd_points
    )
    lambdas = np.linspace(0.0, cfg.lambda_max, cfg.lambda_points)
    predicted_dd, _, _ = diffuse_snapshots(
        diameter_m, x0, lambdas, domain_size=cfg.domain_size
    )

    current = Figure3Distribution.for_diameter(
        diameter_m, cfg.figure3_bandwidth
    )
    integrand = current.pdf(predicted_dd) * prior_density[:, None]
    if slope_deg is not None and sigma_slope_deg is not None:
        integrand *= gaussian_pdf(
            slope_deg, SLOPE_DD_COEFFICIENT * predicted_dd, sigma_slope_deg
        )
    likelihood = np.trapezoid(integrand, x0, axis=0)
    best_index = int(np.argmax(likelihood))
    if best_index == len(lambdas) - 1:
        raise RuntimeError(
            "maximum likelihood is at lambda_max; increase EstimatorConfig.lambda_max"
        )
    lambda_hat = float(lambdas[best_index])

    _, saved, grid = diffuse_snapshots(
        diameter_m,
        x0,
        np.array([lambda_hat]),
        domain_size=cfg.domain_size,
        return_surfaces=True,
    )
    assert saved is not None
    final_surfaces = saved[0]
    posterior = integrand[:, best_index]
    posterior /= np.trapezoid(posterior, x0)
    discrete_weights = posterior / posterior.sum()

    radial_fraction = grid.center_to_edge_radius
    radial_profiles = final_surfaces[:, grid.center_index, grid.center_index :]
    inside = radial_fraction <= 1.5
    radial_fraction = radial_fraction[inside]
    radial_profiles = radial_profiles[:, inside]
    target_fraction = np.linspace(0.0, 1.5, cfg.profile_radius_points)
    interpolated = np.stack(
        [np.interp(target_fraction, radial_fraction, profile) for profile in radial_profiles]
    )

    dd_at_best = predicted_dd[:, best_index]
    return EstimateResult(
        diameter_m=diameter_m,
        slope_deg=slope_deg,
        sigma_slope_deg=sigma_slope_deg,
        lambda_hat=lambda_hat,
        diffusion_amount_m2=lambda_hat * diameter_m * diameter_m,
        current_dd_median=float(_weighted_quantile(dd_at_best[:, None], discrete_weights, 0.5)[0]),
        figure3_bin=current.diameter_bin,
        radius_m=target_fraction * diameter_m / 2.0,
        height_q10_m=_weighted_quantile(interpolated, discrete_weights, 0.10),
        height_q50_m=_weighted_quantile(interpolated, discrete_weights, 0.50),
        height_q90_m=_weighted_quantile(interpolated, discrete_weights, 0.90),
        lambda_values=lambdas,
        likelihood=likelihood,
    )
