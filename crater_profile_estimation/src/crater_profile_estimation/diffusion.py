"""Two-dimensional linear diffusion compatible with synthterrain."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .profile import fresh_profile


@dataclass(frozen=True)
class DiffusionGrid:
    radius_fraction: NDArray[np.float64]
    center_to_edge_radius: NDArray[np.float64]
    center_index: int


def make_grid(domain_size: int) -> DiffusionGrid:
    if domain_size < 21 or domain_size % 2 == 0:
        raise ValueError("domain_size must be an odd integer >= 21")
    axis = np.linspace(-2.0, 2.0, domain_size)
    xx, yy = np.meshgrid(axis, axis, sparse=True)
    center = domain_size // 2
    return DiffusionGrid(
        radius_fraction=np.sqrt(xx**2 + yy**2),
        center_to_edge_radius=axis[center:],
        center_index=center,
    )


def initial_surfaces(
    diameter_m: float, depth_ratios: ArrayLike, domain_size: int
) -> tuple[NDArray[np.float64], DiffusionGrid]:
    grid = make_grid(domain_size)
    ratios = np.atleast_1d(np.asarray(depth_ratios, dtype=float))
    surfaces = np.stack(
        [fresh_profile(grid.radius_fraction, diameter_m, ratio) for ratio in ratios]
    )
    return surfaces, grid


def diffuse_snapshots(
    diameter_m: float,
    depth_ratios: ArrayLike,
    lambda_values: ArrayLike,
    domain_size: int = 81,
    stability_fraction: float = 0.24,
    return_surfaces: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64] | None, DiffusionGrid]:
    """Evaluate g_D for several initial depths and diffusion amounts.

    ``lambda_values`` are the dimensionless values ``(kappa*t)/D**2``.
    The returned d/D array has shape ``(n_depth_ratios, n_lambda_values)``.
    """
    lambdas = np.asarray(lambda_values, dtype=float)
    if lambdas.ndim != 1 or len(lambdas) == 0:
        raise ValueError("lambda_values must be a non-empty one-dimensional array")
    if np.any(lambdas < 0) or np.any(np.diff(lambdas) < 0):
        raise ValueError("lambda_values must be non-negative and sorted")
    if not 0 < stability_fraction <= 0.25:
        raise ValueError("stability_fraction must be in (0, 0.25]")

    surfaces, grid = initial_surfaces(diameter_m, depth_ratios, domain_size)
    dx = 2.0 * diameter_m / domain_size
    max_step = stability_fraction * dx * dx
    targets = lambdas * diameter_m * diameter_m
    relief = np.empty((surfaces.shape[0], len(lambdas)), dtype=float)
    saved = np.empty((len(lambdas), *surfaces.shape), dtype=float) if return_surfaces else None

    current_s = 0.0
    for target_index, target_s in enumerate(targets):
        remaining = target_s - current_s
        steps = int(np.ceil(remaining / max_step)) if remaining > 0 else 0
        step_s = remaining / steps if steps else 0.0
        coefficient = step_s / (dx * dx) if steps else 0.0
        for _ in range(steps):
            old = surfaces
            updated = old.copy()
            updated[:, 1:-1, 1:-1] = old[:, 1:-1, 1:-1] + coefficient * (
                old[:, 2:, 1:-1]
                + old[:, :-2, 1:-1]
                + old[:, 1:-1, 2:]
                + old[:, 1:-1, :-2]
                - 4.0 * old[:, 1:-1, 1:-1]
            )
            surfaces = updated
        current_s = target_s
        relief[:, target_index] = (
            surfaces.max(axis=(1, 2)) - surfaces.min(axis=(1, 2))
        ) / diameter_m
        if saved is not None:
            saved[target_index] = surfaces
    return relief, saved, grid

