"""Offline lookup-table generation and loading."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .diffusion import diffuse_snapshots

DEFAULT_LOOKUP_PATH = Path("data/crater_lookup.npz")
LOOKUP_VERSION = 1


@dataclass(frozen=True)
class LookupBuildConfig:
    diameter_min_m: float = 1.0
    diameter_max_m: float = 177.0
    diameter_points: int = 16
    x0_min: float = 0.04
    x0_max: float = 0.25
    x0_points: int = 22
    lambda_max: float = 0.08
    diffusion_points: int = 81
    radius_points: int = 151
    domain_size: int = 81


@dataclass(frozen=True)
class LookupTable:
    diameter_grid_m: NDArray[np.float64]
    x0_grid: NDArray[np.float64]
    s_grid_m2: NDArray[np.float64]
    current_dd: NDArray[np.float64]
    radius_fraction: NDArray[np.float64]
    profiles_m: NDArray[np.float64]

    @classmethod
    def load(cls, path: Path | str = DEFAULT_LOOKUP_PATH) -> "LookupTable":
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(
                f"lookup table not found: {path}; run build-crater-lookup first"
            )
        with np.load(path, allow_pickle=False) as data:
            version = int(data["version"])
            if version != LOOKUP_VERSION:
                raise ValueError(
                    f"unsupported lookup version {version}; expected {LOOKUP_VERSION}"
                )
            table = cls(
                diameter_grid_m=data["diameter_grid_m"].astype(float),
                x0_grid=data["x0_grid"].astype(float),
                s_grid_m2=data["s_grid_m2"].astype(float),
                current_dd=data["current_dd"].astype(float),
                radius_fraction=data["radius_fraction"].astype(float),
                profiles_m=data["profiles_m"].astype(float),
            )
        table.validate()
        return table

    def validate(self) -> None:
        n_d = len(self.diameter_grid_m)
        n_x = len(self.x0_grid)
        n_s = self.s_grid_m2.shape[1]
        n_r = len(self.radius_fraction)
        if self.s_grid_m2.shape != (n_d, n_s):
            raise ValueError("invalid s_grid_m2 shape")
        if self.current_dd.shape != (n_d, n_x, n_s):
            raise ValueError("invalid current_dd shape")
        if self.profiles_m.shape != (n_d, n_x, n_s, n_r):
            raise ValueError("invalid profiles_m shape")
        if np.any(np.diff(self.diameter_grid_m) <= 0):
            raise ValueError("diameter grid must be strictly increasing")
        if np.any(np.diff(self.x0_grid) <= 0):
            raise ValueError("x0 grid must be strictly increasing")
        if np.any(np.diff(self.s_grid_m2, axis=1) < 0):
            raise ValueError("s grids must be non-decreasing")

    def nearest_diameter_index(self, diameter_m: float) -> int:
        if not self.diameter_grid_m[0] <= diameter_m <= self.diameter_grid_m[-1]:
            raise ValueError(
                "diameter is outside lookup range "
                f"[{self.diameter_grid_m[0]:g}, {self.diameter_grid_m[-1]:g}] m"
            )
        return int(np.argmin(np.abs(np.log(self.diameter_grid_m / diameter_m))))


def _diameter_grid(config: LookupBuildConfig) -> NDArray[np.float64]:
    if config.diameter_min_m <= 0 or config.diameter_max_m <= config.diameter_min_m:
        raise ValueError("invalid diameter range")
    if config.diameter_points < 2:
        raise ValueError("diameter_points must be at least 2")
    grid = np.geomspace(
        config.diameter_min_m, config.diameter_max_m, config.diameter_points
    )
    # Preserve model and Figure 3 boundaries exactly when they are in range.
    boundaries = np.array([10.0, 39.8, 40.0, 100.0, 177.0])
    boundaries = boundaries[
        (boundaries >= config.diameter_min_m)
        & (boundaries <= config.diameter_max_m)
    ]
    return np.unique(np.concatenate((grid, boundaries)))


def build_lookup_table(
    output: Path | str = DEFAULT_LOOKUP_PATH,
    config: LookupBuildConfig | None = None,
) -> LookupTable:
    """Run diffusion offline and save G and H lookup arrays."""
    cfg = config or LookupBuildConfig()
    diameter_grid = _diameter_grid(cfg)
    x0_grid = np.linspace(cfg.x0_min, cfg.x0_max, cfg.x0_points)
    lambda_grid = np.linspace(0.0, cfg.lambda_max, cfg.diffusion_points)
    radius_fraction = np.linspace(0.0, 1.5, cfg.radius_points)

    # Linear diffusion and the FT profile are scale invariant after D
    # normalization. Compute once at D=1, then convert to each D and s grid.
    current_unit, surfaces, grid = diffuse_snapshots(
        1.0,
        x0_grid,
        lambda_grid,
        domain_size=cfg.domain_size,
        return_surfaces=True,
    )
    assert surfaces is not None
    source_radius = grid.center_to_edge_radius
    inside = source_radius <= 1.5
    unit_profiles = np.empty(
        (len(x0_grid), len(lambda_grid), len(radius_fraction)), dtype=float
    )
    for k in range(len(lambda_grid)):
        radial = surfaces[k, :, grid.center_index, grid.center_index :][:, inside]
        for j in range(len(x0_grid)):
            unit_profiles[j, k] = np.interp(
                radius_fraction, source_radius[inside], radial[j]
            )

    s_grid = diameter_grid[:, None] ** 2 * lambda_grid[None, :]
    current_dd = np.broadcast_to(
        current_unit[None, :, :],
        (len(diameter_grid), *current_unit.shape),
    ).copy()
    profiles = (
        diameter_grid[:, None, None, None] * unit_profiles[None, :, :, :]
    )
    table = LookupTable(
        diameter_grid_m=diameter_grid,
        x0_grid=x0_grid,
        s_grid_m2=s_grid,
        current_dd=current_dd,
        radius_fraction=radius_fraction,
        profiles_m=profiles,
    )
    table.validate()

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        version=np.array(LOOKUP_VERSION),
        diameter_grid_m=diameter_grid,
        x0_grid=x0_grid,
        s_grid_m2=s_grid,
        current_dd=current_dd.astype(np.float32),
        radius_fraction=radius_fraction,
        profiles_m=profiles.astype(np.float32),
    )
    return table

