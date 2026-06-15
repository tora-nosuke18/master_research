"""Plot estimated crater profiles."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "crater-profile-matplotlib")
)
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .estimator import EstimateResult


def save_profile_plot(estimate: EstimateResult, output: Path) -> None:
    """Save a symmetric cross-section with its 10--90% uncertainty band."""
    radius = estimate.radius_m
    x = np.concatenate((-radius[:0:-1], radius))

    def symmetric(values):
        return np.concatenate((values[:0:-1], values))

    q10 = symmetric(estimate.height_q10_m)
    q50 = symmetric(estimate.height_q50_m)
    q90 = symmetric(estimate.height_q90_m)

    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    ax.fill_between(x, q10, q90, color="#4C78A8", alpha=0.25, label="10-90% interval")
    ax.plot(x, q50, color="#1F4E79", linewidth=2.2, label="Median profile")
    ax.axhline(0.0, color="0.45", linewidth=1.0, linestyle="--", label="Reference surface")
    ax.axvline(-estimate.diameter_m / 2.0, color="0.65", linewidth=0.8, linestyle=":")
    ax.axvline(estimate.diameter_m / 2.0, color="0.65", linewidth=0.8, linestyle=":")

    observation = f"D = {estimate.diameter_m:g} m"
    if estimate.slope_deg is not None:
        observation += (
            f", S = {estimate.slope_deg:g} deg"
            f" (sigma = {estimate.sigma_slope_deg:g} deg)"
        )
    ax.set_title(f"Estimated Current Crater Profile\n{observation}")
    ax.set_xlabel("Distance from crater center (m)")
    ax.set_ylabel("Elevation relative to surrounding surface (m)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right")
    ax.text(
        0.02,
        0.03,
        (
            f"Estimated d/D = {estimate.current_dd_median:.3f}\n"
            f"kappa*t = {estimate.diffusion_amount_m2:.3g} m^2"
        ),
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
