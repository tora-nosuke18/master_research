"""Command-line interface for crater profile estimation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .estimator import EstimatorConfig, estimate_profile
from .lookup import DEFAULT_LOOKUP_PATH
from .plotting import save_profile_plot

DEFAULT_OUTPUT = Path("output/profile.csv")
DEFAULT_SUMMARY = Path("output/summary.json")
DEFAULT_PLOT = Path("output/profile.png")
DEFAULT_LIKELIHOOD = Path("output/likelihood.csv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--diameter", type=float, required=True, help="crater diameter in m")
    result.add_argument("--slope", type=float, help="observed inner-wall slope in degrees")
    result.add_argument("--sigma-slope", type=float, help="1-sigma slope uncertainty in degrees")
    result.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output profile CSV (default: {DEFAULT_OUTPUT})",
    )
    result.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"JSON summary path (default: {DEFAULT_SUMMARY})",
    )
    result.add_argument(
        "--plot",
        type=Path,
        default=DEFAULT_PLOT,
        help=f"profile image path (default: {DEFAULT_PLOT})",
    )
    result.add_argument(
        "--likelihood-output",
        type=Path,
        default=DEFAULT_LIKELIHOOD,
        help=f"likelihood CSV path (default: {DEFAULT_LIKELIHOOD})",
    )
    result.add_argument(
        "--lookup",
        type=Path,
        default=DEFAULT_LOOKUP_PATH,
        help=f"offline lookup table (default: {DEFAULT_LOOKUP_PATH})",
    )
    result.add_argument(
        "--occluded-inner-wall",
        action="store_true",
        help="set the inner-wall occlusion quality flag",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if (args.slope is None) != (args.sigma_slope is None):
        parser().error("--slope and --sigma-slope must be specified together")
    config = EstimatorConfig()
    estimate = estimate_profile(
        args.diameter,
        args.slope,
        args.sigma_slope,
        config=config,
        lookup=args.lookup,
        occluded_inner_wall=args.occluded_inner_wall,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["radius_m", "height_q10_m", "height_q50_m", "height_q90_m"])
        writer.writerows(
            zip(
                estimate.radius_m,
                estimate.height_q10_m,
                estimate.height_q50_m,
                estimate.height_q90_m,
            )
        )
    args.likelihood_output.parent.mkdir(parents=True, exist_ok=True)
    with args.likelihood_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["diffusion_amount_m2", "log_likelihood", "likelihood"])
        writer.writerows(
            zip(
                estimate.diffusion_values_m2,
                estimate.log_likelihood,
                estimate.likelihood,
            )
        )
    save_profile_plot(estimate, args.plot)
    summary = {
        "diameter_m": estimate.diameter_m,
        "lookup_diameter_m": estimate.lookup_diameter_m,
        "slope_deg": estimate.slope_deg,
        "sigma_slope_deg": estimate.sigma_slope_deg,
        "figure3_bin": estimate.figure3_bin,
        "lambda_hat": estimate.lambda_hat,
        "diffusion_amount_m2": estimate.diffusion_amount_m2,
        "diffusion_q10_m2": estimate.diffusion_q10_m2,
        "diffusion_q50_m2": estimate.diffusion_q50_m2,
        "diffusion_q90_m2": estimate.diffusion_q90_m2,
        "current_dd_median": estimate.current_dd_median,
        "quality_flags": estimate.quality_flags,
        "lookup_table": str(args.lookup),
        "profile_csv": str(args.output),
        "profile_plot": str(args.plot),
        "likelihood_csv": str(args.likelihood_output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
