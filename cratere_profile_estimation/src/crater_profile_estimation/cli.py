"""Command-line interface for crater profile estimation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .estimator import EstimatorConfig, estimate_profile
from .plotting import save_profile_plot

DEFAULT_OUTPUT = Path("output/profile.csv")
DEFAULT_SUMMARY = Path("output/summary.json")
DEFAULT_PLOT = Path("output/profile.png")


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
    result.add_argument("--domain-size", type=int, default=81)
    result.add_argument("--lambda-max", type=float, default=0.08)
    result.add_argument("--lambda-points", type=int, default=81)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if (args.slope is None) != (args.sigma_slope is None):
        parser().error("--slope and --sigma-slope must be specified together")
    config = EstimatorConfig(
        domain_size=args.domain_size,
        lambda_max=args.lambda_max,
        lambda_points=args.lambda_points,
    )
    estimate = estimate_profile(
        args.diameter, args.slope, args.sigma_slope, config=config
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["radius_m", "height_q10_m", "height_q50_m", "height_q90_m"])
        writer.writerows(
            zip(
                estimate.radius_m,
                estimate.height_q10_m,
                estimate.height_q50_m,
                estimate.height_q90_m,
            )
        )
    save_profile_plot(estimate, args.plot)
    summary = {
        "diameter_m": estimate.diameter_m,
        "slope_deg": estimate.slope_deg,
        "sigma_slope_deg": estimate.sigma_slope_deg,
        "figure3_bin": estimate.figure3_bin,
        "lambda_hat": estimate.lambda_hat,
        "diffusion_amount_m2": estimate.diffusion_amount_m2,
        "current_dd_median": estimate.current_dd_median,
        "profile_csv": str(args.output),
        "profile_plot": str(args.plot),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
