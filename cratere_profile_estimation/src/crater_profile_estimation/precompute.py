"""Build the offline crater diffusion lookup table."""

from __future__ import annotations

import argparse
from pathlib import Path

from .lookup import DEFAULT_LOOKUP_PATH, LookupBuildConfig, build_lookup_table


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--output", type=Path, default=DEFAULT_LOOKUP_PATH)
    result.add_argument("--diameter-min", type=float, default=1.0)
    result.add_argument("--diameter-max", type=float, default=177.0)
    result.add_argument("--diameter-points", type=int, default=16)
    result.add_argument("--x0-points", type=int, default=22)
    result.add_argument("--diffusion-points", type=int, default=81)
    result.add_argument("--lambda-max", type=float, default=0.08)
    result.add_argument("--radius-points", type=int, default=151)
    result.add_argument("--domain-size", type=int, default=81)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    config = LookupBuildConfig(
        diameter_min_m=args.diameter_min,
        diameter_max_m=args.diameter_max,
        diameter_points=args.diameter_points,
        x0_points=args.x0_points,
        diffusion_points=args.diffusion_points,
        lambda_max=args.lambda_max,
        radius_points=args.radius_points,
        domain_size=args.domain_size,
    )
    table = build_lookup_table(args.output, config)
    print(
        f"wrote {args.output}: D={len(table.diameter_grid_m)}, "
        f"x0={len(table.x0_grid)}, s={table.s_grid_m2.shape[1]}, "
        f"r={len(table.radius_fraction)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

