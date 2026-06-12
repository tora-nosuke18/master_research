import numpy as np
import pytest

from crater_profile_estimation.estimator import estimate_profile
from crater_profile_estimation.lookup import LookupBuildConfig, build_lookup_table


def test_estimator_returns_finite_profile(lookup_path):
    result = estimate_profile(20.0, lookup=lookup_path)
    assert result.diffusion_amount_m2 >= 0
    assert result.figure3_bin == "10-39.8"
    assert len(result.radius_m) == 31
    assert np.all(np.isfinite(result.height_q50_m))
    assert np.isclose(result.likelihood.sum(), 1.0)
    assert result.quality_flags["flag_no_slope_observation"]


def test_steeper_observation_selects_deeper_current_crater(lookup_path):
    shallow = estimate_profile(
        20.0, slope_deg=8.0, sigma_slope_deg=1.0, lookup=lookup_path
    )
    steep = estimate_profile(
        20.0, slope_deg=14.0, sigma_slope_deg=1.0, lookup=lookup_path
    )
    assert steep.current_dd_median > shallow.current_dd_median
    assert steep.diffusion_amount_m2 <= shallow.diffusion_amount_m2


def test_estimator_rejects_truncated_lookup(tmp_path):
    path = tmp_path / "truncated.npz"
    build_lookup_table(
        path,
        LookupBuildConfig(
            diameter_points=3,
            x0_points=7,
            diffusion_points=3,
            lambda_max=0.001,
            radius_points=21,
            domain_size=31,
        ),
    )
    with pytest.raises(RuntimeError, match="diffusion boundary"):
        estimate_profile(20.0, lookup=path)


def test_estimator_does_not_run_diffusion_online(lookup_path, monkeypatch):
    import crater_profile_estimation.diffusion as diffusion

    def fail(*args, **kwargs):
        raise AssertionError("online diffusion must not run")

    monkeypatch.setattr(diffusion, "diffuse_snapshots", fail)
    result = estimate_profile(20.0, lookup=lookup_path)
    assert np.isfinite(result.current_dd_median)


def test_out_of_lookup_diameter_is_rejected(lookup_path):
    with pytest.raises(ValueError, match="outside lookup range"):
        estimate_profile(200.0, lookup=lookup_path)

