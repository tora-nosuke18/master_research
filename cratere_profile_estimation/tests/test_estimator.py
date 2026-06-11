import numpy as np
import pytest

from crater_profile_estimation.estimator import EstimatorConfig, estimate_profile


FAST = EstimatorConfig(
    initial_dd_points=9,
    lambda_max=0.04,
    lambda_points=17,
    domain_size=41,
    profile_radius_points=31,
)


def test_estimator_returns_finite_profile():
    result = estimate_profile(20.0, config=FAST)
    assert 0 <= result.lambda_hat <= FAST.lambda_max
    assert result.figure3_bin == "10-39.8"
    assert len(result.radius_m) == FAST.profile_radius_points
    assert np.all(np.isfinite(result.height_q50_m))


def test_steeper_observation_selects_deeper_current_crater():
    shallow = estimate_profile(20.0, slope_deg=8.0, sigma_slope_deg=1.0, config=FAST)
    steep = estimate_profile(20.0, slope_deg=14.0, sigma_slope_deg=1.0, config=FAST)
    assert steep.current_dd_median > shallow.current_dd_median
    assert steep.lambda_hat <= shallow.lambda_hat


def test_estimator_rejects_truncated_lambda_search():
    truncated = EstimatorConfig(
        initial_dd_points=9,
        lambda_max=0.001,
        lambda_points=3,
        domain_size=41,
    )
    with pytest.raises(RuntimeError, match="lambda_max"):
        estimate_profile(20.0, config=truncated)

