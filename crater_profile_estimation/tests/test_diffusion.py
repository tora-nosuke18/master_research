import numpy as np

from crater_profile_estimation.diffusion import diffuse_snapshots


def test_diffusion_reduces_relief_monotonically():
    values, _, _ = diffuse_snapshots(
        20.0, [0.11, 0.13], np.linspace(0, 0.03, 7), domain_size=41
    )
    assert np.all(np.diff(values, axis=1) <= 1e-12)


def test_dimensionless_diffusion_is_scale_invariant():
    lambdas = np.array([0.0, 0.01, 0.02])
    small, _, _ = diffuse_snapshots(10.0, [0.12], lambdas, domain_size=41)
    large, _, _ = diffuse_snapshots(100.0, [0.12], lambdas, domain_size=41)
    np.testing.assert_allclose(small, large, atol=1e-12)

