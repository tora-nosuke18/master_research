import numpy as np
import pytest

from crater_profile_estimation.comparison import measure_inner_wall_slope


def test_measure_inner_wall_slope_from_symmetric_profile():
    x = np.linspace(-5.0, 5.0, 101)
    height = 0.2 * np.abs(x)
    slope, mask = measure_inner_wall_slope(x, height, diameter_m=10.0)
    assert slope == pytest.approx(np.degrees(np.arctan(0.2)))
    assert mask.any()
