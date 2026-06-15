import numpy as np
import pytest

from crater_profile_estimation.profile import (
    RIM_HEIGHT_RATIO,
    fresh_profile,
    stopar_fresh_dd,
)


@pytest.mark.parametrize(
    ("diameter", "expected"),
    [(5, 0.10), (10, 0.11), (40, 0.13), (100, 0.15), (200, 0.17), (400, 0.21)],
)
def test_stopar_steps(diameter, expected):
    assert stopar_fresh_dd(diameter) == expected


def test_profile_depth_is_measured_from_rim():
    diameter = 20.0
    dd = 0.12
    heights = fresh_profile(np.array([0.0, 1.0]), diameter, dd)
    assert heights[1] == pytest.approx(RIM_HEIGHT_RATIO * diameter)
    assert (heights[1] - heights[0]) / diameter == pytest.approx(dd)

