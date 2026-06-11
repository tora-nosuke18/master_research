import numpy as np
import pytest

from crater_profile_estimation.distributions import Figure3Distribution, figure3_bin


def test_figure3_bins():
    assert figure3_bin(1) == "1-10"
    assert figure3_bin(10) == "10-39.8"
    assert figure3_bin(39.8) == "39.8-177"
    with pytest.raises(ValueError):
        figure3_bin(200)


def test_digitized_distribution_is_ordered_by_size():
    medians = [
        np.median(Figure3Distribution.for_diameter(d).samples)
        for d in (5, 20, 100)
    ]
    assert medians[0] < medians[1] < medians[2]

