import numpy as np

from crater_profile_estimation.lookup import LookupTable


def test_lookup_shapes_and_physical_constraint(lookup_path):
    table = LookupTable.load(lookup_path)
    assert table.current_dd.shape[:2] == (
        len(table.diameter_grid_m),
        len(table.x0_grid),
    )
    assert table.profiles_m.shape[:3] == table.current_dd.shape
    assert np.all(table.current_dd <= table.x0_grid[None, :, None] + 1e-7)
    assert np.all(np.diff(table.current_dd, axis=2) <= 1e-10)


def test_lookup_uses_nearest_log_diameter(lookup_path):
    table = LookupTable.load(lookup_path)
    index = table.nearest_diameter_index(20.0)
    assert 0 <= index < len(table.diameter_grid_m)
