import pytest

from crater_profile_estimation.lookup import LookupBuildConfig, build_lookup_table


@pytest.fixture(scope="session")
def lookup_path(tmp_path_factory):
    path = tmp_path_factory.mktemp("lookup") / "test_lookup.npz"
    build_lookup_table(
        path,
        LookupBuildConfig(
            diameter_min_m=1.0,
            diameter_max_m=177.0,
            diameter_points=6,
            x0_points=9,
            diffusion_points=17,
            lambda_max=0.04,
            radius_points=31,
            domain_size=41,
        ),
    )
    return path

