import json

from crater_profile_estimation.cli import main


def test_cli_writes_default_results_to_output(tmp_path, monkeypatch, lookup_path):
    monkeypatch.chdir(tmp_path)
    exit_code = main(
        [
            "--diameter",
            "20",
            "--lookup",
            str(lookup_path),
        ]
    )

    profile = tmp_path / "output" / "profile.csv"
    plot = tmp_path / "output" / "profile.png"
    likelihood = tmp_path / "output" / "likelihood.csv"
    summary = tmp_path / "output" / "summary.json"
    assert exit_code == 0
    assert profile.is_file()
    assert plot.is_file()
    assert likelihood.is_file()
    assert summary.is_file()
    assert plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert profile.read_text(encoding="utf-8").startswith("radius_m,height_q10_m")
    summary_data = json.loads(summary.read_text(encoding="utf-8"))
    assert summary_data["profile_csv"] == "output/profile.csv"
    assert "quality_flags" in summary_data
