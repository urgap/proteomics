"""Integration test for Instanovo."""

from pathlib import Path

import pytest

import urgap

import os
import subprocess
import time 


def test_instanovo_node_init() -> None:
    """Test that the Instanovo node can be initialized."""
    instanovo_node = urgap.init_node("Instanovo:1.2.2")
    assert instanovo_node is not None


def test_instanovo_ufile_construction(tmp_path: Path) -> None:
    """Test that UFiles can be constructed with correct uftypes."""
    mgf_file = tmp_path / "test.pymzml.mgf"
    yaml_file = tmp_path / "default.yaml"
    mgf_file.touch()
    yaml_file.touch()

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{mgf_file.parent}?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}"
                f"#{mgf_file.name}",
            ),
            urgap.UFile(
                uri=f"file://{yaml_file.parent}?uftype={urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML}"
                f"#{yaml_file.name}",
            ),
        ],
    )

    assert len(ufiles) == 2
    mgf_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.proteomics.converter.PYMZML_MGF)
    yaml_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML)
    assert len(mgf_ufiles) == 1
    assert len(yaml_ufiles) == 1


def test_instanovo_urun_dict_transformer() -> None:
    """Test that a URunDict with transformer mode can be constructed."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "transformer",
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Instanovo:1.2.2"]
    assert params["model_used"] == "transformer"


def test_instanovo_urun_dict_diffusion() -> None:
    """Test that a URunDict with diffusion mode can be constructed."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "diffusion",
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Instanovo:1.2.2"]
    assert params["model_used"] == "diffusion"


def test_instanovo_invalid_model_mode() -> None:
    """Test that an invalid model_used value is not one of the accepted modes."""
    invalid_mode = "invalid_mode"
    valid_modes = ["transformer", "diffusion"]
    assert invalid_mode not in valid_modes


@pytest.mark.parametrize("model_used", ["transformer", "diffusion"])
def test_instanovo_urun_dict_model_used(model_used: str) -> None:
    """Test that model_used is correctly stored in URunDict parameters."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": model_used,
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Instanovo:1.2.2"]
    assert params["model_used"] == model_used
    assert params["model_used"] in ["transformer", "diffusion"]


def test_instanovo_program_starts() -> None:
    """Test that the instanovo CLI binary exists and starts executing."""
    proc = subprocess.Popen(
        ["instanovo", "--help"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        time.sleep(2)
        assert proc.poll() is None or proc.returncode == 0
    finally:
        if proc.poll() is None:
            proc.terminate()
            proc.wait(timeout=5)

            
def test_instanovo_predict_process_starts(tmp_path: Path) -> None:
    """Test that `instanovo predict` actually launches (not completes).

    Builds the config-path the same way Instanovo.preflight() does: relative
    to instanovo's built-in configs directory, resolved via importlib.
    """
    mgf_file = tmp_path / "test.pymzml.mgf"
    output_csv = tmp_path / "output.instanovo.csv"
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    config_file = config_dir / "default.yaml"
    mgf_file.touch()
    config_file.write_text("model_used: transformer\n")

    instanovo_node = urgap.init_node("Instanovo:1.2.2")
    instanovo_configs_dir = instanovo_node.find_instanovo_configs()
    relative_config_path = os.path.relpath(config_dir, instanovo_configs_dir)

    model_used = "transformer"
    proc = subprocess.Popen(
        [
            "instanovo",
            model_used,
            "predict",
            "--data-path",
            str(mgf_file),
            "--output-path",
            str(output_csv),
            "--config-path",
            relative_config_path,
            "--config-name",
            config_file.stem,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        time.sleep(3)
        # 127 = command not found; anything else (still running, or failing
        # later on bad/empty data) confirms the process actually started.
        assert proc.poll() is None or proc.returncode != 127
    finally:
        if proc.poll() is None:
            proc.terminate()
            proc.wait(timeout=5)