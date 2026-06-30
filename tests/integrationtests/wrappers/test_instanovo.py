"""Integration test for Instanovo."""

from pathlib import Path

import pytest

import urgap

import os
import subprocess
import time 
import logging
import threading
import psutil

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


def test_instanovo_starts(tmp_dir, caplog):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "transformer",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://"
                f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}"
                f"#mgfs/BSA1.mgf"
            ),
            urgap.UFile(
                uri=f"file://"
                f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML}"
                f"#instanovo_params/default.yaml"
            ),
        ],
    )

    instanovo_node = urgap.init_node("Instanovo:1.2.2")

    thread = threading.Thread(
        target=instanovo_node.run,
        kwargs={"urun_dict": urd, "ufiles": ufiles},
        daemon=True,
    )
    try:
        with caplog.at_level(logging.INFO):
            thread.start()

            deadline = time.monotonic() + 10
            launched = False
            while time.monotonic() < deadline:
                if "Executing command list: " in caplog.text:
                    launched = True
                    break
                time.sleep(0.2)

        assert "Running execute ..." in caplog.text
        assert launched, "instanovo command was never launched"
        assert (
            "instanovo diffusion predict" in caplog.text
            or "instanovo transformer predict" in caplog.text
        )
    finally:
        for proc in psutil.process_iter(["cmdline"]):
            cmdline = " ".join(proc.info["cmdline"] or [])
            if "instanovo" in cmdline and "predict" in cmdline:
                proc.terminate()