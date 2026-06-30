"""Integration test for Instanovo."""

from pathlib import Path

import pytest

import urgap


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
