"""Integration test for Instanovo."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest
import urgap

# Clean reference to prevent infinite recursion loop
_original_check_deps = urgap.unode_manager.UNodeManager.check_unode_dependencies


def mock_check_dependencies(self, unode: str) -> tuple:
    """Mock out dependency check to trick the manager into thinking the binary exists."""
    unode_obj, tmp = _original_check_deps(self, unode)
    
    # 1. Force the manager to report that the executable exists
    tmp[unode]["resource_available"] = True
    
    # 2. Use patch.object to bypass the read-only property restriction safely
    if unode_obj.exe_path is None:
        patcher = patch.object(unode_obj, "exe_path", Path("instanovo"))
        patcher.start()
        
    return unode_obj, tmp


def test_instanovo_command_construction_yaml_no_model(tmp_path: Path) -> None:
    """Test that command_list is built correctly with no model_used, using a yaml config file."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {},
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYMZML_MGF}#mgfs/BSA1.mgf",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML}#instanovo_params/default.yaml",
            ),
        ],
    )

    # Patch the manager with our recursion-safe mock wrapper
    with patch("urgap.unode_manager.UNodeManager.check_unode_dependencies", new=mock_check_dependencies):
        instanovo_node = urgap.init_node("Instanovo:1.2.2")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = ""
            with pytest.raises(FileNotFoundError):
                instanovo_node.run(ufiles, urun_dict)

    actual_cmd = [str(c) for c in mock_run.call_args[0][0]]

    assert actual_cmd[0].endswith("instanovo")
    assert actual_cmd[1] == "predict"
    assert actual_cmd[2] == "--data-path"
    assert actual_cmd[3].endswith("BSA1.mgf")
    assert actual_cmd[4] == "--output-path"
    assert actual_cmd[6] == "--config-path"
    assert actual_cmd[8] == "--config-name"
    assert actual_cmd[9] == "default"
    assert len(actual_cmd) == 10


def test_instanovo_command_construction_with_model_and_cli_params(tmp_path: Path) -> None:
    """Test that command_list includes model_used and CLI overrides when no yaml file is provided."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "transformer",
                    "num_beams": 5,
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYMZML_MGF}#mgfs/BSA1.mgf",
            ),
        ],
    )

    # Patch the manager with our recursion-safe mock wrapper
    with patch("urgap.unode_manager.UNodeManager.check_unode_dependencies", new=mock_check_dependencies):
        instanovo_node = urgap.init_node("Instanovo:1.2.2")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0
            mock_run.return_value.stdout = ""
            with pytest.raises(FileNotFoundError):
                instanovo_node.run(ufiles, urun_dict)

    actual_cmd = [str(c) for c in mock_run.call_args[0][0]]

    assert actual_cmd[0].endswith("instanovo")
    assert actual_cmd[1] == "transformer"
    assert actual_cmd[2] == "predict"
    assert actual_cmd[3] == "--data-path"
    assert actual_cmd[4].endswith("BSA1.mgf")
    assert actual_cmd[5] == "--output-path"
    assert actual_cmd[7] == "num_beams=5"
    assert len(actual_cmd) == 8


def test_instanovo_invalid_model_used_logs_error(tmp_path: Path, caplog) -> None:
    """An unrecognized model_used should log an error."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "not_a_real_mode",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYMZML_MGF}#mgfs/BSA1.mgf",
            ),
        ],
    )

    instanovo_node = urgap.init_node("Instanovo:1.2.2")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = ""
        with caplog.at_level(logging.ERROR), pytest.raises(FileNotFoundError):
            instanovo_node.run(ufiles, urun_dict)

    assert "Unknown search mode" in caplog.text


def test_instanovo_param_file_and_cli_params_raises(tmp_path: Path) -> None:
    """Providing both a yaml config file and non-model_used CLI parameters should raise ValueError."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Instanovo:1.2.2": {
                    "model_used": "transformer",
                    "num_beams": 5,
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYMZML_MGF}#mgfs/BSA1.mgf",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML}#instanovo_params/default.yaml",
            ),
        ],
    )

    instanovo_node = urgap.init_node("Instanovo:1.2.2")

    with pytest.raises(ValueError) as excinfo:
        instanovo_node.run(ufiles, urun_dict)

    assert (
        "Both a config yaml file and command-line parameters in the "
        "urun_dict were provided for Instanovo. Please provide only one."
        in str(excinfo.value)
    )