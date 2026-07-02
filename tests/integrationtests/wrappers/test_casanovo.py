"""Integration test for Casanovo."""

import logging

from pathlib import Path
from unittest.mock import patch

import pytest

import urgap


def test_casanovo_command_construction_sequence(tmp_path: Path) -> None:
    """Test that the command_list is built correctly for sequence (de novo) mode."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.2.0": {
                    "search_mode": "sequence",
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
                f"{urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#ms_files/BSA1.mzML",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}#casanovo_yaml/casanovo.yaml",
            ),
        ],
    )

    casanovo_node = urgap.init_node("Casanovo:5.2.0")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = ""
        with pytest.raises(FileNotFoundError):
            casanovo_node.run(ufiles, urun_dict)

    actual_cmd = [str(c) for c in mock_run.call_args[0][0]]

    assert actual_cmd[0] == "casanovo"
    assert actual_cmd[1] == "sequence"
    assert actual_cmd[2].endswith("BSA1.mzML")
    assert actual_cmd[3] == "--config"
    assert actual_cmd[4].endswith("casanovo.yaml")
    assert actual_cmd[5] == "--output_dir"
    assert actual_cmd[7] == "--output_root"
    assert len(actual_cmd) == 9


def test_casanovo_command_construction_sequence_no_param_file(tmp_path: Path) -> None:
    """Test that the command_list omits --config when no param file is provided."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.2.0": {
                    "search_mode": "sequence",
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
                f"{urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#ms_files/BSA1.mzML",
            ),
        ],
    )

    casanovo_node = urgap.init_node("Casanovo:5.2.0")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = ""
        with pytest.raises(FileNotFoundError):
            casanovo_node.run(ufiles, urun_dict)

    actual_cmd = [str(c) for c in mock_run.call_args[0][0]]

    assert actual_cmd[0] == "casanovo"
    assert actual_cmd[1] == "sequence"
    assert actual_cmd[2].endswith("BSA1.mzML")
    assert "--config" not in actual_cmd
    assert actual_cmd[3] == "--output_dir"
    assert actual_cmd[5] == "--output_root"
    assert len(actual_cmd) == 7


def test_casanovo_command_construction_db_search(tmp_path: Path) -> None:
    """Test that the command_list is built correctly for db-search mode, including a FASTA."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#ms_files/BSA1.mzML",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}#casanovo_yaml/casanovo.yaml",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.2.0": {
                    "search_mode": "db-search",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )

    casanovo_node = urgap.init_node("Casanovo:5.2.0")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = ""
        with pytest.raises(FileNotFoundError):
            casanovo_node.run(ufiles, urun_dict)

    actual_cmd = [str(c) for c in mock_run.call_args[0][0]]

    assert actual_cmd[0] == "casanovo"
    assert actual_cmd[1] == "db-search"
    assert actual_cmd[2].endswith("BSA1.mzML")
    assert actual_cmd[3].endswith("BSA1.fasta")
    assert actual_cmd[4] == "--config"
    assert actual_cmd[5].endswith("casanovo.yaml")
    assert actual_cmd[6] == "--output_dir"
    assert actual_cmd[8] == "--output_root"
    assert len(actual_cmd) == 10



def test_casanovo_sequence_mode_with_fasta_logs_error(tmp_dir, caplog):
    """Providing a FASTA in sequence (de novo) mode should log an error and raise."""
    urd = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.2.0": {
                    "search_mode": "sequence",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype="
            f"{urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#ms_files/BSA1.mzML",
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype="
            f"{urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}#casanovo_yaml/casanovo.yaml",
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype="
            f"{urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta",
        ),
    )

    casanovo_node = urgap.init_node("Casanovo:5.2.0")
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError):
        casanovo_node.run(ufiles, urd)

    assert (
        "A fasta file has been provided despite the search mode being set to sequence"
        in caplog.text
    )

def test_casanovo_db_search_missing_fasta_raises(tmp_dir, caplog):
    """db-search with no FASTA should both log an error and raise ValueError."""
    urd = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.2.0": {
                    "search_mode": "db-search",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype="
            f"{urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#ms_files/BSA1.mzML",
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype="
            f"{urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}#casanovo_yaml/casanovo.yaml",
        ),
    )

    casanovo_node = urgap.init_node("Casanovo:5.2.0")
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError):
        casanovo_node.run(ufiles, urd)

    assert "Please input a Fasta file for database searching." in caplog.text
