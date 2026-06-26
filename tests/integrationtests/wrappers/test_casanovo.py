"""Integration test for Casanovo."""

from pathlib import Path

import pytest

import urgap


def test_casanovo_node_init() -> None:
    """Test that the Casanovo node can be initialized."""
    casanovo_node = urgap.init_node("Casanovo:5.1.2")
    assert casanovo_node is not None


def test_casanovo_ufile_construction_sequence(tmp_path: Path) -> None:
    """Test that UFiles can be constructed for sequence (de novo) mode."""
    mzml_file = tmp_path / "test.mzML"
    yaml_file = tmp_path / "default.yaml"
    mzml_file.touch()
    yaml_file.touch()

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{mzml_file.parent}?uftype={urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}"
                f"#{mzml_file.name}",
            ),
            urgap.UFile(
                uri=f"file://{yaml_file.parent}?uftype={urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}"
                f"#{yaml_file.name}",
            ),
        ],
    )

    assert len(ufiles) == 2
    mzml_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML)
    yaml_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML)
    assert len(mzml_ufiles) == 1
    assert len(yaml_ufiles) == 1


def test_casanovo_ufile_construction_db_search(tmp_path: Path) -> None:
    """Test that UFiles can be constructed for db-search mode, including a FASTA."""
    mzml_file = tmp_path / "test.mzML"
    yaml_file = tmp_path / "default.yaml"
    fasta_file = tmp_path / "test.fasta"
    mzml_file.touch()
    yaml_file.touch()
    fasta_file.touch()

    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{mzml_file.parent}?uftype={urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}"
                f"#{mzml_file.name}",
            ),
            urgap.UFile(
                uri=f"file://{yaml_file.parent}?uftype={urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML}"
                f"#{yaml_file.name}",
            ),
            urgap.UFile(
                uri=f"file://{fasta_file.parent}?uftype={urgap.uftypes.proteomics.FASTA}"
                f"#{fasta_file.name}",
            ),
        ],
    )

    assert len(ufiles) == 3
    mzml_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML)
    yaml_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML)
    fasta_ufiles = ufiles.get_path_objects_by_uftype(urgap.uftypes.proteomics.FASTA)
    assert len(mzml_ufiles) == 1
    assert len(yaml_ufiles) == 1
    assert len(fasta_ufiles) == 1


def test_casanovo_urun_dict_sequence() -> None:
    """Test that a URunDict with sequence mode can be constructed."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.1.2": {
                    "search_mode": "sequence",
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Casanovo:5.1.2"]
    assert params["search_mode"] == "sequence"


def test_casanovo_urun_dict_db_search() -> None:
    """Test that a URunDict with db-search mode can be constructed."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.1.2": {
                    "search_mode": "db-search",
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Casanovo:5.1.2"]
    assert params["search_mode"] == "db-search"


def test_casanovo_invalid_search_mode() -> None:
    """Test that an invalid search_mode value is not one of the accepted modes."""
    invalid_mode = "invalid_mode"
    valid_modes = ["sequence", "db-search"]
    assert invalid_mode not in valid_modes


@pytest.mark.parametrize("search_mode", ["sequence", "db-search"])
def test_casanovo_urun_dict_search_mode(search_mode: str) -> None:
    """Test that search_mode is correctly stored in URunDict parameters."""
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Casanovo:5.1.2": {
                    "search_mode": search_mode,
                },
            },
            "unode_parameters": {
                "storage_base_uri": "file:///tmp",
            },
        },
    )
    params = urun_dict.parameters["Casanovo:5.1.2"]
    assert params["search_mode"] == search_mode
    assert params["search_mode"] in ["sequence", "db-search"]
