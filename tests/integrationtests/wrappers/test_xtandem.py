"""Integration test for Xtandem_Alanine."""
from pathlib import Path
from xml.etree import ElementTree

import urgap

import pytest


def test_searchdb_xtandem_alanine_default_params(tmp_path: Path) -> None:
    mgf = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/mgfs?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}#BSA1.mgf")
    fasta = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta")
    params = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/xtandem_params?uftype={urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS}#xtandem_params.xml")
    urd = urgap.URunDict(
        {
            "parameters": {"Xtandem_Alanine:1.0.0": {}},
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )
    search_node = urgap.init_node("Xtandem_Alanine:1.0.0")
    search_result = search_node.run([mgf, params, fasta], urd)
    assert search_result[0].path.exists() is True

    # check that output is not empty
    tree = ElementTree.parse(search_result[0].path)
    model_groups = tree.getroot().findall('group[@type="model"]')
    assert len(model_groups) > 0

    sequences = {
        domain.attrib["seq"]
        for group in model_groups
        for domain in group.iter("domain")
    }
    assert "ECCDKPLLEK" in sequences

  #  search_node.remove_output_folder(output_file=search_result[0])

def test_searchdb_xtandem_alanine__non_default_params(tmp_path: Path) -> None:
    mgf = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/mgfs?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}#BSA1.mgf")
    fasta = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta")
    params = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/xtandem_params?uftype={urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS}#xtandem_params_non_default.xml")
    urd = urgap.URunDict(
        {
            "parameters": {"Xtandem_Alanine:1.0.0": {}},
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )
    search_node = urgap.init_node("Xtandem_Alanine:1.0.0")
    search_result = search_node.run([mgf, params, fasta], urd)
    assert search_result[0].path.exists() is True

    # check that output is not empty
    tree = ElementTree.parse(search_result[0].path)
    model_groups = tree.getroot().findall('group[@type="model"]')
    assert len(model_groups) > 0

    sequences = {
        domain.attrib["seq"]
        for group in model_groups
        for domain in group.iter("domain")
    }
    assert "ECCDKPLLEK" not in sequences



def test_searchdb_xtandem_alanine_no_param_file_raises(tmp_path: Path) -> None:
    """Omitting the param file should raise ValueError."""
    mgf = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/mgfs?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}#BSA1.mgf")
    fasta = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta")
    urd = urgap.URunDict(
        {
            "parameters": {"Xtandem_Alanine:1.0.0": {}},
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )
    search_node = urgap.init_node("Xtandem_Alanine:1.0.0")
    with pytest.raises(ValueError):
        search_node.run([mgf, fasta], urd)


def test_searchdb_xtandem_alanine_multiple_param_files_raises(tmp_path: Path) -> None:
    """Supplying more than one param file should raise ValueError."""
    mgf = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/mgfs?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}#BSA1.mgf")
    fasta = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta")
    params_1 = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/xtandem_params?uftype={urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS}#xtandem_params.xml")
    params_2 = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/xtandem_params?uftype={urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS}#xtandem_params_2.xml")
    urd = urgap.URunDict(
        {
            "parameters": {"Xtandem_Alanine:1.0.0": {}},
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )
    search_node = urgap.init_node("Xtandem_Alanine:1.0.0")
    with pytest.raises(ValueError):
        search_node.run([mgf, fasta, params_1, params_2], urd)