"""Integration test for Pyiohat."""


import pandas as pd

import pytest 

import urgap




def test_pyiohat(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )


    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
    assert unified_files[0].path.exists() is True




def test_pyiohat_metadata_file(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )


    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
    assert unified_files[1].path.exists() is True




def test_pyiohat_wo_xml_file(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )


    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
    assert unified_files[0].path.exists() is True




def test_pyiohat_wo_xml_file_metadata_file(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )


    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
    assert unified_files[1].path.exists() is True
def test_pyiohat_with_param_file(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {},
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_JSON}#pyiohat_json"
            f"/pyiohat_params.json"
        ),
    )
 
    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
    assert unified_files[0].path.exists() is True
 
 
def test_pyiohat_param_file_and_urun_dict_params_raises(tmp_dir):
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_JSON}#pyiohat_json"
            f"/pyiohat_params.json"
        ),
    )
 
    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    with pytest.raises(ValueError):
        pyiohat_node.run(ufiles, urd)
import pandas as pd
import pytest
 
import urgap
 
 
def test_pyiohat_output_content(tmp_dir):
    """Check that the unified pyiohat output csv actually contains the
    expected PSM-level columns and data, not just that the file exists.
    """
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[KR])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )
 
    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
 
    df = pd.read_csv(unified_files[0].path)
 
    # Make sure we actually got rows back, not an empty file.
    assert len(df) == 91
    assert df["modifications"].str.contains("Carbamidomethyl:").sum() == 51
    assert pytest.approx(df["ucalc_mz"].mean(), abs=1e-3) == 468.709
    assert pytest.approx(df["exp_mz"].mean(), abs=1e-3) == 486.985
    assert (df["raw_data_location"] == "mzMLs/BSA1.mzML").all()
    assert df["enzn"].all()
    assert df["enzc"].all()
 


 
 
def test_pyiohat_altered_enzyme_param_changes_output(tmp_dir):
    """Changing the enzyme cleavage-site regex should change which PSMs/
    peptides survive into the unified output, demonstrating that the
    parameter is actually being passed through and respected.
    """
    urd = urgap.URunDict(
        {
            "parameters": {
                "Pyiohat:1.9.0": {
                    "enzyme": "(?<=[C])(?![P])",
                    "terminal_cleavage_site_integrity": "any",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    ufiles = (
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID}#mzids/BSA1_msgfplus_2021_03_22.mzid"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
        ),
        urgap.UFile(
            uri=f"file://"
            f"{urgap._test_folder}/data?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods/usermods.xml"
        ),
        urgap.UFile(
            uri=f"file://{urgap._test_folder}/data?uftype={urgap.uftypes.ms.SPECTRA_META_CSV}#meta_data"
            f"/5832277fafa758daf43584d38502eb47_1.spectra_meta.csv"
        ),
    )
 
    pyiohat_node = urgap.init_node("Pyiohat:1.9.0")
    unified_files = pyiohat_node.run(ufiles, urd)
 
    df = pd.read_csv(unified_files[0].path)
 
    # Make sure we actually got rows back, not an empty file.
    assert len(df) == 91
    assert not df["enzn"].any()
    assert df["enzc"].sum() == 2
