"""Integration test for MSFragger."""


import pandas as pd


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


