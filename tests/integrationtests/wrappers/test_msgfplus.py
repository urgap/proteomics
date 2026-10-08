import xml.etree.ElementTree as ETree

import urgap


def test_searchdb_msgfplus_carbamidomethyl_in_usermods():
    mgf = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/mgfs?uftype={urgap.uftypes.proteomics.converter.PYMZML_MGF}#BSA1.mgf"
    )
    fasta = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta"
    )
    unimod = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/usermods?uftype={urgap.uftypes.proteomics.MODS_XML}#unimod_no_carbamidomethyl.xml"
    )
    usermod = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/usermods?uftype={urgap.uftypes.proteomics.MODS_XML}#usermods_carbamidomethyl.xml"
    )
    params1 = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/params?uftype={urgap.uftypes.proteomics.params.MSGFPLUS_TXT}#MSGFPlus_params_no_Carbamidomethyl.txt"
    )
    params2 = urgap.UFile(
        uri=f"file:///{urgap._test_folder}/data/params?uftype={urgap.uftypes.proteomics.params.MSGFPLUS_TXT}#MSGFPlus_params.txt"
    )

    urd = urgap.URunDict(
        {
            "parameters": {
                "MsgfPlus:20240326": {},
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{urgap._test_folder}/data",
            },
        }
    )

    search_node = urgap.init_node("MsgfPlus:20240326")
    search_result = search_node.run(
        [mgf, fasta, params1],
        urd,
    )
    assert search_result[0].path.exists() is True

    tree = ETree.parse(search_result[0].path)
    root = tree.getroot()

    for ele in root.findall(".//{*}Modification"):
        assert ele[0].get("name") != "Carbamidomethyl"
    search_node.remove_output_folder(output_file=search_result[0])

    search_result = search_node.run(
        [mgf, fasta, params2],
        urd,
    )
    assert search_result[0].path.exists() is True
    tree = ETree.parse(search_result[0].path)
    root = tree.getroot()

    C_DETECTED = False
    for ele in root.findall(".//{*}Modification"):
        if ele[0].get("name") == "Carbamidomethyl":
            C_DETECTED = True
    assert C_DETECTED is True
    search_node.remove_output_folder(output_file=search_result[0])
