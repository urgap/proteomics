import urgap


def test_wrapper_peptide_forest_3_simple():
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data/peptide_forest_input?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_CSV}"
                f"#a5fadd6fc4b64af6a2137a0a281397a0_1_of_1.pyiohat.csv"
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data/peptide_forest_input?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_JSON}"
                f"#a5fadd6fc4b64af6a2137a0a281397a0_1_of_1.pyiohat.json"
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data/peptide_forest_input?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_CSV}"
                f"#ba82cff4ccb4f14ad47a66867989791b_1_of_1.pyiohat.csv"
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data/peptide_forest_input?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_JSON}"
                f"#ba82cff4ccb4f14ad47a66867989791b_1_of_1.pyiohat.json"
            ),
        ]
    )
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "PeptideForest:3.1.1": {
                    "q_cut": 0.01,
                    "q_cut_train": 0.1,
                    "sensitivity": 0.9,
                    "n_train": 1,
                    "n_test": 1,
                    "random_seed": 42,
                    "initial_engine": "msgfplus_2021_03_22",
                },
            },
            "unode_parameters": {
                "force": True,
            },
        }
    )
    peptideforest = urgap.init_unode("PeptideForest:3.1.1")
    pf_out = peptideforest.run([ufiles], urun_dict)
    assert pf_out[0].path.exists() is True
