"""Integration test for Sage."""
import pandas as pd
from pathlib import Path
import urgap

import pytest


def test_searchdb_sage(tmp_path: Path) -> None:
    params = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/params?uftype={urgap.uftypes.proteomics.dbsearch.SAGE_PARAMS}#sage_config.json")
    params_alt = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/params?uftype={urgap.uftypes.proteomics.dbsearch.SAGE_PARAMS}#sage_config_alternate.json")
    mzml = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/ms_files?uftype={urgap.uftypes.ms.converter.mzml.ANY}#BSA1.mzML")
    fasta = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/fastas?uftype={urgap.uftypes.proteomics.FASTA}#BSA1.fasta")
    urd = urgap.URunDict(
        {
            "parameters": {"Sage:0.14.7": {
                },
            },
            "unode_parameters": { 
                "storage_base_uri": f"file://{tmp_path}",
            },
        },
    )
    search_node = urgap.init_node("Sage:0.14.7")
    
    for parameters in (params, params_alt):
        search_result = search_node.run([parameters, mzml, fasta], urd)
        assert search_result[0].path.exists() is True
    
        # check that output is not empty
        df = pd.read_csv(search_result[0].path, delimiter="\t")
        if parameters == params:
            assert df.shape[0] == 103
            
        elif parameters == params_alt:
            assert df.shape[0] == 3
            
        search_node.remove_output_folder(output_file=search_result[0])


