"""Integration test for Sage."""
import pandas as pd
from pathlib import Path
import urgap

import pytest


def test_searchdb_sage(tmp_path: Path) -> None:
    params = urgap.UFile(uri=f"file:///{urgap._test_folder}/data/sage_params?uftype={urgap.uftypes.proteomics.params.SAGE_JSON}#sage.params")
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
    search_result = search_node.run([params], urd)
    assert search_result[0].path.exists() is True

    # check that output is not empty
    df = pd.read_csv(search_result[0].path, delimiter="\t")
    assert df.shape[0] > 13

    search_node.remove_output_folder(output_file=search_result[0])


