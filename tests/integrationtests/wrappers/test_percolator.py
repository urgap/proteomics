"""Tests for the Percolator:3.7.1 unode."""

from pathlib import Path

import pandas as pd

import urgap


def test_wrapper_percolator_3_7_1_concatenated(tmp_dir: Path) -> None:
    """Test using a concatenated search-input (referring to the target-decoy database)."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYIOHAT_CSV}"
                "#pyiohat_data/pyiohat.csv",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYIOHAT_JSON}"
                "#pyiohat_data/pyiohat_results.json",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Percolator:3.7.1": {
                    "delimiter": "<|>",
                    "search-input": "concatenated",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
                "remove_temporary_files": False,
            },
        },
    )
    node = urgap.init_unode("Percolator:3.7.1")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True
    df = pd.read_csv(output_files[0].path, index_col=0)

    assert len(df) > 1
    assert "q-value" in df.columns
    assert "posterior_error_prob" in df.columns

    # concatenated-specific: the PEP distribution reaches 1.0
    # (separate tops out at ~0.918)
    assert df["posterior_error_prob"].max() >= 0.999
    # concatenated-specific: exact count of PSMs at 1% FDR (separate: 15346)
    assert (df["q-value"] <= 0.01).sum() == 1398


def test_wrapper_percolator_3_7_1_separate(tmp_dir: Path) -> None:
    """Test using separate as search input resulting in different q-values and PEP."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYIOHAT_CSV}"
                "#pyiohat_data/pyiohat.csv",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYIOHAT_JSON}"
                "#pyiohat_data/pyiohat_results.json",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "Percolator:3.7.1": {
                    "delimiter": "<|>",
                    "search-input": "separate",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
                "remove_temporary_files": False,
            },
        },
    )
    node = urgap.init_unode("Percolator:3.7.1")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True
    df = pd.read_csv(output_files[0].path, index_col=0)

    assert len(df) > 1
    assert "q-value" in df.columns
    assert "posterior_error_prob" in df.columns

    # separate-specific: the PEP distribution never reaches 1.0
    # (concatenated hits 1.0)
    assert df["posterior_error_prob"].max() < 0.95
    # separate-specific: exact count of PSMs at 1% FDR (concatenated: 15328)
    assert (df["q-value"] <= 0.01).sum() == 1420
