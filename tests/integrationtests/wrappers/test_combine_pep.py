"""Tests for the CombinePEP:1.0.0 unode."""

import csv

from pathlib import Path

import pytest

import urgap


def test_wrapper_combine_pep_1_0_0_two_engines(tmp_dir: Path) -> None:
    """Test that CombinePEP correctly combines PEPs across Casanovo and X!Tandem."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.validator.PERCOLATOR_CSV}#percolator_files/percolator_casanovo.csv",
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.validator.PERCOLATOR_CSV}#percolator_files/percolator_xtandem.csv",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "CombinePEP:1.0.0": {
                    "psm_defining_colnames": [
                        "spectrum_id", "sequence", "modifications", "charge",
                    ],
                    "window_size": 100,
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    node = urgap.init_unode("CombinePEP:1.0.0")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True

    with output_files[0].path.open(newline="") as test_file:
        rows = list(csv.DictReader(test_file))

    assert "Bayes PEP" in rows[0]
    assert "combined PEP" in rows[0]
    assert "combined PEP engines" in rows[0]

    # 2971 casanovo-only + 2977 xtandem-only + 43 shared PSMs written once per
    # engine (86 rows) = 6034 data rows
    assert len(rows) == 6034

    # PSM (spectrum_id=7498, sequence=DGTITETDGSTR, modifications="", charge=2)
    # is one of the 43 PSMs shared by both engines. Its two output rows should
    # both carry the same Bayes PEP, computed from casanovo PEP 0.000193267 and
    # xtandem PEP 0.180512.
    shared_psm_rows = [
        row for row in rows
        if row["spectrum_id"] == "7498" and row["sequence"] == "DGTITETDGSTR"
    ]
    assert len(shared_psm_rows) == 2

    for row in shared_psm_rows:
        assert row["combined PEP engines"] == "casanovo_5_0;xtandem_alanine"
        assert float(row["Bayes PEP"]) == pytest.approx(4.2578e-05, rel=1e-3)


def test_wrapper_combine_pep_1_0_0_single_engine(tmp_dir: Path) -> None:
    """Test that CombinePEP runs correctly with only one engine's PSMs."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.validator.PERCOLATOR_CSV}#percolator_files/percolator_xtandem.csv",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "CombinePEP:1.0.0": {
                    "psm_defining_colnames": [
                        "spectrum_id", "sequence", "modifications", "charge",
                    ],
                    "window_size": 100,
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    node = urgap.init_unode("CombinePEP:1.0.0")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True

    with output_files[0].path.open(newline="") as test_file:
        rows = list(csv.DictReader(test_file))

    # single-engine run: 3020 unique PSM keys (dict-collapsed from 3023 raw rows)
    assert len(rows) == 3020

    # single-engine test: every row's "combined PEP engines" value should be
    # exactly the one engine name
    assert rows[0]["combined PEP engines"] == "xtandem_alanine"
