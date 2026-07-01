"""Tests for the GenerateTargetDecoyFasta:2.0.0 unode."""

from pathlib import Path

import urgap


def test_wrapper_generate_target_decoy_fasta_2_0_simple(tmp_dir: Path) -> None:
    """Test that reverse_protein mode reverses each tryptic peptide in place."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.FASTA}#fastas/Hfvol_prot_230328_bridgeCDS_cRAP_target-decoy_trypsin.fasta",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "GenerateTargetDecoyFasta:2.0.0": {
                    "--seed": "TEST_SEED",
                    "--enzyme_pattern": "trypsin",
                    "--mode": "reverse_protein",
                    "--decoy_tag": "decoy_",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    node = urgap.init_unode("GenerateTargetDecoyFasta:2.0.0")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True
    with output_files[0].path.open() as test_file:
        lines = list(test_file)
        assert len(lines) == 85552
    with output_files[1].path.open() as immutable_file:
        immutable_lines = list(immutable_file)

    # reverse_protein-specific: row 11 is the first line of the first decoy
    # peptide, which should be the target's tail peptide reversed
    assert lines[10].startswith("NDFRAQVF")

    assert len(immutable_lines) == 0


def test_wrapper_generate_target_decoy_fasta_2_0_params(tmp_dir: Path) -> None:
    """Test that shuffle_peptide mode does not produce a reversal of the target."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.FASTA}#fastas/Hfvol_prot_230328_bridgeCDS_cRAP_target-decoy_trypsin.fasta",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "GenerateTargetDecoyFasta:2.0.0": {
                    "--seed": "TEST_SEED",
                    "--enzyme_pattern": "trypsin",
                    "--mode": "shuffle_peptide",
                    "--decoy_tag": "decoy_",
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
            },
        },
    )
    node = urgap.init_unode("GenerateTargetDecoyFasta:2.0.0")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True
    with output_files[0].path.open() as test_file:
        lines = list(test_file)
        assert len(lines) == 85552
    with output_files[1].path.open() as immutable_file:
        immutable_lines = list(immutable_file)

    # shuffle_peptide test:
    assert len(immutable_lines) == 0
    # shuffle_peptide-specific: row 11 should NOT match the reverse_protein
    # output for the same peptide, confirming --mode actually changed behavior
    assert not lines[10].startswith("NDFRAQVF")
