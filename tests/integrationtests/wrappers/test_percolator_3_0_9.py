from pathlib import Path
import urgap

def test_wrapper_percolator_3_7_1_simple(tmp_dir: Path) -> None:
    """Test that Percolator successfully processes unified CSVs and infers proteins."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.converter.PYIOHAT_CSV}#unified_csvs/human_ecoli_sample_pyiohat.csv"
            ),
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.FASTA}#fastas/BSA1.fasta",
            ),
        ],
    )

    urun_dict = urgap.URunDict(
        {
            "parameters": {
                # Nesting parameters under the node identifier to match the new style
            "Percolator:3.7.1": {
                "infer_proteins": True,
                "delimiter": "<|>",
                "bigger_scores_better": {"omssa_2_1_9": False},
                "validation_score_field": {"omssa_2_1_9": "omssa:evalue"},
            },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
                "remove_temporary_files": False,
            },
        },
    )

    # Note: If your registry still uses underscores, change this to "percolator_3_7_1"
    node = urgap.init_node("Percolator:3.7.1")
    
    output_files = node.run(
        ufiles,
        urun_dict,
    )

    # Verify the output target file was created
    assert output_files[0].path.exists() is True
    
    # Read the file to ensure content was generated
    with output_files[0].path.open() as test_file:
        lines = list(test_file)
        assert len(lines) > 0  # Swap with an exact line check once you know the target coun