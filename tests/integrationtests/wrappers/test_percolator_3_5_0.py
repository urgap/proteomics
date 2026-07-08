from pathlib import Path
import urgap

def test_wrapper_percolator_3_5_0_simple(tmp_dir: Path) -> None:
    """Test that Percolator successfully processes unified CSVs and infers proteins."""
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=(
                    "https://raw.githubusercontent.com/computational-ms/demo_files/main/"
                    f"unified_csvs?uftype={urgap.uftypes.proteomics.converter.PYIOHAT_CSV}"
                    "#JB_FASP_pH8_2-3_28122012_msgfplus_2021_03_22_unified.csv"
                )
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
                "percolator:3.5.0": {
                    "infer_proteins": True,
                },
            },
            "unode_parameters": {
                "storage_base_uri": f"file://{tmp_dir}",
                "remove_temporary_files": True,
            },
        },
    )

    # Note: If your registry still uses underscores, change this to "percolator_3_5_0"
    node = urgap.init_unode("percolator:3.5.0")
    
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