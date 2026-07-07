#!/usr/bin/env python3

# test_node = ursgal.init_node("percolator_3_5_0").download_engine_on_demand()


# @pytest.mark.slow
# @pytest.mark.skipif(
#     not (
#         ursgal.home
#         / "resources"
#         / sys.platform
#         / platform.architecture()[0]
#         / "percolator_3_5_0"
#         / "percolator"
#     ).exists(),
#     reason=f"Binary not " f"available for {sys.platform} {platform.architecture()[0]}",
# )
# @pytest.mark.parametrize(
#     "provide_clean_node_dirs",
#     [
#         (
#             (
# ursgal.UFile(
#     uri=f"https://raw.githubusercontent.com/computational-ms/demo_files/main/unified_csvs?uftype={ursgal.uftypes.proteomics.converter.PYIOHAT_CSV}#JB_FASP_pH8_2-3_28122012_msgfplus_2021_03_22_unified.csv"
# ),
#                 ursgal.UFile(
#                     uri=f"file://{ursgal._test_folder}/data?uftype={ursgal.uftypes.proteomics.FASTA}#fastas/BSA1.fasta"
#                 ),
#             ),
#             ursgal.URunDict(
#                 {
#                     "parameters": {
#                         "infer_proteins": True,
#                     },
#                     "unode_parameters": {
#                         "remove_temporary_files": True,
#                     },
#                 }
#             ),
#             ["percolator_3_5_0"],
#         ),
#         (
#             (
# ursgal.UFile(
#     uri=f"https://raw.githubusercontent.com/computational-ms/demo_files/main/unified_csvs?uftype={ursgal.uftypes.proteomics.converter.PYIOHAT_CSV}#JB_FASP_pH8_2-3_28122012_msgfplus_2021_03_22_unified.csv"
# ),
#             ),
#             ursgal.URunDict(
#                 {
#                     "parameters": {
#                         "infer_proteins": False,
#                     },
#                     "unode_parameters": {
#                         "remove_temporary_files": True,
#                     },
#                 }
#             ),
#             ["percolator_3_5_0"],
#         ),
#     ],
#     indirect=["provide_clean_node_dirs"],
# )
# def test_wrapper_percolator_3_5_0(provide_clean_node_dirs):
#     nodes, ufiles, urun_dict = provide_clean_node_dirs
#     node_name = "percolator_3_5_0"
#     percolator_node = nodes[node_name]

#     output_files = percolator_node.run(
#         ufiles=ufiles,
#         urun_dict=urun_dict,
#     )
#     for output_file in output_files:
#         assert output_file.path.exists() is True
