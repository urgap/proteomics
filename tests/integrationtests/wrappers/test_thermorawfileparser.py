#!/usr/bin/env python3
import pytest

import urgap

#
# @pytest.mark.slow
# @pytest.mark.skipif(
#     urgap.init_node("thermo_raw_file_parser_2_0_0").has_all_required_installations()
#     is False,
#     reason="Node Requires {0}".format(
#         urgap.init_node("thermo_raw_file_parser_2_0_0").META_INFO["requires"]
#     ),
# )
#
# pytest.mark.parametrize(
#     "provide_clean_node_dirs",
#     [
#         (
#             urgap.UFile(
#                 uri=f"file://{urgap._test_folder}/data?uftype="
#                 f"{urgap.uftypes.proteomics.THERMO_RAW}#mzML/BSA1.raw"
#             ),
#             urgap.URunDict(
#                 {
#                     "parameters": {
#                         "ThermoRawFileParser:2.0.0": {},
#                     },
#                     "unode_parameters": {
#                         "force": True,
#                     },
#                 }
#             ),
#             ["ThermoRawFileParser:2.0.0"],
#         )
#     ],
#     indirect=["provide_clean_node_dirs"],
# )


def test_wrapper_thermorawfileparser_2_0_0_simple(tmp_dir):
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.THERMO_RAW}#raw/small.raw"
            )
        ]
    )
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "ThermoRawFileParser:2.0.0": {},
            },
            "unode_parameters": {
                "force": True,
            },
        }
    )
    node = urgap.init_unode("ThermoRawFileParser:2.0.0")
    output_files = node.run(
        ufiles,
        urun_dict,
    )
    assert output_files[0].path.exists() is True
    with open(output_files[0].path, "r") as test_file:
        lines = [l for l in test_file]
        assert len(lines) == 2620
