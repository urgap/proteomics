#!/usr/bin/env python3
import pytest

import urgap


@pytest.mark.parametrize("parameters", [{}, {"--msLevel": 2}])
def test_wrapper_thermorawfileparser_2_0_0_simple(tmp_dir, parameters):
    ufiles = urgap.UFileList(
        [
            urgap.UFile(
                uri=f"file://{urgap._test_folder}/data?uftype="
                f"{urgap.uftypes.proteomics.THERMO_RAW}#ms_files/small.raw"
            )
        ]
    )
    urun_dict = urgap.URunDict(
        {
            "parameters": {
                "ThermoRawFileParser:2.0.0": parameters,
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
    if parameters == {}:
        assert len(lines) == 2620
    elif parameters == {"--msLevel": 2}:
        assert len(lines) == 2066
