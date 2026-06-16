#!/usr/bin/env python3
import pytest

import urgap

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
