import tempfile

import numpy as np
import pandas as pd
import pytest

import urgap


def test_spectrum_meta_data_wrapper_bsa_mzml():
    with tempfile.TemporaryDirectory() as tmp_dir:
        ufiles = urgap.UFileList(
            [
                urgap.UFile(
                    uri=f"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/examples/BSA"
                    f"?uftype={urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML}#BSA1.mzML"
                )
            ]
        )
        urd = urgap.URunDict(
            {
                "parameters": {
                    "ExtractSpectrumMetaData:1.0.0": {"-tf": "%Y-%m-%d %H:%M:%S"},
                },
                "unode_parameters": {
                    "force": True,
                    "storage_base_uri": f"file://{tmp_dir}",
                },
            }
        )
        test_node = urgap.init_unode("ExtractSpectrumMetaData:1.0.0")
        response = test_node.run(ufiles, urd)
        expected = {"run_shape": (1, 6), "size": 13.64, "median_mz": 532.68}

        # 1. Assert the run metadata
        assert response[0].path.exists()
        df = pd.read_csv(response[0].path)
        assert df.shape == expected["run_shape"]
        assert df["file"][0] == ufiles[0].object_name == df["lineage_root"][0]
        assert pytest.approx(df["file_size_mb"][0], 0.01) == expected["size"]
        # assert np.isnan(df.stop_time[0])

        # 2. Assert the spectrum metadata
        assert response[1].path.exists()
        df2 = pd.read_csv(response[1].path)
        assert pytest.approx(df2.precursor_mz.median(), 0.01) == expected["median_mz"]
