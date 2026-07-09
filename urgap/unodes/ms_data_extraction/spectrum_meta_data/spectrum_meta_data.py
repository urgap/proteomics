"""Urgap ExtractSpectrumMetaData wrapper."""

import logging

import urgap


class ExtractSpectrumMetaData(urgap.unode.UNodeBase):
    """Urgap wrapper for the ExtractSpectrumMetaData resource.

    This wrapper calls the extract_meta_data resource from simepy package to extract
    peak information from an input raw or mzml file. Note: for raw_file extraction the
    cz_xcalibur package is required.
    """

    META_INFO = {
        "name": "ExtractSpectrumMetaData",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "1.0.0", "exe_path": "ms_data_extraction/1_0_0/extract_meta_data.py"},
        ],
        "parameters_not_triggering_rerun": [],
        # ^-- TODO: this does not seem to work. All tests are skipped, even if working
        #       with mzml files only
        "input_uftypes": {
            urgap.uftypes.proteomics.THERMO_RAW: {
                "min": 0,
                "max": 1,
            },
            urgap.uftypes.any.MZML: {
                "min": 0,
                "max": 1,
            },
            urgap.uftypes.ms.converter.mzml.ANY: {
                "min": 0,
                "max": 1,
            },
        },
        "output_uftypes": {
            urgap.uftypes.ms.RUN_META_CSV: {"min": 1, "max": 1},
            urgap.uftypes.ms.SPECTRA_META_CSV: {"min": 1, "max": 1},
            urgap.uftypes.ms.SPECTRA_NOISE_CSV: {"min": 1, "max": 1},
            urgap.uftypes.ms.INSTRUMENT_UNIT_CSV: {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("data_extractor",),
        "citation": "Urgap team (2021)",
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize ExtractSpectrumMetaData class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for ExtractSpectrumMetaData wrapper.

        Checks that only one file at a time is processed and if so, provides the file
        path to the main function of the resource to extract metadata from.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        index_dict = utrace.input_files.get_index_groups_by_uftypes()
        uftype_to_process = None
        for uftype, index_list in index_dict.items():
            if len(index_list) > 0:
                if uftype_to_process is None:
                    uftype_to_process = uftype
                else:
                    logging.warning(
                        "Meta node received multiple mass spec files, "
                        "however node can only process one at the time!"
                        "Therefore will process %s but not %s!",
                        uftype_to_process,
                        uftype,
                    )

        idx_ms = index_dict[uftype_to_process][0]
        input_file = utrace.input_files[idx_ms]

        run_meta_output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.RUN_META_CSV,
        )[0]
        spec_meta_output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.SPECTRA_META_CSV,
        )[0]
        spec_noise_output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.SPECTRA_NOISE_CSV,
        )[0]
        instrument_unit_output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.INSTRUMENT_UNIT_CSV,
        )[0]
        object_name = input_file.object_name
        utrace.urun_dict.command_list = [
            "python",
            str(self.exe_path),
            "-i",
            str(input_file.path),
            "-ro",
            str(run_meta_output_file),
            "-so",
            str(spec_meta_output_file),
            "-sno",
            str(spec_noise_output_file),
            "-iuo",
            str(instrument_unit_output_file),
            "-on",
            str(object_name),
            "-lr",
            str(object_name),
        ]
        for parameter_key, parameter_value in utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ].items():
            utrace.urun_dict.command_list.extend([parameter_key, parameter_value])

        return utrace
