"""Urgap pyiohat_1_9_0 unode."""

import json
import logging

from pathlib import Path

import pandas as pd

import urgap


class Pyiohat(urgap.unode.UNodeBase):
    """Urgap unode for the pyiohat_1_9_0 resource.

    This unode calls the main resource to pyProtista csv files coming from different
    proteomics search engines. The purpose is to bring the search engine specific
    names into a unified format for further processing and merging of data.
    """

    META_INFO = {
        "name": "Pyiohat",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "1.9.0", "exe_path": "pyiohat/1_9_0/pyiohat_resource.py"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.dbsearch.ANY: {
                "min": 0,
                "max": 1,
            },
            urgap.uftypes.proteomics.quantification.FLASHLFQ_PSM_TSV: {
                "min": 0,
                "max": 1,
            },
            urgap.uftypes.proteomics.MODS_XML: {"min": 0, "max": -1},
            urgap.uftypes.ms.SPECTRA_META_CSV: {"min": 1, "max": -1},
            urgap.uftypes.proteomics.FASTA: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.converter.PYIOHAT_JSON: {"min": 1, "max": 1},
        },
        "engine_type": ("converter", "proteomics"),
        "citation": "Urgap team (2021)",
        "release_date": "12.10.2022",
    }

    def __init__(self, *args: str, **kwargs: str) -> None:
        """Initialize Pyiohat class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for Pyiohat unode.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        search_input_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.dbsearch.ANY,
        )
        quant_input_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.quantification.FLASHLFQ_PSM_TSV,
        )
        if len(search_input_file) == 1 and len(quant_input_file) == 0:
            input_file = search_input_file[0]
        elif len(search_input_file) == 0 and len(quant_input_file) == 1:
            input_file = quant_input_file[0]
        else:
            logging.warning("Needs at least one quant input or one search input file")

        fasta_files = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.FASTA,
        )
        fasta_file = fasta_files[0]
        meta_data_files = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.ms.SPECTRA_META_CSV,
        )

        # will probably be overridden for multiple runs
        # make temp file or add hash?
        concatenated_meta = pd.concat([pd.read_csv(file) for file in meta_data_files])
        tmp_md_file = str(utrace.output_files[0].path) + "_md.csv"
        concatenated_meta.to_csv(tmp_md_file, index=False)
        self.tmp_files.append(tmp_md_file)

        output_file = utrace.output_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.converter.PYIOHAT_CSV,
        )[0]
        metadata_file = utrace.output_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.converter.PYIOHAT_JSON,
        )[0]

        params_dict = utrace.urun_dict.parameters[f"{self.META_INFO['unode_full_identifier']}"]


        use_param_file = utrace.urun_dict.unode_parameters.get("use_param_file", False)
        if use_param_file:
            tmp_params_file = str(utrace.output_files[0].path) + "_params.json"
            with Path(tmp_params_file).open("w", encoding="utf-8") as f:
                json.dump(params_dict, f)
            self.tmp_files.append(tmp_params_file)
            param_arg = tmp_params_file
        else:
            param_arg = json.dumps(params_dict)

        utrace.urun_dict.command_list = [
            "python",
            str(self.exe_path),
            "--input_file",
            str(input_file),
            "--spectrum_meta_data",
            str(tmp_md_file),
            "--output_file",
            str(output_file),
            "--fasta_file",
            str(fasta_file),
            "--metadata_file",
            str(metadata_file),
            "--parameters",
            param_arg,
        ]


        return utrace
