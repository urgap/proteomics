"""Urgap sage wrapper."""

import json

from pathlib import Path

import urgap


class Sage(urgap.unode.UNodeBase):
    """Sage wrapper for the Sage search engine.

    Sage is a proteomics database search engine - a tool that transforms raw mass spectra from
    proteomics experiments into peptide identifications via database searching & spectral matching.
    """

    META_INFO = {
        "name": "Sage",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "0.14.7", "exe_path": "Sage/0_14_7/sage"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.ms.converter.mzml.ANY: {"min": 1, "max": -1},
            urgap.uftypes.proteomics.dbsearch.SAGE_PARAMS: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.FASTA: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.dbsearch.SAGE_TSV: {"min": 1, "max": 1},

        },
        "engine": None,
        "engine_type": ("db_search", "proteomics"),
        "citation": """Lazear, M. R. (2023). Sage: An open-source tool for fast proteomics searching and quantification at scale. Journal of Proteome Research, 22(11), 3652-3659. https://doi.org/10.1021/acs.jproteome.3c00486""",
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize Sage class."""
        super().__init__(*args, **kwargs)

    def preflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Preflight routine for Sage wrapper.

        During preflight,
            - parameters are formatted
            - mods are mapped and formatted

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.

        Raises:
            ValueError: If both a parameter file and command-line parameters
                in the urun_dict are provided. Please provide only one.
        """
        param_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.dbsearch.SAGE_PARAMS,
        )

        mzml_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.converter.mzml.ANY,
        )[0]
        fasta_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.FASTA,
        )[0]

        config_dict = json.loads(param_files[0].read_text())
        config_dict["output_directory"] = str(mzml_file.parent)

        config_path = param_files[0].parent / "sage_config.json"
        config_path.write_text(json.dumps(config_dict))
        utrace.urun_dict.command_list = [
            str(self.exe_path),
            str(config_path),
            "--fasta",
            str(fasta_file),
            str(mzml_file),


        ]
        return utrace



    def postflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Postflight routine for Sage wrapper.

        During postflight the Sage native .tsv output file is converted into the
        pre-defined urgap output file, which is of csv format.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        full_path = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.converter.mzml.ANY,
        )[0]
        sage_tsv = full_path.parent / "results.sage.tsv"
        self.tmp_files.append(sage_tsv)
        Path(sage_tsv).rename(utrace.output_files[0].path)
        return utrace
