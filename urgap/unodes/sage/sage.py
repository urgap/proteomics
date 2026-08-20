"""Urgap sage wrapper."""


import urgap

import json

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
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML: {"min": 0, "max": -1},
            urgap.uftypes.proteomics.params.SAGE_JSON: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.FASTA: {"min": 0, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.dbsearch.SAGE_TSV: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.params.SAGE_JSON: {"min": 0, "max": 1},
            
        },
        "engine": None,
        "engine_type": ("identification",),
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
            - param file is written

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.

        Raises:
            ValueError: If both a parameter file and command-line parameters
                in the urun_dict are provided. Please provide only one.
        """
        params_dict = utrace.urun_dict.parameters[
            f"{self.META_INFO['name']}:{self.META_INFO['versions'][0]['version']}"
        ]
        
        param_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.params.SAGE_JSON,
        )
        param_file_provided = len(param_files) == 1
        cmdline_json_provided = "-json" in params_dict
        
          
        if param_file_provided and cmdline_json_provided: 
            raise ValueError(
                "Both parameter and command-line parameter provided."
                "Please only use one",
            )
        elif param_file_provided:
            param_json_path = param_files[0]
        elif cmdline_json_provided:
            param_json_path = params_dict["-json"]
        else:
            raise ValueError(
                "No parameters provided."
                "Please provide parameters for Sage",
            )            
          
        mzml_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML,
        )[0]

        config_dict = json.loads(param_json_path.read_text())
        config_dict["output_directory"] = str(mzml_file.parent)

        config_path = param_json_path.parent / "sage_config.json"
        config_path.write_text(json.dumps(config_dict))
        utrace.urun_dict.command_list = [
            str(self.exe_path),
            str(config_path),

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
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML,
        )[0]
        Sage_tsv = full_path.parent / "results.sage.tsv"
        self.tmp_files.append(Sage_tsv)

        with (
            Sage_tsv.open() as fin,
            utrace.output_files[0].path.open("w") as fout,
        ): 
            for line in fin:
                fout.write(line)
        return utrace