"""Urgap sage wrapper."""


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
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML: {"min": 0, "max": -1},
            urgap.uftypes.proteomics.MZML: {"min": 0, "max": -1},
            urgap.uftypes.proteomics.FASTA: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.params.SAGE_JSON: {"min": 0, "max": 1},
        },
        "output_uftypes": {
            proteomics.dbsearch.SAGE_JSON: {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("identification",),
        "citation": """Lazear, M. R. (2023). Sage: An open-source tool for fast proteomics searching and quantification at scale. Journal of Proteome Research, 22(11), 3652–3659. https://doi.org/10.1021/acs.jproteome.3c00486""",
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
        json = params_dict["-json"]

        mzml_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.MZML,
        )[0]

        fasta_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.FASTA,
        )[0]

        param_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.dbsearch.SAGE_PARAMS,
        )

        cli_params_dict = {k: v for k, v in params_dict.items() if k != "-json"}
        param_file_provided = len(param_files) == 1
        cli_params_provided = bool(cli_params_dict)

        if param_file_provided and cli_params_provided:
            msg = (
                "Both a parameter file and command-line parameters in the "
                "urun_dict were provided for Sage. Please provide only one."
            )
            raise ValueError(msg)

        cli_param_args = []
        for key, value in cli_params_dict.items():
            value_str = "true" if value is True else "false" if value is False else str(value)
            cli_param_args.append(f"--{key}")
            cli_param_args.append(value_str)

        utrace.urun_dict.command_list = [
            "rust",
            f"-Json{json}",
            str(self.exe_path),
        ]

        if param_file_provided:
            original_param_file = param_files[0]
            patched_param_file = original_param_file.parent / (
                "patched_" + original_param_file.name
            )
            self.tmp_files.append(patched_param_file)

            found_database_name = False
            found_output_format = False
            with (
                original_param_file.open() as fin,
                patched_param_file.open("w") as fout,
            ):
                for line in fin:
                    stripped = line.strip()
                    if stripped.startswith("database_name"):
                        fout.write(f"database_name = {fasta_file}\n")
                        found_database_name = True
                    elif stripped.startswith("output_format"):
                        fout.write("output_format = tsv\n")
                        found_output_format = True
                    else:
                        fout.write(line)
                if not found_database_name:
                    fout.write(f"database_name = {fasta_file}\n")
                if not found_output_format:
                    fout.write("output_format = tsv\n")

            utrace.urun_dict.command_list.append(str(patched_param_file))
        else:
            cli_param_args.append("--database_name")
            cli_param_args.append(str(fasta_file))
            cli_param_args.append("--output_format")
            cli_param_args.append("tsv")

        utrace.urun_dict.command_list.extend(cli_param_args)
        utrace.urun_dict.command_list.append(str(mzml_file))

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
        Sage_tsv = full_path.parent / (full_path.stem + ".tsv")
        self.tmp_files.append(Sage_tsv)

        with (
            Sage_tsv.open() as fin,
            utrace.output_files[0].path.open("w") as fout,
        ):
            for line in fin:
                fout.write(line)
        return utrace