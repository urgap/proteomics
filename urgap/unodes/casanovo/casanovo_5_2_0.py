"""Urgap casanonvo_5_2_0 wrapper."""

import logging

import urgap


class Casanovo(urgap.unode.UNodeBase):
    """Urgap wrapper for the casanovo_5_2_0 search engine.

    Casanovo is a de novo peptide sequencing program that sequences peptides from an  mzML, mzXML, or mgf file into an mztab file.
    This requires a gpu to run, or else the program will not be able to complete and continue going forever. See publication provided under META_INFO[
    "citation"] for further info.
    """

    META_INFO = {
        "name": "Casanovo",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {
                "version": "5.2.0",
                "exe_path": "$casanovo",
            },
        ],
        "parameters_not_triggering_rerun": [],
        "engine": None,
        "input_uftypes": {
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.FASTA: {"min": 0, "max": 1},
        },

        "output_uftypes": {
            urgap.uftypes.proteomics.denovosearch.CASANOVO_MZTAB: {"min": 1, "max": 1},
        },
        "engine_type": ("identification",),
        "citation": """Yilmaz, M., Fondrie, W. E., Bittremieux, W., Melendez, C.F., Nelson, R., Ananth, V., Oh, S. & Noble, W. S. Sequence-to-sequence translation from mass spectra to peptides with a transformer model. in Nature Communications 15, 6427 (2024). doi:10.1038/s41467-024-49731-x""",
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize casanovo_5_2_0 class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for casanovo_5_2_0 wrapper.

        During preflight,
            - parameters are formatted
            - param file is written
            -make sure only 1 file type is chosesn and file type is correct one

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        input_params = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]
        mzml_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML,
        )
        if len(mzml_file) == 0:
            logging.error(
                "Provide an mzML file as input for casanovo. No mzML file has been provided",
            )
        elif len(mzml_file) == 1:
            spec_file = mzml_file[0]

        param_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.denovosearch.CASANOVO_YAML,
        )[0]

        output_mztab = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.denovosearch.CASANOVO_MZTAB,
        )[0]

        # the mode in which casanovo is in run
        search_mode = input_params["search_mode"]
        if search_mode not in ["sequence", "db-search"]:
            logging.error(
                "Unknown search mode %s. Search mode has to be either 'sequence' or 'db-search'",
                search_mode,
            )
        utrace.urun_dict.command_list = [
            "casanovo",
             f"{search_mode}",
            str(spec_file),
            "--config",
            param_file,
            "--output_dir",
            output_mztab.parent,
            "--output_root",
            output_mztab.name.replace(".mztab", ""),
        ]

        if search_mode == "db-search":
            fasta_file = utrace.input_files.get_path_objects_by_uftype(
                urgap.uftypes.proteomics.FASTA,
            )[0]
            utrace.urun_dict.command_list.insert(3, fasta_file)
            # check if this is behaving the way you'd like it to
            if len(fasta_file) == 0:
                logging.error("Please input a Fasta file for database searching.")
        elif search_mode == "sequence":
            # check ig this is the correct way of assessing the fasta files or the one above
            fasta_file_list = utrace.input_files.get_path_objects_by_uftype(
                urgap.uftypes.proteomics.FASTA,
            )
            if len(fasta_file_list) != 0:
                logging.error(
                    "A fasta file has been provided despite the search mode being set to sequence, i.e. de novo.",
                )

        return utrace

    def postflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Postflight routine for casanovo_5_1_2 wrapper. Currently a no-op."""
        return utrace
