"""urgap CombinePEP 1.0.0 wrapper."""

import pandas as pd

import urgap


class CombinePEP(urgap.unode.UNodeBase):
    """urgap wrapper for the combine_pep_1_0_0 executable.

    Combines posterior error probabilities (PEPs) across multiple
    Percolator result CSVs (one per search engine) into unified Bayes
    PEP and windowed combined PEP estimates.
    """

    META_INFO = {
        "name": "CombinePEP",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "1.0.0", "exe_path": "combine_pep/1_0_0/combine_pep_1_0_0.py"},
        ],
        "parameters_not_triggering_rerun": [],
        "engine_type": ("validation", "proteomics"),
        "platform_independent": True,

        "input_uftypes": {
            urgap.uftypes.proteomics.validator.PERCOLATOR_CSV: {"min": 1, "max": -1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.validator.COMBINEPEP_CSV: {"min": 1, "max": 1},
        },
        "utranslation_style": "combine_pep_style_1",
        "citation": "Kremer, L. P. M., Leufken, J., Oyunchimeg, P., Schulze, S., & Fufezan, C. (2016). Ursgal, Universal Python Module Combining Common Bottom-Up Proteomics Tools for Large-Scale Analysis."
        " Journal of Proteome Research, 15(3), 788-794. https://doi.org/10.1021/acs.jproteome.5b00860 ",
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize CombinePEP node."""
        super().__init__(*args, **kwargs)

    def preflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Build the combine_pep_1_0_0 command list from run parameters.

        Reads grouping columns and window size from the run parameters,
        infers the search engine for each input Percolator CSV, and
        composes the CLI command used to invoke combine_pep_1_0_0.py.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        params_dict = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]

        config = {
            "input_csvs": [],
            "input_engines": [],
            "columns_for_grouping": params_dict["psm_defining_colnames"],
            "pep_colname": "posterior_error_prob",
            "input_sep": ",",
            "output_sep": ",",
            "join_sep": ";",
            "window_size": params_dict["window_size"],
        }

        for infile in utrace.input_files:
            config["input_csvs"].append(str(infile.path))
            df = pd.read_csv(infile.path, nrows=1)
            config["input_engines"].append(df.iloc[0]["search_engine"])

        output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.validator.COMBINEPEP_CSV,
        )[0]

        utrace.urun_dict.command_list = [
            "python", str(self.exe_path),
            "-i", *config["input_csvs"],
            "-e", *config["input_engines"],
            "-c", *config["columns_for_grouping"],
            "-o", str(output_file),
            "-is", config["input_sep"],
            "-os", config["output_sep"],
            "-js", config["join_sep"],
            "-w", str(config["window_size"]),
        ]
        return utrace
