"""Urgap PeptideForest wrapper."""

import json
import os

import urgap


class PeptideForest(urgap.unode.UNodeBase):
    """Urgap wrapper for the PeptideForest module.

    Peptide Forest is a machine learning based tool for semisupervised integration of
    multiple peptide identification search engines. See publication
    provided under META_INFO["citation"] for further info.

    Requires:
        unimod_mapper
        peptide_forest
    """

    META_INFO = {
        "name": "PeptideForest",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {
                "version": "3.1.1",
                "exe_path": "peptide_forest/3_1_1/peptide_forest_3_1_1.py",
            }
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV: {"min": 1, "max": -1},
            urgap.uftypes.proteomics.converter.PYIOHAT_JSON: {"min": 1, "max": -1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.validator.PEPTIDEFOREST_CSV: {
                "min": 1,
                "max": 1,
            }
        },
        "engine": None,
        "engine_type": ("validation", "proteomics"),
        "citation": "Urgap team (2021)",
    }

    def __init__(self, *args: str, **kwargs: str):
        """Initialize PeptideForest class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for PeptideForest UNode.

        During preflight,
            - peptideforest config_json is written
            - command list is composed

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        config_json = self._write_config_json(utrace)
        output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.validator.PEPTIDEFOREST_CSV
        )[0]
        utrace.urun_dict.command_list = [
            "python",
            str(self.exe_path),
            "-c",
            str(config_json),
            "-o",
            str(output_file),
        ]
        return utrace

    def _write_config_json(
        self,
        utrace: urgap.UTrace,
    ) -> os.PathLike:
        """Write config_json required by peptideforest.

        The function formats user input parameters into peptideforest style and writes
        them out into a config_json file required for peptideforest execution.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            Path to config_json.
        """
        data = {"input_files": {}}
        initial_engine = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]["initial_engine"]
        data.update({"initial_engine": initial_engine, "column_mapping": {}})
        input_csvs = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV
        )
        input_metadata_jsons = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYIOHAT_JSON
        )
        for file in input_csvs:
            file_dict = {str(file): {}}
            file_name = file.name.rstrip(urgap.uftypes.proteomics.converter.PYIOHAT_CSV)
            for json_file in input_metadata_jsons:
                if json_file.name.startswith(file_name):
                    with open(json_file, "r") as fh:
                        json_content = json.load(fh)
                engine = json_content["Parser"].split("/")[-1].rstrip("_parser.py")
                if engine == "msfragger_4":
                    engine = "msfragger_4_2"
                file_dict[str(file)].update(
                    {
                        "engine": engine,
                        "score_col": json_content["validation_score_field"],
                        "bigger_score_better": json_content["bigger_scores_better"],
                    }
                )
            data["input_files"].update(file_dict)
        config_path = utrace.output_files[0].path.parent / "config.json"
        with open(config_path, "w") as fh:
            json.dump(data, fh)
        return config_path
