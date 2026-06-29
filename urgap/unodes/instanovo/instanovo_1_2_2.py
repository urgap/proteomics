"""Urgap Instanovo_1_2_2 wrapper."""

import importlib.util
import logging
import os
import sys

from pathlib import Path

import urgap


class Instanovo(urgap.unode.UNodeBase):
    """Urgap wrapper for the Instanovo de-novosearch engine.
 
    Instanovo is a de novo peptide sequencing program that sequences peptides from an mgf file into an csv file.
    This requires a gpu to run, or else the program will not be able to complete and continue going forever. See publication provided under META_INFO[
    "citation"] for further info.
    """

    META_INFO = {
        "name": "Instanovo",
        "wrapper_version": {"major": 1, "minor": 2, "patch": 2},
        "versions": [
            {"version": "1.2.2", "exe_path": "$instanovo"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.converter.PYMZML_MGF: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.denovosearch.INSTANOVO_CSV : {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("identification",),
        "citation": """Eloff, K., Kalogeropoulos, K., Mabona, A. et al. InstaNovo enables diffusion-powered de novo peptide sequencing in large-scale proteomics experiments.
        at Mach Intell 7, 565-579 (2025). https://doi.org/10.1038/s42256-025-01019-5""",
    }

    def __init__(self, *args: str, **kwargs: str) -> None:
        """Initialize Instanovo_1_2_2 class."""
        super().__init__(*args, **kwargs)

    def find_instanovo_configs(self) -> Path:
        """Locate the instanovo package's built-in configs directory.
 
        Returns:
            Path: Resolved path to the instanovo/configs directory.
 
        Raises:
            FileNotFoundError: If the instanovo/configs directory could not
                be located via the package spec or sys.path.
        """
        # try to find the package spec first (works if package is installed/importable)
        spec = importlib.util.find_spec("instanovo")
        if spec and spec.submodule_search_locations:
            pkg_path = Path(spec.submodule_search_locations[0])
            cfg = pkg_path / "configs"
            if cfg.exists():
                return cfg.resolve()

        # fallback: scan sys.path entries for a possible instanovo/configs folder
        for entry in map(Path, sys.path):
            candidate = entry / "instanovo" / "configs"
            if candidate.exists():
                return candidate.resolve()

        msg = "Could not locate instanovo/configs on this system"
        raise FileNotFoundError(msg)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for Instanovo_1 wrapper.
 
        During preflight,
            - parameters are formatted
            - mods are mapped and formatted
            - param file is written
 
        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.
 
        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        input_params = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]
        mgf_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYMZML_MGF,
        )[0]
        param_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.denovosearch.INSTANOVO_YAML,
        )[0]

        output_csv = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.denovosearch.INSTANOVO_CSV,
        )[0]

        # Get the directory containing YOUR custom config file
        config_dir = Path(param_file).parent.resolve()

        # InstaNovo executes from its built-in configs directory, so we need
        # to calculate the relative path from there to your custom config directory
        instanovo_configs_dir = self.find_instanovo_configs()

        # Calculate relative path from InstaNovo's configs dir to your config dir
        relative_config_path = os.path.relpath(config_dir, instanovo_configs_dir)
        model_used = input_params["model_used"]
        if model_used not in ["transformer", "diffusion"]:
            logging.error(
                "Unknown search mode %s. Search mode has to be either "
                "'transformer' or 'diffusion'",
                model_used,
            )
        utrace.urun_dict.command_list = [
            str(self.exe_path),
            f"{model_used}",
            "predict",
            "--data-path",
            mgf_file,
            "--output-path",
            output_csv,
            "--config-path",
            relative_config_path,
            "--config-name",
            param_file.stem,
        ]
        return utrace

    def postflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Postflight routine for Instanovo_1 wrapper.
 
        During postflight the Instanovo native .csv output file is converted into the
        pre-defined urgap output file, which is of csv format.
 
        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.
 
        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        return utrace
