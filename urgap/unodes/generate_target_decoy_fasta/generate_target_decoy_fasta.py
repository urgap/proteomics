"""Urgap generate_target_decoy_fasta wrapper."""

import urgap


class GenerateTargetDecoyFasta(urgap.unode.UNodeBase):
    """Urgap wrapper for the generate_target_decoy_fasta resource.

    Allows to create a target decoy database from one (or more) fasta files.
    """

    META_INFO = {
        "name": "GenerateTargetDecoyFasta",
        "wrapper_version": {"major": 2, "minor": 0, "patch": 0},
        "versions": [
            {
                "version": "2.0.0",
                "exe_path": "generate_target_decoy_fasta_2_0_0/generate_target_decoy_fasta_2_0_0.py",
            },
        ],
        "parameters_not_triggering_rerun": [],
        "engine_type": ("proteomics"),
        "engine": None,
        "requires": {
            "other_uftypes": {
                "python_packages": [
                    "pyahocorasick",
                ],
            },
        },
        "input_uftypes": {
            urgap.uftypes.proteomics.FASTA: {
                "min": 1,
                "max": -1,
            },
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.FASTA: {
                "min": 1,
                "max": 1,
            },
            urgap.uftypes.ms.IMMUTABLE_PEPTIDES: {
                "min": 1,
                "max": 1,
            },
        },
        "citation": "Kremer, L. P. M., Leufken, J., Oyunchimeg, P., Schulze, S., & Fufezan, C. (2016). Ursgal, Universal Python Module Combining Common Bottom-Up Proteomics Tools for Large-Scale Analysis."
        " Journal of Proteome Research, 15(3), 788-794. https://doi.org/10.1021/acs.jproteome.5b00860 ",
    }

    def __init__(self, *args: str, **kwargs: str) -> None:
        """Initialize generate_target_decoy_fasta class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for generate_target_decoy_fasta wrapper.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        utrace.urun_dict.command_list = [
            "python",
            str(self.exe_path),
            "--output",
            str(
                utrace.output_files.get_path_objects_by_uftype(
                    urgap.uftypes.proteomics.FASTA,
                )[0],
            ),
            "--immutable_file",
            str(
                utrace.output_files.get_path_objects_by_uftype(
                    urgap.uftypes.ms.IMMUTABLE_PEPTIDES,
                )[0],
            ),

        ]
        for parameter_key, parameter_value in utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ].items():
            if parameter_value is not None:
                utrace.urun_dict.command_list.extend([parameter_key, parameter_value])
        if len(utrace.input_files) > 0:
            utrace.urun_dict.command_list.append("--input")
        for ufile in utrace.input_files:
            utrace.urun_dict.command_list.append(str(ufile.path))
        return utrace
