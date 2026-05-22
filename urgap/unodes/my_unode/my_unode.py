"""Your new unode goes here."""

import urgap


class MyUnode(urgap.unode.UNodeBase):
    """Urgap wrapper for <your tool here>."""

    META_INFO = {
        "name": "MyUnode",
        "wrapper_version": {"major": 0, "minor": 0, "patch": 1},
        "versions": [
            {
                "version": "0.0.1",
                "exe_path": "$my_unode_exe",
            },
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            # urgap.uftypes.my_uftype.SOME_INPUT: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            # urgap.uftypes.my_uftype.SOME_OUTPUT: {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("my_unode",),
        "citation": "",
        "parameter_examples": "",
    }

    def __init__(self) -> None:
        """Initialize MyUnode class."""
        super().__init__()

    def preflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Preflight routine: format parameters and compose the command list."""
        # TODO: build utrace.urun_dict.command_list
        return utrace

    def postflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Postflight routine: move output files to the right uftypes."""
        # TODO: utrace.move_output_files(files=[...], uftype=...)
        return utrace
