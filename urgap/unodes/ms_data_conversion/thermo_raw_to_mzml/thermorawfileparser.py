"""Urgap ThermoRawFileParser wrapper."""

import sys

import urgap


class ThermoRawFileParser(urgap.unode.UNodeBase):
    """Urgap wrapper for the ThermoRawFileParser executable.

    ThermoRawFileParser can be used to convert Thermo.raw files to open XML-based
    format, mzml,  for encoding mass spectrometer data. See publication provided under
    META_INFO["citation"] for further info.

    The exe can be downloaded from:
        https://github.com/CompOmics/ThermoRawFileParser/releases/tag/v.2.0.0-dev

    Requires:
        Mono
    """

    META_INFO = {
        "name": "ThermoRawFileParser",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {
                "version": "2.0.0",
                "exe_path": "ms_data_conversion/thermorawfileparser/2.0.0/ThermoRawFileParser.exe",
            }
        ],
        "parameters_not_triggering_rerun": [],
        "engine_type": ("converter",),
        "engine": None,
        "input_uftypes": {
            urgap.uftypes.proteomics.THERMO_RAW: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML: {
                "min": 1,
                "max": 1,
            },
        },
        "citation": """
        Hulstaert, N., Shofstahl, J., Sachsenberg, T., Walzer, M., Barsnes, H., Martens, L., & Perez-Riverol, Y. (2019). ThermoRawFileParser: Modular, Scalable, and Cross-Platform RAW File Conversion.
        In Journal of Proteome Research (Vol. 19, Issue 1, pp. 537-542). American Chemical Society (ACS). https://doi.org/10.1021/acs.jproteome.9b00328
        """,
    }

    def __init__(self, *args: str, **kwargs: str):
        """Initialize ThermoRawFileParser class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for ThermoRawFileParser wrapper.

        During preflight,
            - params are extracted from urun_dict
            - command list is composed

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        if sys.platform in ["win32"]:
            utrace.urun_dict.command_list = []
        else:
            utrace.urun_dict.command_list = ["mono"]

        input_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.THERMO_RAW
        )[0]
        output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.ms.converter.mzml.THERMORAWPARSER_MZML
        )[0]
        utrace.urun_dict.command_list.extend(
            [
                str(self.exe_path),
                f"-i={input_file!s}",
                f"-b={output_file!s}",
                "-f=1",
            ]
        )
        params = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]
        for key, value in params.items():
            if value is not None:
                utrace.urun_dict.command_list.append(f"{key}={value}")
        return utrace
