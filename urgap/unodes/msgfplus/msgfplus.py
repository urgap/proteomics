"""Urgap MsgfPlus wrapper."""


import multiprocessing

from pathlib import Path

import urgap


class UnrecognizedFragmentationMethodError(ValueError):
    """Raised when a FragmentationMethodID is not recognized."""

    def __init__(self, method_id: str) -> None:
        """Initialize the error with the unrecognized FragmentationMethodID."""
        super().__init__(f"FragmentationMethodID '{method_id}' is not a recognized value.")

class MsgfPlus(urgap.unode.UNodeBase):
    """Urgap wrapper for the MsgfPlus search engine.

    MS-GF+ (aka MSGF+ or MSGFPlus) performs peptide identification by scoring MS/MS
    spectra against peptides derived from a protein sequence database. MS-GF+ is
    optimized for a variety of spectral types, i.e., combinations of fragmentation
    method, instrument, enzyme, and experimental protocols. See publication provided
    under META_INFO["citation"] for further info.

    Requires:
        java
        unimod_mapper
    """

    META_INFO = {
        "name": "MsgfPlus",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "20240326", "exe_path": "MsgfPlus/2024_03_26/MSGFPlus.jar"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.converter.PYMZML_MGF: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.FASTA: {"min": 1, "max": 1},
            # urgap.uftypes.proteomics.MODS_XML: {"min": 0, "max": -1},
            urgap.uftypes.proteomics.params.MSGFPLUS_TXT: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID: {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("db_search", "proteomics"),
        "citation": """
        Kim, S., Mischerikow, N., Bandeira, N., Navarro, J. D., Wich, L., Mohammed, S., Heck, A. J. R., & Pevzner, P. A. (2010). The Generating Function of CID, ETD, and CID/ETD Pairs of Tandem Mass Spectra: Applications to Database Search.
        In Molecular Cellular Proteomics (Vol. 9, Issue 12, pp. 2840-2852). Elsevier BV. https://doi.org/10.1074/mcp.m110.003731
        """,
    }

    def __init__(self, *args: str, **kwargs: str) -> None:
        """Initialize MsgfPlus class."""
        super().__init__(*args, **kwargs)
        self.mgf_new_input_file = None


    def reformat_mgf_input(
        self,
        utrace: urgap.UTrace,
    ) -> None:
        """Reformat the mgf input file, to be injested into the msgfplus search.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.
        """
        param_to_activation_type = {
            0: "CID",
            1: "CID",
            2: "ETD",
            3: "HCD",
        }
        conf_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.params.MSGFPLUS_TXT,
        )[0]
        mgf_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.converter.PYMZML_MGF,
        )[0]
        with Path.open(conf_file) as file:
            for line in file:
                if line.strip().startswith("FragmentationMethodID="):
                    id_string = line.strip().split("=")[1]
                    method_id = int(id_string)
                    if method_id in param_to_activation_type:
                        activation_type = param_to_activation_type[method_id]
                    else:
                        raise UnrecognizedFragmentationMethodError(method_id)
        with mgf_file.open(encoding="UTF-8") as mgf_org_input_file:
            lines = mgf_org_input_file.readlines()
        self.mgf_new_input_file = str(mgf_file.parent / mgf_file.stem) + "_tmp.mgf"
        with Path(self.mgf_new_input_file).open("w", encoding="UTF-8") as mgf_new_input_file:
            for line in lines:
                if line.startswith("CHARGE"):
                    print(line, file=mgf_new_input_file)
                    print(
                        f"ACTIVATIONMETHOD={activation_type}",
                        file=mgf_new_input_file,
                    )
                else:
                    print(line, file=mgf_new_input_file)
        self.tmp_files.append(self.mgf_new_input_file)

    def create_command_list(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Create the command list from input parameters.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        input_params = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]

        utrace.urun_dict.command_list = [
            "java",
            "-jar",
            str(self.exe_path),
        ]
        clist = utrace.urun_dict.command_list

        for key, value in input_params.items():
            if key == "-thread":
                if value == "max - 1":

                    cpu_value = multiprocessing.cpu_count() - 1
                    clist.extend((key, cpu_value))
                else:
                    clist.extend((key, value))
            elif key == "-Xmx":
                clist.insert(
                    1,
                    f"{key}{value}",
                )

        mgf_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.converter.PYMZML_MGF,
        )[0]
        fasta_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.FASTA,
        )[0]
        output_file = utrace.output_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.dbsearch.MSGFPLUS_MZID,
        )[0]
        conf_file = utrace.input_files.get_path_objects_by_uftype(
            uftype=urgap.uftypes.proteomics.params.MSGFPLUS_TXT,
        )[0]

        clist.extend(
            [
                "-s",
                str(mgf_file),
                "-d",
                str(fasta_file),
                "-o",
                str(output_file),

                "-conf",
                str(conf_file),
            ],
        )
        for key, value in input_params.items():
            if key not in ["-thread", "-Xmx"]:
                clist.extend([key, value])
        return utrace

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for MsgfPlus wrapper.

        During preflight,
            - parameters are formatted
            - mods are mapped
            - helper files are written
            - command list is composed

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        self.reformat_mgf_input(utrace=utrace)
        utrace = self.create_command_list(utrace=utrace)

        enzyme_txt_path = utrace.output_files[0].path.parent / "params" / "enzymes.txt"
        if enzyme_txt_path.exists() is False:
            (utrace.output_files[0].path.parent/"params").symlink_to(
              Path(self.exe_path).parent / "Docs" / "Examples",
            )

        self.tmp_files.append(utrace.output_files[0].path.parent / "params")
        return utrace
