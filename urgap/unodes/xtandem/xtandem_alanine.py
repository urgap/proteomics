"""Urgap xtandem_alanine wrapper."""

import logging

from xml.etree import ElementTree as ET

import urgap


class XtandemAlanine(urgap.unode.UNodeBase):
    """Urgap wrapper for the xtandem_alanine search engine.

    X! Tandem is an open source software that can match tandem mass spectra with
    peptide sequences, in a process that has come to be known as protein identification.
    See publication provided under META_INFO["citation"] for further info.
    """

    META_INFO = {
        "name": "Xtandem_Alanine",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "1.0.0", "exe_path": "xtandem/alanine/tandem.exe"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.converter.PYMZML_MGF: {"min": 1, "max": 1},
            urgap.uftypes.proteomics.FASTA: {"min": 1, "max": 1},
            # NOTE: confirm this uftype name/location exists in your urgap install.
            urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS: {"min": 1, "max": 1},
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.dbsearch.XTANDEM_XML: {"min": 1, "max": 1},
        },
        "engine": None,
        "engine_type": ("identification",),
        "citation": """
        Craig, R., & Beavis, R. C. (2004). TANDEM: matching proteins with tandem mass spectra.
        In Bioinformatics (Vol. 20, Issue 9, pp. 1466-1467). Oxford University Press (OUP). https://doi.org/10.1093/bioinformatics/bth092
        """,
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize xtandem_alanine class."""
        super().__init__(*args, **kwargs)

    def preflight(self, utrace: urgap.UTrace) -> urgap.UTrace:
        """Preflight routine for xtandem_alanine wrapper.

        During preflight,
            - the supplied param file is validated and copied in as default_input.xml
            - taxonomy.xml and input.xml are generated with run-specific paths
            - the command list is composed

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.

        Raises:
            ValueError: If no param file (or more than one) was supplied, or if
                the param file does not specify "protein, taxon".
        """
        param_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.dbsearch.XTANDEM_PARAMS,
        )
        if len(param_files) != 1:
            msg = (
                "Xtandem_Alanine requires exactly one param file to be "
                "supplied via the XTANDEM_PARAMS input uftype."
            )
            raise ValueError(msg)

        mgf_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYMZML_MGF,
        )[0]
        fasta_file = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.FASTA,
        )[0]
        output_file = utrace.output_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.dbsearch.XTANDEM_XML,
        )[0]
        output_dir = output_file.parent

        xml_paths = {
            "default_input": output_dir / "default_input.xml",
            "taxonomy": output_dir / "taxonomy.xml",
            "input": output_dir / "input.xml",
        }

        # The user's param file is already valid X!Tandem bioml - copy it in as-is.
        with param_files[0].open() as fin, xml_paths["default_input"].open("w") as fout:
            param_file_content = fin.read()
            fout.write(param_file_content)
        logging.info(
            "Wrote input file %s (copied from supplied param file)",
            "default_input.xml",
        )

        # X!Tandem requires "protein, taxon" in default_input.xml to exactly
        # match a <taxon label="..."> entry in taxonomy.xml - read it from the
        # supplied param file rather than hardcoding a value here.
        taxon_root = ET.fromstring(param_file_content)
        taxon = None
        for note in taxon_root.findall("note"):
            if note.attrib.get("label") == "protein, taxon":
                taxon = (note.text or "").strip()
                break
        if not taxon:
            msg = (
                'Supplied param file must specify "protein, taxon" so it can '
                "be matched in the generated taxonomy.xml."
            )
            raise ValueError(msg)

        taxonomy_content = f"""<?xml version='1.0' encoding='iso-8859-1'?>
    <bioml label="x! taxon-to-file matching list">
    <taxon label="{taxon}">
    <file URL="{fasta_file}" format="peptide" />
    </taxon>
    </bioml>
    """
        with xml_paths["taxonomy"].open("w") as out:
            print(taxonomy_content, file=out)
            logging.info("Wrote input file %s", "taxonomy.xml")

        input_content = """<?xml version='1.0' encoding='iso-8859-1'?>
    <bioml>
    <note label="list path, default parameters" type="input">{default_input}</note>
    <note label="list path, taxonomy information" type="input">{taxonomy}</note>
    <note label="spectrum, path" type="input">{mgf_input_file}</note>
    <note label="output, path" type="input">{output_file_incl_path}</note>
        </bioml>""".format(
            default_input=xml_paths["default_input"],
            taxonomy=xml_paths["taxonomy"],
            mgf_input_file=mgf_file,
            output_file_incl_path=output_file,
        )
        with xml_paths["input"].open("w") as out:
            print(input_content, file=out)
            logging.info("Wrote input file %s", "input.xml")

        utrace.urun_dict.command_list = [
            str(self.exe_path),
            str(xml_paths["input"]),
        ]
        return utrace
