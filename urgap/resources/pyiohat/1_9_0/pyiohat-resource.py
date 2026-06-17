"""pyProstista resource."""

import argparse
import json
from pathlib import Path

# !/usr/bin/env python
from pyiohat import Unify


def main(
    search_result,
    fasta_file,
    spectrum_meta_data,
    xml_file_list,
    output_file,
    metadata_file,
    parameters,
):
    """Convert engine format to unified format using pyProtista-idents.

    Args:
        search_result (str, Path): path to engine output file
        fasta_file (str, Path): path to query database
        spectrum_meta_data (str, Path): path to meta data file
        xml_file_list (list): list of modification xml files
        output_file (str, Path): output file path
        parameters (dict): pyProtista parameter collection
    """
    params = json.loads(parameters)
    # params = {p["translated_key"]: p["translated_value"] for p in parameters.values()}
    params["database"] = fasta_file
    params["rt_pickle_name"] = spectrum_meta_data
    params["xml_file_list"] = [Path(f) for f in xml_file_list]
    u = Unify(search_result, params=params)
    df = u.get_dataframe()
    df.to_csv(output_file, index=False)
    metadata = u.get_metadata()
    with open(metadata_file, 'w', encoding='utf‑8') as f:
        json.dump(metadata, f, indent=2, sort_keys=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        help="input files to be aligned",
    )
    parser.add_argument(
        "-f",
        "--fasta_file",
        dest="fasta_file",
        help="fasta file",
    )
    parser.add_argument(
        "-md",
        "--spectrum_meta_data",
        dest="spectrum_meta_data",
        help="spectrum meta data file",
    )
    parser.add_argument(
        "-x",
        "--xml_file_list",
        dest="xml_file_list",
        action="append",
        default=[],
        help="xml file list",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        help="output file",
    )
    parser.add_argument(
        "-m",
        "--metadata_file",
        dest="metadata_file",
        help="metadata file",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        dest="parameters",
        help="json encoded parameters string",
    )

    args = parser.parse_args()

    main(
        search_result=args.input_file,
        fasta_file=args.fasta_file,
        spectrum_meta_data=args.spectrum_meta_data,
        xml_file_list=args.xml_file_list,
        output_file=args.output_file,
        metadata_file=args.metadata_file,
        parameters=args.parameters,
    )
