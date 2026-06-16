"""pymzml resource."""

import argparse
import logging

from pathlib import Path
from typing import TYPE_CHECKING

import pymzml

if TYPE_CHECKING:
    from io import TextIOWrapper

    from pymzml.run import Reader
    from pymzml.spec import Spectrum


def _parse_rt_unit(
    spec: "Spectrum",
    scan_time: float,
    was_warned: bool,
) -> tuple[float, bool]:
    """Calculate normalized retention time in seconds based on metadata."""
    _, scan_time_unit = spec.scan_time
    unit_upper = scan_time_unit.upper()

    if unit_upper in ["MINUTE", "MINUTES"]:
        return scan_time * 60, was_warned

    if unit_upper not in ["SECOND", "SECONDS"]:
        if not was_warned:
            logging.warning(
                "[Warning] The retention time unit is not recognized "
                "or not specified. It is assumed to be minutes and "
                "continues with that.",
            )
            was_warned = True
        return scan_time * 60, was_warned

    return scan_time, was_warned


def _filter_spectrum(
    spec_ms_level: int,
    ms_level: int,
    spectrum_id: int,
    scan_inclusion_list: list[int] | None,
    scan_exclusion_list: list[int],
) -> bool:
    """Evaluate filtration matrices to confirm scan processing eligibility."""
    if spec_ms_level != ms_level:
        return False
    if scan_inclusion_list is not None and spectrum_id not in scan_inclusion_list:
        return False
    return spectrum_id not in scan_exclusion_list


def _process_peaks(
    peaks: list[tuple[float, float]],
    mz_correction_factor: float,
    number_of_mz_decimals: int,
    number_of_i_decimals: int,
) -> tuple[list[str], bool]:
    """Format and offset-correct individual isotopic mass peaks."""
    mz_i_list = []
    for raw_mz, intensity in peaks:
        corrected_mz = raw_mz + (raw_mz * mz_correction_factor)
        mz_i_list.append(
            f"{corrected_mz:<10.{number_of_mz_decimals}f} "
            f"{intensity:<10.{number_of_i_decimals}f}",
        )
    return mz_i_list, len(mz_i_list) > 0


def _build_charge_range(
    min_charge: int | None,
    max_charge: int | None,
) -> str | None:
    """Construct printable formatted strings indicating spectrum limits."""
    if min_charge is not None and max_charge is not None:
        precursor_charge_range = f"{min_charge}"
        for charge in range(min_charge + 1, max_charge + 1):
            precursor_charge_range += f" and {charge}"
        return precursor_charge_range
    return None


def _write_ion_entry(
    fout: "TextIOWrapper",
    mgf_name_base: str,
    spectrum_id: int,
    precursor_charge: int | None,
    precursor_mz: float,
    scan_time: float,
    ion_mode: str,
    precursor_charge_range: str | None,
    mz_i_string: str,
) -> None:
    """Construct and serialize the text block representation of a mass scan."""
    if ion_mode in ["pos", "+"]:
        ion_sign = "+"
    elif ion_mode in ["neg", "-"]:
        ion_sign = "-"
    else:
        ion_sign = ion_mode
        logging.warning("[Warning] Unknown ion mode: %s", ion_mode)

    if precursor_charge is not None:
        c_string = f"CHARGE={precursor_charge}{ion_sign}"
    elif precursor_charge_range is not None:
        c_string = f"CHARGE={precursor_charge_range}{ion_sign}"
    else:
        c_string = "CHARGE="

    fout.write(
        f"BEGIN IONS\n"
        f"TITLE={mgf_name_base}.{spectrum_id}.{spectrum_id}."
        f"{precursor_charge}\n"
        f"SCANS={spectrum_id}\n"
        f"RTINSECONDS={round(scan_time, 11)}\n"
        f"PEPMASS={precursor_mz}\n"
        f"{c_string}\n"
        f"{mz_i_string}\n"
        f"END IONS\n\n",
    )


def main(
    mzml: str | Path | None = None,
    mgf: str | Path | None = None,
    number_of_i_decimals: int = 5,
    number_of_mz_decimals: int = 5,
    machine_offset_in_ppm: float | None = None,
    scan_exclusion_list: list[str | int] | None = None,
    scan_inclusion_list: list[int] | None = None,
    scan_skip_modulo_step: int | None = None,
    ms_level: int = 2,
    precursor_min_charge: int = 1,
    precursor_max_charge: int = 5,
    ion_mode: str = "pos",
    signal_to_noise_threshold: float | None = None,
    **kwargs: str | float | None,
) -> Path:
    """Convert mzML to mgf.

    A new mgf file will be created at the same location as the mzML file

    Usage:
    ./pymzml2mgf.py <mzML_file_name> <mgf_file_name>
    """
    _ = kwargs
    logging.info("Converting file:\n\tmzml : %s\n\tto\n\tmgf : %s", mzml, mgf)
    if not mzml or not mgf:
        err_msg = "Both mzml and mgf paths must be provided."
        raise ValueError(err_msg)

    mzml_path = Path(mzml)
    mgf_path = Path(mgf)
    mgf_name_base = mgf_path.name.split(".mgf")[0]
    run: Reader = pymzml.run.Reader(str(mzml_path))

    mgf_entries = 0
    written_specs = 0
    exclusion_ints = (
        [int(spec_id) for spec_id in scan_exclusion_list] if scan_exclusion_list else []
    )

    mz_correction_factor = machine_offset_in_ppm * 1e-6 if machine_offset_in_ppm else 0.0

    precursor_charge_range = _build_charge_range(
        precursor_min_charge,
        precursor_max_charge,
    )

    mzml_basename = mzml_path.name
    was_warned = False

    with mgf_path.open("w") as fout:
        for n, raw_spec in enumerate(run):
            if n % 5000 == 0:
                logging.info(
                    "File : %40s : Processing spectrum %d",
                    mzml_basename,
                    n,
                )

            scan_time, _ = raw_spec.scan_time
            scan_time, was_warned = _parse_rt_unit(
                raw_spec,
                scan_time,
                was_warned,
            )
            spectrum_id = int(raw_spec.ID)

            if not _filter_spectrum(
                raw_spec.ms_level,
                ms_level,
                spectrum_id,
                scan_inclusion_list,
                exclusion_ints,
            ):
                continue

            working_spec = (
                raw_spec.remove_noise(
                    mode="median",
                    signal_to_noise_threshold=signal_to_noise_threshold,
                )
                if signal_to_noise_threshold is not None
                else raw_spec
            )

            peaks_2_write = working_spec.peaks("centroided")
            precursor_mz = working_spec.selected_precursors[0]["mz"]
            precursor_charge = working_spec.selected_precursors[0].get(
                "charge",
                None,
            )
            precursor_mz += precursor_mz * mz_correction_factor

            if len(peaks_2_write) == 0:
                continue

            if (
                scan_skip_modulo_step is not None
                and mgf_entries % scan_skip_modulo_step != 0
            ):
                mgf_entries += 1
                continue

            mgf_entries += 1

            mz_i_list, valid_peaks = _process_peaks(
                peaks_2_write,
                mz_correction_factor,
                number_of_mz_decimals,
                number_of_i_decimals,
            )
            if not valid_peaks:
                continue
            mz_i_string = "\n".join(mz_i_list)

            _write_ion_entry(
                fout,
                mgf_name_base,
                spectrum_id,
                precursor_charge,
                precursor_mz,
                scan_time,
                ion_mode,
                precursor_charge_range,
                mz_i_string,
            )
            written_specs += 1

    logging.info("Wrote %d mgf entries to file %s", written_specs, mgf_path)
    return mgf_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        dest="input_file",
        help="input file to be converted",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        dest="output_file",
        help="output file",
    )
    parser.add_argument(
        "-ni",
        "--num_i_decimals",
        type=int,
        dest="num_i_decimals",
        help="number of decimal points for intensity",
        default=5,
    )
    parser.add_argument(
        "-nm",
        "--num_mz_decimals",
        type=int,
        dest="num_mz_decimals",
        help="number of decimal points for mz",
        default=5,
    )
    parser.add_argument(
        "-of",
        "--offset",
        type=int,
        dest="offset",
        help="machine offset in ppm",
        default=0,
    )
    parser.add_argument(
        "-el",
        "--exclusion_list",
        action="append",
        dest="exclusion_list",
        help="list of excluded spectra",
        default=None,
    )
    parser.add_argument(
        "-il",
        "--inclusion_list",
        action="append",
        dest="inclusion_list",
        help="list of excluded spectra",
        default=None,
    )
    parser.add_argument(
        "-s",
        "--skip_step",
        type=int,
        dest="skip_step",
        help="scan skip modulo step",
        default=None,
    )
    parser.add_argument(
        "-ms",
        "--ms_level",
        type=int,
        dest="ms_level",
        help="ms level",
        default=2,
    )
    parser.add_argument(
        "-cmin",
        "--min_charge",
        type=int,
        dest="min_charge",
        help="minimum precursor charge",
        default=1,
    )
    parser.add_argument(
        "-cmax",
        "--max_charge",
        type=int,
        dest="max_charge",
        help="maximum precursor charge",
        default=5,
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        dest="mode",
        help="ion mode",
        default="pos",
    )
    parser.add_argument(
        "-sn",
        "--signal_noise_threshold",
        type=float,
        dest="signal_noise_threshold",
        help="signal to noise threshold",
        default=None,
    )
    args = parser.parse_args()
    main(
        mzml=args.input_file,
        mgf=args.output_file,
        number_of_i_decimals=args.num_i_decimals,
        number_of_mz_decimals=args.num_mz_decimals,
        machine_offset_in_ppm=args.offset,
        scan_exclusion_list=args.exclusion_list,
        scan_inclusion_list=args.inclusion_list,
        scan_skip_modulo_step=args.skip_step,
        ms_level=args.ms_level,
        precursor_min_charge=args.min_charge,
        precursor_max_charge=args.max_charge,
        ion_mode=args.mode,
        signal_to_noise_threshold=args.signal_noise_threshold,
    )
