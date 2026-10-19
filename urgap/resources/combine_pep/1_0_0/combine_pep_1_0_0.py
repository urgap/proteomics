"""Combine per-engine PEP scores into unified Bayes and windowed PEP estimates."""

import csv
import functools
import itertools
import operator

from collections.abc import Iterator
from pathlib import Path


def adjust_window_size(desired_window_size: int, iter_len: int, minimum: int = 29) -> int:
    """Adjust the sliding window size based on the total number of values.

    When there are few values (below 1/5 of the window size), the window
    size is decreased.
    """
    if desired_window_size < iter_len // 5:
        adjusted_window_size = desired_window_size
    else:
        adjusted_window_size = desired_window_size // 5
        adjusted_window_size = max(adjusted_window_size, minimum)
    if adjusted_window_size != desired_window_size:
        pass
    return adjusted_window_size


def sliding_window_slow(
    iterable: list, window_size: int, flexible: bool = True,
) -> Iterator[tuple]:
    """Generate a sliding window over `iterable`.

    Slow but readable version using list slicing. Currently not used.
    """
    if flexible:
        window_size = adjust_window_size(window_size, len(iterable))

    if window_size % 2 == 0:
        window_size += 1

    half_window_size = int((window_size - 1) / 2)
    for center_i, center_value in enumerate(iterable):
        start_i = center_i - half_window_size
        start_i = max(start_i, 0)
        stop_i = center_i + half_window_size + 1
        yield iterable[start_i:stop_i], center_value


def sliding_window(
    elements: list, window_size: int, flexible: bool = True,
) -> Iterator[tuple]:
    """Generate a sliding window over `elements` without container types.

    Gives you sliding window functionality without using container types
    (list, deque etc.) to speed it up. Only works for lists of numbers.
    Yields the sum of all numbers in the sliding window (= the number of
    decoys in the sliding window in our case), the central number of the
    sliding window (required for the test only), and the current length
    of the sliding window (= total number of PSMs in the sliding window).
    Used for PEP calculation:
    PEP_of_PSM = (n_decoys_in_window * 2) / n_total_PSMs_in_window.
    """
    if flexible:
        window_size = adjust_window_size(window_size, len(elements))

    if window_size % 2 == 0:
        window_size += 1

    half_window_size = int((window_size - 1) / 2)

    start_gen, stop_gen = itertools.tee(elements)  # get 2 generators, one that tells
    # us which number to subtract from the back and another one that tells us
    # which number to add at the front of the sliding window

    n_decoys = 0  # keep track of the number of decoys in current sliding window
    current_win_size = 0  # keep track of current window size
    previous_start_i, previous_stop_i = 0, 0  # remember where our sliding window
    # was one iteration earlier (start and stop positions of the sliding window)

    for center_i, center_value in enumerate(elements):
        start_i = center_i - half_window_size
        start_i = max(start_i, 0)
        stop_i = center_i + half_window_size + 1

        if start_i != previous_start_i:
            # I have to substract a number! (the number behind (=left) of the
            # sliding window)
            n_decoys -= next(start_gen)
            current_win_size -= 1

        if stop_i != previous_stop_i:
            # I have to add a number! (the next number, the one in front
            # (=right) of the sliding window)
            for _i in range(stop_i - previous_stop_i):
                try:
                    n_decoys += next(stop_gen)
                    current_win_size += 1
                except StopIteration:
                    break  # cause StopIteration silently ends for-loops, will be fixed in py3.6 :)

        previous_start_i, previous_stop_i = start_i, stop_i
        yield n_decoys, center_value, current_win_size


class PSMDecoyStatusError(Exception):
    """Raised when a PSM's target/decoy status cannot be determined."""


class CombinedPEP:
    """Combine per-engine PEP scores into a unified Bayes PEP and windowed PEP."""

    def __init__(self) -> None:
        """Initialize empty PSM/score lookup structures."""
        self.psm_dicts = {}
        self.score_dict = {}
        self.input_csv_fieldnames = set()

    @staticmethod
    def row_is_decoy(row: dict) -> bool:
        """Check if a unified CSV row is a target or decoy PSM.

        Returns True if decoy and False if target.
        """
        if row["is_decoy"].lower() == "true":
            is_decoy = True
        elif row["is_decoy"].lower() == "false":
            is_decoy = False
        else:
            msg = "Could not determine whether PSM is decoy or not."
            raise PSMDecoyStatusError(msg)
        if "decoy_" in row.get("proteinacc_start_stop_pre_post_;", ""):
            is_decoy = True
        return is_decoy

    @staticmethod
    def naive_bayes(prob_list: list[float]) -> float:
        """Combine independent probabilities a, b, c using naive Bayes.

                                a*b*c
        combined_prob = -------------------------
                        a*b*c + (1-a)*(1-b)*(1-c)

        For a straightforward explanation, see
        http://www.paulgraham.com/naivebayes.html
        """
        # multiplying all probabilities: a*b*c
        multiplied_probs = functools.reduce(operator.mul, prob_list, 1)
        # multiplying the opposite of all probabilities: (1-a)*(1-b)*(1-c)
        multiplied_opposite_probs = functools.reduce(
            operator.mul, (1 - p for p in prob_list), 1,
        )
        return multiplied_probs / (multiplied_probs + multiplied_opposite_probs)

    def add_engine_result_csv(self, input_csv_path: str, input_engine_name: str) -> None:
        """Parse one engine's result CSV and buffer it as a PSM-to-row dict."""
        with Path(input_csv_path).open(encoding="utf8") as csv_obj:
            reader = csv.DictReader(csv_obj)
            for fieldname in reader.fieldnames:
                self.input_csv_fieldnames.add(fieldname)
            self.psm_dicts[input_engine_name] = self.reader_to_psm_dict(reader)

    def reader_to_psm_dict(self, reader: csv.DictReader) -> dict:
        """Turn a CSV reader object into a dictionary.

        Dict key is a tuple of column values that can be specified by the
        user (usually a combination of Seq/Mods/Spectrum). Dict value is
        the DictReader row.
        """
        psm_dict = {}
        for row in reader:
            row_key = self.get_row_key(row, self.columns_for_grouping)
            psm_dict[row_key] = row
        return psm_dict

    @staticmethod
    def list_to_sorted_tuple(values: list) -> tuple:
        """Return a sorted, de-duplicated tuple built from `values`."""
        return tuple(sorted(set(values)))

    @staticmethod
    def get_row_key(row: dict, columns_for_grouping: tuple) -> tuple:
        """Build the grouping key for a row from the configured columns."""
        return tuple(row[c] for c in columns_for_grouping)

    @staticmethod
    def all_combinations(some_iterable: list) -> Iterator[tuple]:
        """Yield every non-empty combination of `some_iterable`, smallest first."""
        for i in range(len(some_iterable)):
            for combo in itertools.combinations(some_iterable, i + 1):
                yield tuple(sorted(combo))

    def generate_psm_to_scores_dict(self, input_engines: list[str]) -> None:
        """Compute Bayes PEP and windowed combined PEP for every engine combination."""
        for engine_combo in self.all_combinations(input_engines):

            engines_not_in_combo = {e for e in input_engines if e not in engine_combo}

            if engines_not_in_combo:
                pass
            else:
                pass

            self.score_dict[engine_combo] = {}

            # get all shared PSMs:
            all_psms_of_combo_engines = set.intersection(
                *(set(self.psm_dicts[engine].keys()) for engine in engine_combo),
            )

            # remove all PSMs that are found by other engines:
            for other_eng in engines_not_in_combo:
                all_psms_of_combo_engines -= set(self.psm_dicts[other_eng].keys())

            if not all_psms_of_combo_engines:  # nothing to do here...
                continue

            # For every PSM, use naive Bayes to calculate the combined PEP
            # ('Bayes PEP') among all engines and add it to self.score_dict:

            decoy_count_of_intersection = 0
            for _psm_count_of_intersection, psm_key in enumerate(
                all_psms_of_combo_engines, start=1,
            ):
                engine_peps = []
                engine_decoy_bools = set()
                for eng in engine_combo:
                    d = self.psm_dicts[eng]
                    engine_peps.append(float(d[psm_key]["posterior_error_prob"]))
                    engine_decoy_bools.add(self.row_is_decoy(d[psm_key]))

                if len(engine_decoy_bools) != 1:
                    psm_is_decoy = True
                else:
                    psm_is_decoy = next(iter(engine_decoy_bools))

                if psm_is_decoy:
                    decoy_count_of_intersection += 1

                bayes_pep = self.naive_bayes(engine_peps)

                self.score_dict[engine_combo][psm_key] = {
                    "Bayes PEP": bayes_pep,
                    "Is decoy": psm_is_decoy,
                }


            # Use the Bayes PEP to sort all PSMs that are shared by the current
            # combination of engines. Compute the PEP (similar to localized FDR)
            # for each PSM and store it in score_dict:

            psms_sorted_by_bayes_pep = sorted(
                self.score_dict[engine_combo].items(), key=lambda kv: kv[1]["Bayes PEP"],
            )
            sorted_decoy_bools = [
                kv_tuple[1]["Is decoy"] for kv_tuple in psms_sorted_by_bayes_pep
            ]

            for i, (n_decoys, __, current_win_size) in enumerate(
                sliding_window(sorted_decoy_bools, self.window_size),
            ):

                n_false_positives = 2 * n_decoys
                n_false_positives = min(n_false_positives, current_win_size)
                intersection_pep = n_false_positives / current_win_size
                current_psm_key = psms_sorted_by_bayes_pep[i][0]
                self.score_dict[engine_combo][current_psm_key][
                    "combined PEP"
                ] = intersection_pep

            assert i + 1 == len(sorted_decoy_bools), (
                "sliding_window() "
                "did not return a sliding window for all PSMs! This "
                "should never happen!"
            )

    def write_output_csv(self, output_csv_path: str) -> None:
        """Write the merged CSV with added Bayes PEP and combined PEP columns."""
        new_scores = ["combined PEP", "Bayes PEP"]
        fieldnames = (
            list(self.input_csv_fieldnames) + new_scores + ["combined PEP engines"]
        )
        with Path(output_csv_path).open("w", encoding="utf8") as out_obj:
            writer = csv.DictWriter(out_obj, fieldnames=fieldnames)
            writer.writeheader()

            for engine_combo, combo_score_dict in self.score_dict.items():

                for psm_key, score_dict_val in combo_score_dict.items():

                    psm_rows_from_all_engines = [
                        self.psm_dicts[engine][psm_key] for engine in engine_combo
                    ]
                    for out_row in psm_rows_from_all_engines:
                        # add column that lists all engines that found the PSM:
                        out_row["combined PEP engines"] = self.join_sep.join(
                            engine_combo,
                        )
                        # add columns with Bayes PEP and combined PEP:
                        for score_field in new_scores:
                            out_row[score_field] = score_dict_val[score_field]
                        writer.writerow(out_row)


def main(
    columns_for_grouping: list[str] | None = None,
    input_csvs: list[str] | None = None,
    output_csv: str | None = None,
    input_sep: str | None = None,
    output_sep: str | None = None,
    join_sep: str | None = None,
    pep_colname: str | None = None,
    input_engines: list[str] | None = None,
    window_size: int | None = None,
) -> None:
    """Combine per-engine PEP scores into a merged CSV with Bayes and combined PEP.

    Steps performed:
        1. Set parsed attributes from command line as class attributes.
        2. Parse unified input CSV files (with PEPs) and buffer them as
           PSM-to-row dictionaries.
        3. Calculate the Bayes PEP for each PSM.
        4. Retrieve all possible engine combinations (like a Venn diagram).
        5. For each engine combination:
            a) retrieve list of shared PSMs
            b) sort PSMs by Bayes PEP
            c) loop over sorted PSMs and calculate combined PEP (sliding window)
        6. Write merged CSV output file with added Bayes PEP and combined
           PEP columns.
    """
    c = CombinedPEP()

    # Set parsed attributes from command line as class attribute:
    c.columns_for_grouping = c.list_to_sorted_tuple(columns_for_grouping)
    c.input_sep = input_sep
    c.output_sep = output_sep
    c.join_sep = join_sep
    c.pep_colname = pep_colname
    c.window_size = window_size

    # Parse input CSV files and buffer them as PSM-to-row dicts:
    for input_csv, input_engine in zip(input_csvs, input_engines, strict=False):
        c.add_engine_result_csv(input_csv, input_engine)

    # Calculate Bayes PEP and combined PEP for each PSM:
    c.generate_psm_to_scores_dict(input_engines)

    # Write merged CSV with Bayes PEP and combined PEP column:
    c.write_output_csv(output_csv)


def parse_args(verbose: bool = True) -> dict:
    """Parse command line arguments and return them as a dict.

    Only used when executing this script from command line.
    """
    import argparse

    # no need to import argparse when script is executed by importing
    # main function.

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_csvs",
        type=str,
        nargs="+",
        required=True,
        help="Paths to unified input CSV files (2 or more)",
    )
    parser.add_argument(
        "-c",
        "--columns_for_grouping",
        type=str,
        nargs="+",
        required=True,
        help="Column names by which the rows should be grouped",
    )
    parser.add_argument("-o", "--output_csv", type=str, help="Output CSV name")
    parser.add_argument(
        "-is",
        "--input_sep",
        type=str,
        default=",",
        help="Input file column delimiter character",
    )
    parser.add_argument(
        "-os",
        "--output_sep",
        type=str,
        default=",",
        help="Output file column delimiter character",
    )
    parser.add_argument(
        "-js",
        "--join_sep",
        type=str,
        default=";",
        help="Delimiter for multiple values in the same field",
    )
    parser.add_argument(
        "-e",
        "--input_engines",
        type=str,
        nargs="+",
        required=True,
        help="The search engines of each input file (must be same order as input_csvs)",
    )
    parser.add_argument(
        "-w",
        "--window_size",
        type=int,
        default=251,
        help="The size of the sliding window for PEP calculation.",
    )

    args = vars(parser.parse_args())  # convert to dict
    if verbose:
        for _arg, _val in sorted(args.items()):
            pass
    return args


if __name__ == "__main__":
    command_line_args = parse_args()
    main(**command_line_args)
