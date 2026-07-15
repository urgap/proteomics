"""urgap percolator_3_7_1 wrapper."""

import os
import shutil

from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd

from chemical_composition.chemical_composition_kb import PROTON

import urgap


class PercolatorEngineMismatchError(Exception):
    """Raised when input PSMs originate from more than one search engine."""


def _apply_parallel(
    df_grouped: pd.core.groupby.generic.DataFrameGroupBy,
    func: object,
    threads: int = -1,
) -> pd.DataFrame:
    """Apply a function to grouped dataframe rows across worker processes.

    Args:
        df_grouped: Result of a pandas groupby call.
        func: Function to apply to each group.
        threads: Number of worker processes; -1 uses all available CPUs.

    Returns:
        Concatenated dataframe of all processed groups.
    """
    if threads == -1:
        threads = cpu_count()
    with Pool(threads) as p:
        ret_list = p.map(func, [group for _, group in df_grouped])
    return pd.concat(ret_list)


class Percolator(urgap.unode.UNodeBase):
    """urgap wrapper for the percolator_3_7_1 executable.

    Percolator uses a semi-supervised machine learning to discriminate correct from
    incorrect peptide-spectrum matches, and calculates accurate statistics such as
    q-value (FDR) and posterior error probabilities. See publication
    provided under META_INFO["citation"] for further info.
    """

    META_INFO = {
        "name": "Percolator",
        "wrapper_version": {"major": 1, "minor": 0, "patch": 0},
        "versions": [
            {"version": "3.7.1", "exe_path": "percolator/3_7_1/percolator"},
        ],
        "parameters_not_triggering_rerun": [],
        "input_uftypes": {
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV: {
                "min": 1,
                "max": -1,
            },
            urgap.uftypes.proteomics.FASTA: {
                "min": 0,
                "max": 1,
            },
        },
        "output_uftypes": {
            urgap.uftypes.proteomics.validator.PERCOLATOR_CSV: {"min": 1, "max": 2},
        },
        "engine_type": ("validation", "proteomics"),
        "citation": """
        The, M., MacCoss, M. J., Noble, W. S., & Kall, L. (2016). Fast and Accurate Protein False Discovery Rates on Large-Scale Proteomics Data Sets with Percolator 3.0.
        In Journal of the American Society for Mass Spectrometry (Vol. 27, Issue 11, pp. 1719-1727). American Chemical Society (ACS). https://doi.org/10.1007/s13361-016-1460-7
        """,
    }

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize percolator_3_7_1 class."""
        super().__init__(*args, **kwargs)

    def preflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Preflight routine for percolator_3_7_1 wrapper.

        During preflight,
            - input file aligned with percolator style is formatted and created
            - command list is composed

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        self.output_type_dict = utrace.output_files.get_index_groups_by_uftypes()
        psm_file_indices = utrace.output_files.get_indices_by_uftype(
            urgap.uftypes.proteomics.validator.PERCOLATOR_CSV,
        )
        self.result_psms = (
            str(utrace.output_files[psm_file_indices[0]].path) + "targets_broken"
        )

        self.decoy_psms = (
            str(utrace.output_files[psm_file_indices[0]].path) + "decoys_broken"
        )

        input_tsv = self.create_input_file(utrace)
        utrace.urun_dict.command_list = [
            self.exe_path,
            "--only-psms",
            input_tsv,
            "--results-psms",
            self.result_psms,
            "--decoy-results-psms",
            self.decoy_psms,
        ]
        return self.create_command_list(utrace)

    def postflight(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Postflight routine for percolator_3_7_1 wrapper.

        During postflight the individual fixed and decoy dataframes coming from the
        percolator tool, are read and merged together before they are stored into the
        pre-defined urgap output file.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        fixed_targets_path = utrace.output_files[0].path.parent / Path(
            utrace.output_files[0].path.name + "_percolator_out_fixed.tsv",
        )

        fixed_decoys_path = utrace.output_files[0].path.parent / Path(
            utrace.output_files[0].path.name + "_percolator_out_fixed_decoys.tsv",
        )

        self.tmp_files.append(self.result_psms)
        self.tmp_files.append(self.decoy_psms)
        self.tmp_files.append(fixed_targets_path)
        self.tmp_files.append(fixed_decoys_path)

        for broken_path, fixed_path in [
            (self.decoy_psms, fixed_decoys_path),
            (self.result_psms, fixed_targets_path),
        ]:
            with Path(broken_path).open() as fin, Path(fixed_path).open("w") as fout:
                for i, line in enumerate(fin):
                    if i == 0:
                        fout.write(line)
                        continue
                    line_fields = line.split("\t")[:5]
                    prot = " ".join(line.split("\t")[5:])
                    line_fields.append(prot)
                    fout.write("\t".join(line_fields))
        # rename files again

        output_decoys = pd.read_csv(fixed_decoys_path, sep="\t", index_col=False)
        output_decoys = output_decoys[["PSMId", "q-value", "posterior_error_prob"]]

        output_targets = pd.read_csv(fixed_targets_path, sep="\t", index_col=False)
        output_targets = output_targets[["PSMId", "q-value", "posterior_error_prob"]]

        qvals = pd.concat([output_targets, output_decoys])

        unified_df = pd.read_csv(self.merged_frame)

        final_df = unified_df.merge(
            qvals, left_on="PSMId", right_on="PSMId", how="left",
        )
        final_df = final_df[~final_df["q-value"].isna()]
        idx = self.output_type_dict[".percolator.csv"][0]
        final_df.to_csv(utrace.output_files[idx].path)

        # Part specific for only version 3.7.1
        if self.META_INFO["unode_version"] == "3.7.1" and (
            utrace.output_files[0].path.parent / "target_protein_qvals.tsv"
        ).exists():
            utrace.extend_output_files_by_uftype(
                urgap.uftypes.proteomics.validator.PERCOLATOR_CSV,
            )
            protein_targets = (
                utrace.output_files[0].path.parent / "target_protein_qvals.tsv"
            )
            protein_decoys = (
                utrace.output_files[0].path.parent / "decoy_protein_qvals.tsv"
            )
            targets = pd.read_csv(
                protein_targets,
                sep="\t",
            )
            decoys = pd.read_csv(
                protein_decoys,
                sep="\t",
            )
            self.tmp_files.extend([protein_decoys, protein_targets])
            td_df = pd.concat([targets, decoys]).sort_values("ProteinGroupId")
            output_path = utrace.output_files[1]
            td_df.to_csv(str(output_path.path), index=False)

        return utrace

    def create_command_list(
        self,
        utrace: urgap.UTrace,
    ) -> urgap.UTrace:
        """Create the command list to execute percolator executable.

        Based on the input parameters, the command list is created, which will be used
        during execute step to run percolator.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            UTrace object, combination of urun_dict, ufile_list and unode.meta.
        """
        # Percolator-specific mapping from internal param name -> CLI flag.
        # Extend this as new percolator params are supported.
        cli_flag_map = {
            "infer_proteins": "--picked-protein",
            "percolator_post_processing": None,  # positional, handled below
        }

        params_dict = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]

        # Params consumed elsewhere (create_input_file) or not CLI-relevant.
        skip_keys = {
            "bigger_scores_better",
            "validation_score_field",
            "delimiter",
            "enzyme",
            "database",
            "cpus",
        }

        for key, value in params_dict.items():
            if key in skip_keys:
                continue

            if key == "infer_proteins":
                if value is True:
                    utrace.urun_dict.command_list.append(cli_flag_map["infer_proteins"])
                    utrace.urun_dict.command_list.append(utrace.input_files[1].path)

                    target_proteins = (
                        utrace.output_files[0].path.parent / "target_protein_qvals.tsv"
                    )
                    decoy_proteins = (
                        utrace.output_files[0].path.parent / "decoy_protein_qvals.tsv"
                    )

                    utrace.urun_dict.command_list.append("-l")
                    utrace.urun_dict.command_list.append(f"{target_proteins}")
                    utrace.urun_dict.command_list.append("-L")
                    utrace.urun_dict.command_list.append(f"{decoy_proteins}")
                continue

            if key == "percolator_post_processing":
                if value is not None:
                    utrace.urun_dict.command_list.append(value)
                continue

            if value is True:
                flag = cli_flag_map.get(key, f"--{key}")
                utrace.urun_dict.command_list.append(flag)
            elif value is False or value is None:
                continue
            else:
                flag = cli_flag_map.get(key, f"--{key}")
                utrace.urun_dict.command_list.append(flag)
                utrace.urun_dict.command_list.append(value)

        return utrace

    def _load_input_dataframe(
        self,
        utrace: urgap.UTrace,
    ) -> tuple:
        """Load PSM csv files and drop rows with no target/decoy status.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            Tuple of (dataframe, original column list, node parameters, delimiter).
        """
        params_dict = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]
        delimiter = params_dict["delimiter"]

        unified_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV,
        )
        dfs = [pd.read_csv(f) for f in unified_files]
        df = pd.concat(dfs)
        old_columns = df.columns

        # Drop PSMs that couldn't be mapped to any protein (no target/decoy status
        # possible, expected for de novo callers like Instanovo where not every
        # sequence has a protein match)
        df = df.dropna(subset=["is_decoy"])
        df = df.sort_values(["spectrum_id", "rank"])

        df.loc[df["is_decoy"], "Label"] = "-1"
        df.loc[~df["is_decoy"], "Label"] = "1"

        return df, old_columns, params_dict, delimiter

    def _add_charge_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """One-hot encode charge states and derive charge-dependent features.

        Args:
            df: Combined PSM dataframe.

        Returns:
            Dataframe with charge one-hot columns and peptide length/mass added.
        """
        df = df.merge(
            pd.get_dummies(df.charge, prefix="charge").astype(int),
            left_index=True,
            right_index=True,
        )
        df = df.loc[1:]

        empty_charges = [
            f"charge_{i}" for i in range(11) if i not in df["charge"].unique()
        ]
        df.loc[:, empty_charges] = 0

        df["peplen"] = df["sequence"].str.len()
        df["mass"] = (df["exp_mz"] * df["charge"]) - (df["charge"] - 1) * PROTON
        df["dm"] = df["ucalc_mz"] - df["exp_mz"]
        df["absdm"] = abs(df["dm"])
        return df

    def _add_mass_and_delta_features(
        self,
        df: pd.DataFrame,
        params_dict: dict,
    ) -> pd.DataFrame:
        """Compute score, rank, delta-score, and enzyme-derived features.

        Args:
            df: PSM dataframe with charge features already added.
            params_dict: Node parameters for this percolator run.

        Returns:
            Dataframe with score, rank, delta score, and enzyme features added.
        """
        if len(df["search_engine"].unique()) != 1:
            msg = "Multiple engines detected in dataframe. Percolator can only handle one search engine at a time."
            raise PercolatorEngineMismatchError(msg)
        se = df["search_engine"].iloc[0]

        bigger_scores_better = params_dict["bigger_scores_better"][se]
        validate_score_field = params_dict["validation_score_field"][se]

        df["score"] = df[validate_score_field]
        if bigger_scores_better is False:
            df["score"] = -np.log10(df["score"])

        df["sp"] = df["score"].rank(method="max")
        df["lnrsp"] = np.log(df["sp"])

        threads = params_dict.get("cpus", 1)
        df = _apply_parallel(
            df.groupby("spectrum_id"), self.delta_score, threads=threads,
        ).reset_index(drop=True)

        df["enzn"] = df["enzn"].astype(int)
        df["enzc"] = df["enzc"].astype(int)
        df["enzint"] = df["missed_cleavages"].astype(int)
        return df

    def _add_peptide_and_protein_columns(
        self,
        df: pd.DataFrame,
        delimiter: str,
    ) -> pd.DataFrame:
        """Build the percolator-style Peptide, Proteins, and ScanNr columns.

        Args:
            df: PSM dataframe with score and delta features already added.
            delimiter: Delimiter used to split flanking-residue strings.

        Returns:
            Dataframe with Peptide, Proteins, and ScanNr columns populated.
        """
        df["modifications"] = df["modifications"].fillna("")
        df.loc[df["modifications"] == "", "Peptide"] = (
            df["sequence_pre_aa"].str.split(delimiter).str[0]
            + "." + df["sequence"] + "."
            + df["sequence_post_aa"].str.split(delimiter).str[0]
        )
        df.loc[df["modifications"] != "", "Peptide"] = (
            df["sequence_pre_aa"].str.split(delimiter).str[0]
            + "." + df["sequence"] + "[#" + df["modifications"] + "]."
            + df["sequence_post_aa"].str.split(delimiter).str[0]
        )

        df["Proteins"] = df["protein_id"]
        df["ScanNr"] = df["spectrum_id"]
        return df

    def _write_feature_file(
        self,
        df: pd.DataFrame,
        old_columns: object,
        utrace: urgap.UTrace,
    ) -> os.PathLike:
        """Write the percolator input tsv and the target/decoy merge frame.

        Args:
            df: Fully-featured PSM dataframe.
            old_columns: Column names from the original input dataframe(s).
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            Path to the generated percolator input tsv file.
        """
        features = [
            "PSMId", "Label", "ScanNr", "lnrsp", "deltlcn", "deltcn",
            "score", "sp", "mass", "peplen",
            "charge_1", "charge_2", "charge_3", "charge_4", "charge_5",
            "charge_6", "charge_7", "charge_8", "charge_9", "charge_10",
            "enzn", "enzc", "enzint", "dm", "absdm", "Peptide", "Proteins",
        ]

        df = df.reset_index()
        df = df.rename(columns={"index": "PSMId"})

        feature_df = df[features]
        fname = utrace.output_files[0].path.parent / "percolator_input.tsv"
        self.tmp_files.append(fname)
        feature_df = feature_df.sort_values("ScanNr")
        feature_df["Proteins"] = (
            feature_df["Proteins"].str.split(r"<\|>").str.join("\t")
        )

        feature_df.to_csv(fname, sep="\t", index=False)
        self.remove_quotes(fname)

        new_columns = [*list(old_columns), "PSMId"]
        self.merged_frame = utrace.output_files[0].path.parent / "merge_frame.csv"
        self.tmp_files.append(self.merged_frame)
        df[new_columns].reset_index().to_csv(self.merged_frame, index=False)
        return fname

    def create_input_file(
        self,
        utrace: urgap.UTrace,
    ) -> os.PathLike:
        """Create the input file following percolator convention.

        Args:
            utrace: Combination of urun_dict, ufile_list and unode.meta.

        Returns:
            Path to input file.
        """
        df, old_columns, params_dict, delimiter = self._load_input_dataframe(utrace)
        df = self._add_charge_features(df)
        df = self._add_mass_and_delta_features(df, params_dict)
        df = self._add_peptide_and_protein_columns(df, delimiter)
        return self._write_feature_file(df, old_columns, utrace)

    def remove_quotes(self, file: os.PathLike) -> None:
        """Remove quotes from each line within the input file.

        Args:
            file: Path to the file to strip quotes from, in place.
        """
        no_quotes = file.parent / "no_quotes.txt"
        with file.open() as fin, no_quotes.open("w") as fout:
            for line in fin:
                cleaned_line = line.replace('"', "")
                fout.write(cleaned_line)
        shutil.move(no_quotes, file)

    def delta_score(self, grp: pd.DataFrame) -> pd.DataFrame:
        """Calculate the delta score.

        Args:
            grp: Dataframe of PSMs for a single spectrum, grouped by spectrum_id.

        Returns:
            Input dataframe + delta score columns.
        """
        grp = grp.sort_values("rank")
        grp["deltcn"] = (grp["score"] - grp.shift(-1)["score"]).fillna(0) / grp["score"]
        grp["deltlcn"] = (grp["score"] - min(grp["score"])) / grp["score"]
        grp["deltcn"] = grp["deltcn"].replace([np.inf, -np.inf], 0.0)
        grp["deltlcn"] = grp["deltlcn"].replace([np.inf, -np.inf], 0.0)
        return grp
