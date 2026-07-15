"""urgap percolator_3_7_1 wrapper."""

import os
import shutil
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
from chemical_composition.chemical_composition_kb import PROTON

import urgap


class percolator_3_7_1(urgap.unode.UNodeBase):
    """
    urgap wrapper for the percolator_3_7_1 executable.

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
        The, M., MacCoss, M. J., Noble, W. S., & Käll, L. (2016). Fast and Accurate Protein False Discovery Rates on Large-Scale Proteomics Data Sets with Percolator 3.0. 
        In Journal of the American Society for Mass Spectrometry (Vol. 27, Issue 11, pp. 1719–1727). American Chemical Society (ACS). https://doi.org/10.1007/s13361-016-1460-7
        """,
    }

    def __init__(self, *args, **kwargs):
        """Initialize percolator_3_7_1 class."""
        super(percolator_3_7_1, self).__init__(*args, **kwargs)
        pass

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
            urgap.uftypes.proteomics.validator.PERCOLATOR_CSV
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
        utrace = self.create_command_list(utrace)
        return utrace

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
            utrace.output_files[0].path.name + "_percolator_out_fixed.tsv"
        )

        fixed_decoys_path = utrace.output_files[0].path.parent / Path(
            utrace.output_files[0].path.name + "_percolator_out_fixed_decoys.tsv"
        )

        self.tmp_files.append(self.result_psms)
        self.tmp_files.append(self.decoy_psms)
        self.tmp_files.append(fixed_targets_path)
        self.tmp_files.append(fixed_decoys_path)

        for broken_path, fixed_path in [
            (self.decoy_psms, fixed_decoys_path),
            (self.result_psms, fixed_targets_path),
        ]:
            with open(str(broken_path)) as fin, open(str(fixed_path), "wt") as fout:
                for i, line in enumerate(fin):
                    if i == 0:
                        fout.write(line)
                        continue
                    l = line.split("\t")[:5]
                    prot = " ".join(line.split("\t")[5:])
                    l.append(prot)
                    fout.write("\t".join(l))
        # rename files again

        output_decoys = pd.read_csv(fixed_decoys_path, sep="\t", index_col=False)
        output_decoys = output_decoys[["PSMId", "q-value", "posterior_error_prob"]]

        output_targets = pd.read_csv(fixed_targets_path, sep="\t", index_col=False)
        output_targets = output_targets[["PSMId", "q-value", "posterior_error_prob"]]

        qvals = pd.concat([output_targets, output_decoys])

        unified_df = pd.read_csv(self.merged_frame)

        final_df = pd.merge(
            unified_df, qvals, left_on="PSMId", right_on="PSMId", how="left"
        )
        final_df = final_df[~final_df["q-value"].isna()]
        idx = self.output_type_dict[".percolator.csv"][0]
        final_df.to_csv(utrace.output_files[idx].path)

        # Part specific for only version 3.7.1
        if self.META_INFO["unode_version"] == "3.7.1":
            if (
                utrace.output_files[0].path.parent / "target_protein_qvals.tsv"
            ).exists():
                utrace.extend_output_files_by_uftype(
                    urgap.uftypes.proteomics.validator.PERCOLATOR_CSV
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
        CLI_FLAG_MAP = {
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
                    utrace.urun_dict.command_list.append(CLI_FLAG_MAP["infer_proteins"])
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
                flag = CLI_FLAG_MAP.get(key, f"--{key}")
                utrace.urun_dict.command_list.append(flag)
            elif value is False or value is None:
                continue
            else:
                flag = CLI_FLAG_MAP.get(key, f"--{key}")
                utrace.urun_dict.command_list.append(flag)
                utrace.urun_dict.command_list.append(value)

        return utrace

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
        req_headers = ["PSMId", "Label", "ScanNr", "Peptide", "Proteins"]
        features = [
            "PSMId", "Label", "ScanNr", "lnrsp", "deltlcn", "deltcn",
            "score", "sp", "mass", "peplen",
            "charge_1", "charge_2", "charge_3", "charge_4", "charge_5",
            "charge_6", "charge_7", "charge_8", "charge_9", "charge_10",
            "enzn", "enzc", "enzint", "dm", "absdm", "Peptide", "Proteins",
        ]

        all_headers = req_headers + features
        default_directions_features = {col: 0 for col in features}
        default_directions_features.update({col: "-" for col in req_headers})
        default_directions_features["PSMId"] = "DefaultDirection"

        params_dict = utrace.urun_dict.parameters[
            f"{self.META_INFO['unode_full_identifier']}"
        ]

        delimiter = params_dict["delimiter"]

        unified_files = utrace.input_files.get_path_objects_by_uftype(
            urgap.uftypes.proteomics.converter.PYIOHAT_CSV
        )
        dfs = []
        for f in unified_files:
            _df = pd.read_csv(f)
            dfs.append(_df)
        df = pd.concat(dfs)

        old_columns = df.columns

        # Drop PSMs that couldn't be mapped to any protein (no target/decoy status possible,
        # expected for de novo callers like Instanovo where not every sequence has a protein match)
        df = df.dropna(subset=["is_decoy"])

        df = df.sort_values(["spectrum_id", "rank"])

        df.loc[df["is_decoy"] == True, "Label"] = "-1"
        df.loc[df["is_decoy"] == False, "Label"] = "1"

        # One hot encode charges
        df = pd.merge(
            df,
            pd.get_dummies(df.charge, prefix="charge").astype(int),
            left_index=True,
            right_index=True,
        )
        df = df.loc[1:]

        empty_charges = [
            f"charge_{i}" for i in range(0, 11) if i not in df["charge"].unique()
        ]
        df.loc[:, empty_charges] = 0

        df["peplen"] = df["sequence"].str.len()
        df["mass"] = (df["exp_mz"] * df["charge"]) - (df["charge"] - 1) * PROTON
        df["dm"] = df["ucalc_mz"] - df["exp_mz"]
        df["absdm"] = abs(df["dm"])

        if len(df["search_engine"].unique()) == 1:
            se = df["search_engine"].iloc[0]
        else:
            print(df["search_engine"].unique())
            raise Exception(
                "Multiple engines detected in dataframe. Percolator can only handle one search engine at a time."
            )

        bigger_scores_better = params_dict["bigger_scores_better"][se]
        validate_score_field = params_dict["validation_score_field"][se]

        df["score"] = df[validate_score_field]
        if bigger_scores_better is False:
            df["score"] = -np.log10(df["score"])

        df["sp"] = df["score"].rank(method="max")
        df["lnrsp"] = np.log(df["sp"])

        def applyParallel(dfGrouped, func, threads=-1):
            if threads == -1:
                threads = cpu_count()
            with Pool(threads) as p:
                ret_list = p.map(func, [group for name, group in dfGrouped])
            return pd.concat(ret_list)

        threads = params_dict.get("cpus", 1)
        df = applyParallel(
            df.groupby("spectrum_id"), self.delta_score, threads=threads
        ).reset_index(drop=True)

        df["enzn"] = df["enzn"].astype(int)
        df["enzc"] = df["enzc"].astype(int)
        df["enzint"] = df["missed_cleavages"].astype(int)

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

        # TEMPORARY DEBUG: copy the file somewhere permanent before percolator runs
        import shutil as _shutil
        _shutil.copy(fname, "/shared/rc/proteome/urgap/connor_example_scripts/debug_percolator_input.tsv")

        _new = list(old_columns) + ["PSMId"]
        self.merged_frame = utrace.output_files[0].path.parent / "merge_frame.csv"
        self.tmp_files.append(self.merged_frame)
        df[_new].reset_index().to_csv(self.merged_frame, index=False)
        return fname

    def remove_quotes(self, file: os.PathLike):
        """Remove quotes from each line within the input file.

        Args:
            Path to file.
        """
        no_quotes = file.parent / "no_quotes.txt"
        with open(file) as fin, open(no_quotes, "wt") as fout:
            for line in fin:
                line = line.replace('"', "")
                fout.write(line)
        shutil.move(no_quotes, file)

    def delta_score(self, grp: pd.DataFrame) -> pd.DataFrame:
        """Calculate the delta score.

        Args:
            Input dataframe.

        Returns:
            Input dataframe + delta score columns.
        """
        grp = grp.sort_values("rank")
        grp["deltcn"] = (grp["score"] - grp.shift(-1)["score"]).fillna(0) / grp["score"]
        grp["deltlcn"] = (grp["score"] - min(grp["score"])) / grp["score"]
        grp["deltcn"] = grp["deltcn"].replace([np.inf, -np.inf], 0.0)
        grp["deltlcn"] = grp["deltlcn"].replace([np.inf, -np.inf], 0.0)
        return grp
