"""Generate Target Decoy Fasta resource."""
import argparse
import itertools
import logging
import random

from collections import defaultdict
from pathlib import Path

import ahocorasick
import regex as re
import tqdm


def digest(pattern: str, peptide: str) -> list[str]:
    """Digest a peptide using an enzyme regex expression.

    Args:
        pattern (str): regex pattern indicating enzyme specific cleavage site
        peptide (str): peptide sequence

    Returns:
        (list): cleaved peptide fragments
    """
    return list(filter(None, re.split(pattern, peptide)))


def convert_fasta_to_dict(fasta_path: Path) -> dict[str, str]:
    """Convert fasta file to python dictionary.

    Args:
        fasta_path (Path): path object of fasta file

    Returns:
        (dict): dict where keys are protein identifiers and values are sequences
    """
    fasta_dict = {}
    sequence = None
    protein_id = None
    with Path(fasta_path).open() as fasta:
        for line in fasta:
            if line.startswith(">"):
                if sequence is not None and protein_id is not None:
                    fasta_dict[protein_id] = sequence
                protein_id = line[1:-1]
                sequence = ""
            else:
                sequence += line[:-1]
        if protein_id is not None:
            fasta_dict[protein_id] = sequence
    return fasta_dict


def _merge_fasta_id(fasta_id_list: list[str]) -> str:
    n_ids = len(fasta_id_list)
    if n_ids == 1:
        merged_id = fasta_id_list[0]
    else:
        merged_id = f"{fasta_id_list[0]}_NumberOfIdenticalSequences_{n_ids}"
    return merged_id


def _generate_shuffled_candidate(peptide: str) -> str:
    """Generate a single shuffled candidate sequence, preserving N/C termini.

    Args:
        peptide (str): peptide sequence to shuffle

    Returns:
        (str): candidate shuffled sequence
    """
    preserve_n = peptide[0]
    preserve_c = peptide[-1]

    if len(peptide) < 8:
        shuffled_mid_sequence = "".join(
            random.choice(
                list(itertools.permutations(peptide[1:-1], len(peptide) - 2)),
            ),
        )
    else:
        mid_seq = list(peptide[1:-1])
        mid_seq = sorted(mid_seq, key=lambda _: random.random())
        shuffled_mid_sequence = "".join(mid_seq)

    return preserve_n + shuffled_mid_sequence + preserve_c


def _has_long_target_overlap(auto: ahocorasick.Automaton, candidate: str) -> bool:
    """Check whether a candidate sequence overlaps too long with any target peptide.

    Args:
        auto (ahocorasick.Automaton): automaton built from target peptides
        candidate (str): candidate shuffled sequence

    Returns:
        (bool): True if an overlap longer than 5 residues is found
    """
    return any(len(match) > 5 for _ind, match in auto.iter_long(candidate))


def _shuffle_single_peptide(
    peptide: str,
    target_peptides: set,
    enzyme_pattern: str,
    auto: ahocorasick.Automaton,
) -> str:
    """Attempt to generate a shuffled decoy for a single peptide.

    Args:
        peptide (str): peptide sequence to shuffle
        target_peptides (set): set of all target peptides
        enzyme_pattern (str): regex pattern indicating enzyme specific cleavage site
        auto (ahocorasick.Automaton): automaton built from target peptides

    Returns:
        (str): shuffled sequence, the original peptide (if too short to shuffle),
            or "<IMMUTABLE>" if no valid shuffle could be found
    """
    if len(peptide) < 5:
        return peptide

    internal_cleavages = re.findall(enzyme_pattern, peptide)

    for _ in range(10000):
        candidate = _generate_shuffled_candidate(peptide)

        if candidate in target_peptides:
            continue

        if _has_long_target_overlap(auto, candidate):
            continue

        new_internal_cleavages = re.findall(enzyme_pattern, candidate)
        if len(internal_cleavages) == len(new_internal_cleavages):
            return candidate

    return "<IMMUTABLE>"


def shuffle_peptides(target_peptides: set, enzyme_pattern: str) -> tuple[list[str], list[str]]:
    """Shuffle set of target_peptides to generate decoy peptides.

    Args:
        target_peptides (set): set of all to-be-shuffled peptides
        enzyme_pattern (str): regex pattern indicating enzyme specific cleavage site

    Returns:
        decoys (list): generated decoy peptides
        immutable_peptides (list): list of peptides for which no decoy could be generated
    """
    decoys = []
    immutable_peptides = []
    auto = ahocorasick.Automaton()
    for seq in target_peptides:
        auto.add_word(seq, seq)
    auto.make_automaton()

    for peptide in tqdm.tqdm(target_peptides):
        shuffled_seq = _shuffle_single_peptide(
            peptide, target_peptides, enzyme_pattern, auto,
        )
        if shuffled_seq == "<IMMUTABLE>":
            immutable_peptides.append(peptide)
            decoys.append(peptide)
        else:
            decoys.append(shuffled_seq)

    return decoys, immutable_peptides


def main(
    input_files: list[Path],
    mode: str,
    enzyme_pattern: str,
    decoy_tag: str,
    output_file: Path,
    immutable_peptides_file: Path | None = None,
    seed: str | None = None,
) -> None:
    """Create a target decoy database from one (or more) fasta files.

    Args:
        input_files (list): list of path object to fasta files
        mode (str): mode of decoy generation ("reverse_protein" or "shuffle_peptide")
        enzyme_pattern (str): regex pattern indicating enzyme specific cleavage site
        decoy_tag (str): prefix protein identifiers that are decoys are tagged with
        output_file (Path): output file path
        immutable_peptides_file (Path, optional): output file path for immutable peptides
        seed (str): seed for random module
    """
    if seed is not None:
        logging.info("Target decoy generation randomization is using seed: %s", seed)
        random.seed(seed)

    sequence_dicts = []
    for fasta in input_files:
        fasta_dict = convert_fasta_to_dict(fasta)
        sequence_dicts.append(
            {seq: protein_id for protein_id, seq in fasta_dict.items() if "REVERSED" not in protein_id},
        )

    sequence_dict = defaultdict(list)
    for d in sequence_dicts:
        for seq, protein_id in d.items():
            sequence_dict[seq].append(protein_id)

    sequence_dict = {seq: _merge_fasta_id(ids) for seq, ids in sequence_dict.items()}
    target_decoy_database = {}
    if mode == "reverse_protein":
        target_decoy_database = {
            protein_id: {"target_sequence": seq, "decoy_sequence": seq[::-1]}
            for seq, protein_id in sequence_dict.items()
        }
        immutable_peptides = []

    elif mode == "shuffle_peptide":
        digested_seqs = [digest(enzyme_pattern, seq) for seq in sequence_dict]
        unique_fragments = set(itertools.chain.from_iterable(digested_seqs))
        unique_fragments_shuffled, immutable_peptides = shuffle_peptides(
            target_peptides=unique_fragments, enzyme_pattern=enzyme_pattern,
        )
        unique_fragments = dict(zip(unique_fragments, unique_fragments_shuffled, strict=False))
        shuffled_protein_seqs = [
            "".join(unique_fragments[fragment] for fragment in protein)
            for protein in digested_seqs
        ]
        target_decoy_database = {
            protein_id: {
                "target_sequence": seq,
                "decoy_sequence": shuffle_seq,
            }
            for (seq, protein_id), shuffle_seq in zip(
                sequence_dict.items(), shuffled_protein_seqs, strict=False,
            )
        }

    else:
        logging.error("Unknown target decoy generation mode")
        raise ValueError

    logging.info(
        "Generated decoys for %d of %d database proteins",
        len(target_decoy_database),
        len(sequence_dict),
    )
    logging.info("%d immutable peptides recorded.", len(immutable_peptides))

    with Path(output_file).open("w") as f:
        logging.info("Writing target decoy database to file...")
        for protein_id, seqs in target_decoy_database.items():
            # Write target
            f.write(f">{protein_id}\n")
            f.write(
                re.sub(
                    pattern="(.{80})",
                    repl="\\1\n",
                    string=seqs["target_sequence"],
                    flags=re.DOTALL,
                )
                + "\n",
            )
            # Write decoy
            if seqs["decoy_sequence"] != "<IMMUTABLE>":
                f.write(f">{decoy_tag}{protein_id}\n")
                f.write(
                    re.sub(
                        pattern="(.{80})",
                        repl="\\1\n",
                        string=seqs["decoy_sequence"],
                        flags=re.DOTALL,
                    )
                    + "\n",
                )

    if immutable_peptides_file is not None:
        with Path(immutable_peptides_file).open("w") as f:
            f.writelines(pep + "\n" for pep in immutable_peptides)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", dest="input_files", nargs="+")
    parser.add_argument(
        "--output",
        dest="output_file",
    )
    parser.add_argument(
        "--immutable_file",
        required=False,
        dest="immutable_file",
    )
    parser.add_argument(
        "--enzyme_pattern",
        dest="enz_pat",
    )
    parser.add_argument(
        "--decoy_tag",
        dest="decoy_tag",
    )
    parser.add_argument(
        "--seed",
        required=False,
        dest="seed",
    )
    parser.add_argument(
        "--mode",
        dest="mode",
    )
    known_args = parser.parse_args()
    main(
        input_files=known_args.input_files,
        mode=known_args.mode,
        enzyme_pattern=known_args.enz_pat,
        decoy_tag=known_args.decoy_tag,
        output_file=known_args.output_file,
        immutable_peptides_file=known_args.immutable_file,
        seed=known_args.seed,
    )
