"""Build training-ready antigenicity inputs from the FluAntigen H1N1 pairs.

Reproduces the official pipeline stages (data/antigenicity/prepare.py output
format + data/antigenicity/split.py random split), using accession IDs as the
join key throughout (both the pairs table and HIN1_MERGED already key on HA
EPI accessions, so no name-alias step is needed here):

  1. Export the virus-side and reference-side HA sequences referenced by the
     matched pairs into separate FASTA files.
  2. Run MMseqs2 (query=virus, target=reference) to get pairwise alignments,
     exactly like data/antigenicity/build_test_pairs.py / pipeline/run.sh.
  3. Join alignments onto the HI pairs to produce the aligned
     virus_seq/reference_seq columns required by --data_module
     hi_regression_aln (see data/data_modules/protein_classification_dm.py).
  4. Random 80/10/10 train/valid/test split with a fixed seed, matching
     data/antigenicity/split.py's ratios (not the original paper's split,
     which was never published).
"""
import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from Bio import SeqIO

EPI_RE = re.compile(r"^EPI\d+$")


def ha_id(description):
    values = [item.strip() for item in description.split("|")]
    matches = [item for item in values if EPI_RE.fullmatch(item)]
    return matches[-1] if matches else None


def load_sequences(merged_dir):
    accession_to_sequence = {}
    for record in SeqIO.parse(str(merged_dir / "ha.fasta"), "fasta"):
        accession = ha_id(record.description)
        if accession is not None:
            accession_to_sequence[accession] = str(record.seq)
    return accession_to_sequence


def write_fasta(path, accession_to_sequence, accessions):
    with path.open("w") as handle:
        for accession in sorted(accessions):
            handle.write(">%s\n%s\n" % (accession, accession_to_sequence[accession]))


def read_m8(path):
    pair_to_alignment = {}
    with path.open() as handle:
        for line in handle:
            query, target, qaln, taln, qstart, qend, tstart, tend, mismatch = line.strip().split()
            pair_to_alignment.setdefault((query, target), (qaln, taln))
    return pair_to_alignment


def save_split(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["virus", "reference", "virus_seq", "reference_seq", "hi"])
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-dataset", type=Path, required=True)
    parser.add_argument("--pairs-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1005)
    args = parser.parse_args()

    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    accession_to_sequence = load_sequences(args.merged_dataset)

    with args.pairs_csv.open(newline="") as handle:
        pairs = list(csv.DictReader(handle))
    hi_folds_path = output / "a_h1n1_hi_folds.csv"
    with hi_folds_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["virus", "reference", "hi"])
        writer.writeheader()
        for row in pairs:
            writer.writerow({"virus": row["virus"], "reference": row["reference"], "hi": row["hi"]})

    virus_accessions = {row["virus"] for row in pairs}
    reference_accessions = {row["reference"] for row in pairs}
    virus_fasta = output / "a_h1n1_virus.fasta"
    reference_fasta = output / "a_h1n1_vaccine.fasta"
    write_fasta(virus_fasta, accession_to_sequence, virus_accessions)
    write_fasta(reference_fasta, accession_to_sequence, reference_accessions)

    m8_path = output / "a_h1n1.m8"
    subprocess.run(
        ["mmseqs", "easy-search", str(virus_fasta), str(reference_fasta), str(m8_path),
         str(output / "mmseqs_tmp"), "--format-output",
         "query,target,qaln,taln,qstart,qend,tstart,tend,mismatch", "--max-seqs", "2000"],
        check=True,
    )
    pair_to_alignment = read_m8(m8_path)

    aligned_rows = []
    missing_alignment = 0
    for row in pairs:
        key = (row["virus"], row["reference"])
        alignment = pair_to_alignment.get(key)
        if alignment is None:
            missing_alignment += 1
            continue
        qaln, taln = alignment
        aligned_rows.append({
            "virus": row["virus"], "reference": row["reference"],
            "virus_seq": qaln, "reference_seq": taln, "hi": row["hi"],
        })

    rng = np.random.RandomState(args.seed)
    index = rng.permutation(len(aligned_rows))
    n_train = int(len(index) * 0.8)
    n_valid = int(len(index) * 0.9)
    train_index, valid_index, test_index = index[:n_train], index[n_train:n_valid], index[n_valid:]

    split_dir = output / "random_split"
    split_dir.mkdir(exist_ok=True)
    save_split(split_dir / "train.csv", [aligned_rows[i] for i in train_index])
    save_split(split_dir / "valid.csv", [aligned_rows[i] for i in valid_index])
    save_split(split_dir / "test.csv", [aligned_rows[i] for i in test_index])

    report = {
        "input_pairs": len(pairs),
        "unique_virus_accessions": len(virus_accessions),
        "unique_reference_accessions": len(reference_accessions),
        "pairs_with_mmseqs_alignment": len(aligned_rows),
        "pairs_without_alignment_dropped": missing_alignment,
        "seed": args.seed,
        "split_sizes": {"train": len(train_index), "valid": len(valid_index), "test": len(test_index)},
        "note": "Random 80/10/10 split (data/antigenicity/split.py ratios); "
                "this is NOT the original paper's train/valid/test split, which was never published.",
    }
    (output / "CHECK_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
