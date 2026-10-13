"""Reproduce the official VaxSeer prepare.py logic on FluAntigen H1N1 data.

Mirrors data/antigenicity/prepare.py exactly:
  1. Build accession/name/sequence indices from the HA FASTA (min length 553).
  2. Resolve virus/reference names with the official two-digit-year and
     parenthetical alias rules (no additional abbreviation dictionary).
  3. Convert to log2 titer, then to relative-to-self-titer HI (clipped at 0),
     using name-level self-pairs exactly like the official script.
  4. Expand every matched name to the FULL cross product of its accessions
     (no manual disambiguation), then aggregate identical (virus_seq, ref_seq)
     pairs by mean, exactly like prepare.py.
"""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

MIN_SEQ_LENGTH = 553


def normalize(name):
    return str(name).replace(" ", "_").replace("-", "_").lower()


def read_fasta(path, max_date="9999-12-31"):
    max_year, max_month = int(max_date.split("-")[0]), int(max_date.split("-")[1])
    accid2seq = {}
    id2records = defaultdict(set)
    id2accids = defaultdict(set)
    seq2accid = defaultdict(list)
    for record in SeqIO.parse(str(path), "fasta"):
        descs = record.description.split("|")
        try:
            if len(str(record.seq)) < MIN_SEQ_LENGTH:
                continue
            ha_id = None
            date = None
            for x in descs:
                if "EPI" in x and "EPI_ISL" not in x and re.fullmatch(r"EPI\d+", x):
                    ha_id = x
                if re.match(r"[0-9]{4}-[0-9]{2}(-[0-9]{2})?$", x):
                    date = x
            if ha_id is None or date is None:
                continue
            year, month = int(date.split("-")[0]), int(date.split("-")[1])
            if year > max_year or (year == max_year and month >= max_month):
                continue
            strain_name = normalize(descs[0])
            id2records[strain_name].add(str(record.seq))
            id2accids[strain_name].add(ha_id)
            accid2seq[ha_id] = str(record.seq)
            seq2accid[str(record.seq)].append(ha_id)
        except Exception:
            continue
    return accid2seq, id2records, id2accids, seq2accid


def search_id(names, query):
    return [name for name in names if query in name]


def alias_name(all_names, query_name):
    if len(query_name.split("/")[-1]) == 2:
        query_name = "/".join(query_name.split("/")[:-1]) + "/20" + query_name.split("/")[-1]
    if "(" in query_name:
        groups = re.findall(r"[(](.*?)[)]", query_name)
        if not groups:
            return query_name
        short = groups[0]
        if len(short.split("/")[-1]) == 2:
            if int(short.split("/")[-1][0]) <= 2:
                short = "/".join(short.split("/")[:-1]) + "/20" + short.split("/")[-1]
            else:
                short = "/".join(short.split("/")[:-1]) + "/19" + short.split("/")[-1]
        prefix = query_name.split("(")[0]
        matches = search_id(all_names, short)
        new_name = None
        for match in matches:
            if prefix in match:
                new_name = match
        query_name = new_name if new_name is not None else short
    return query_name


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fluantigen-xlsx", type=Path, required=True)
    parser.add_argument("--sheet", default="H1N1")
    parser.add_argument("--merged-dataset", type=Path, required=True)
    parser.add_argument("--max-date-exclude", default="9999-12-31")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    accid2seq, id2records, id2accids, seq2accid = read_fasta(
        args.merged_dataset / "ha.fasta", max_date=args.max_date_exclude)
    all_names = set(id2records.keys())

    frame = pd.read_excel(args.fluantigen_xlsx, sheet_name=args.sheet)
    frame = frame.rename(columns={"VirusStrain": "virus", "SerumStrain": "reference", "HI titre": "hi"})
    frame = frame[["virus", "reference", "hi"]].dropna(subset=["virus", "reference"])
    frame["hi"] = frame["hi"].apply(to_float)
    total_rows = len(frame)
    censored_or_missing_hi = int(frame["hi"].isna().sum())
    frame = frame.dropna(subset=["hi"])
    frame = frame[frame["hi"] > 0]

    frame["virus"] = frame["virus"].apply(normalize)
    frame["reference"] = frame["reference"].apply(normalize)

    virus_set = list(set(frame["virus"]))
    reference_set = list(set(frame["reference"]))

    virus_name_to_accids = {}
    for name in virus_set:
        resolved = alias_name(all_names, name)
        if resolved in id2accids:
            virus_name_to_accids[name] = id2accids[resolved]
    reference_name_to_accids = {}
    for name in reference_set:
        resolved = alias_name(all_names, name)
        if resolved in id2accids:
            reference_name_to_accids[name] = id2accids[resolved]

    frame["hi"] = np.log2(frame["hi"])
    ref2ref_hi = {}
    for ref, virus, hi in zip(frame["reference"], frame["virus"], frame["hi"]):
        if virus == ref:
            ref2ref_hi[ref] = hi

    normalized_rows = []
    for ref, virus, hi in zip(frame["reference"], frame["virus"], frame["hi"]):
        if ref in ref2ref_hi:
            normalized_rows.append((virus, ref, max(ref2ref_hi[ref] - hi, 0.0)))
    dropped_no_self_titer = len(frame) - len(normalized_rows)

    new_pairs = defaultdict(list)
    unmatched_virus_names = defaultdict(int)
    unmatched_reference_names = defaultdict(int)
    matched_name_rows = 0
    for virus, ref, hi in normalized_rows:
        virus_ok = virus in virus_name_to_accids
        ref_ok = ref in reference_name_to_accids
        if not virus_ok:
            unmatched_virus_names[virus] += 1
        if not ref_ok:
            unmatched_reference_names[ref] += 1
        if not (virus_ok and ref_ok):
            continue
        matched_name_rows += 1
        for virus_accid in virus_name_to_accids[virus]:
            for ref_accid in reference_name_to_accids[ref]:
                virus_seq = accid2seq[virus_accid]
                ref_seq = accid2seq[ref_accid]
                new_pairs[(virus_seq, ref_seq)].append(hi)

    output_rows = []
    for (virus_seq, ref_seq), values in new_pairs.items():
        output_rows.append({
            "virus": seq2accid[virus_seq][0],
            "reference": seq2accid[ref_seq][0],
            "hi": np.mean(values),
            "n_source_records": len(values),
        })
    output_frame = pd.DataFrame(output_rows)
    output_frame.to_csv(output / "a_h1n1_pairs_fluantigen_official_style.csv", index=False)

    with (output / "unmatched_virus_names.tsv").open("w") as handle:
        handle.write("name\tcount\n")
        for name, count in sorted(unmatched_virus_names.items(), key=lambda x: -x[1]):
            handle.write("%s\t%d\n" % (name, count))
    with (output / "unmatched_reference_names.tsv").open("w") as handle:
        handle.write("name\tcount\n")
        for name, count in sorted(unmatched_reference_names.items(), key=lambda x: -x[1]):
            handle.write("%s\t%d\n" % (name, count))

    report = {
        "fluantigen_total_rows": total_rows,
        "censored_or_non_numeric_hi_dropped": censored_or_missing_hi,
        "rows_after_hi_cleaning": len(frame),
        "rows_dropped_no_self_titer_reference": dropped_no_self_titer,
        "rows_with_self_titer": len(normalized_rows),
        "rows_with_both_names_matched": matched_name_rows,
        "unique_unmatched_virus_names": len(unmatched_virus_names),
        "unique_unmatched_reference_names": len(unmatched_reference_names),
        "final_unique_sequence_pairs": len(output_rows),
        "note": "Full cross-product expansion per matched name, then mean over identical (virus_seq, reference_seq) pairs, exactly as data/antigenicity/prepare.py.",
    }
    (output / "CHECK_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
