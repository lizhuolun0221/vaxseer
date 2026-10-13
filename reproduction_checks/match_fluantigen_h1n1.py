"""Match FluAntigen H1N1 HI records to HIN1_MERGED HA sequences by strain name.

Uses the same alias-resolution rules as the official prepare.py / repository
audits: exact normalized match, two-digit year expansion, and parenthetical
short-name lookup. No fuzzy or biological sequence matching is performed.
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


EPI_RE = re.compile(r"^EPI\d+$")


def normalize(name):
    return str(name).strip().replace(" ", "_").replace("-", "_").lower()


def resolve(query, available_names):
    if query in available_names:
        return [query], "exact"
    converted = query
    tail = converted.split("/")[-1]
    if len(tail) == 2 and tail.isdigit():
        converted = "/".join(converted.split("/")[:-1]) + "/20" + tail
    if converted in available_names and converted != query:
        return [converted], "short_year"
    if "(" in query:
        groups = re.findall(r"[(](.*?)[)]", query)
        if not groups:
            return [], "malformed_parenthesis"
        short = groups[0]
        short_tail = short.split("/")[-1]
        if len(short_tail) == 2 and short_tail.isdigit():
            century = "20" if int(short_tail[0]) <= 2 else "19"
            short = "/".join(short.split("/")[:-1]) + "/" + century + short_tail
        prefix = query.split("(")[0]
        candidates = sorted(n for n in available_names if short in n and prefix in n)
        if candidates:
            return candidates, "parenthetical_prefix"
        if short in available_names:
            return [short], "parenthetical_short"
        return [], "parenthetical_unmatched"
    return [], "unmatched"


def ha_id(description):
    values = [item.strip() for item in description.split("|")]
    matches = [item for item in values if EPI_RE.fullmatch(item)]
    return matches[-1] if matches else None


def load_merged(merged_dir):
    accession_to_sequence = {}
    for record in SeqIO.parse(str(merged_dir / "ha.fasta"), "fasta"):
        accession = ha_id(record.description)
        if accession is not None:
            accession_to_sequence[accession] = str(record.seq)
    name_to_accessions = defaultdict(set)
    with (merged_dir / "metadata.csv").open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("HA Segment_Id") or "").split("|")[0].strip()
            name = row.get("Isolate_Name")
            if accession and name and accession in accession_to_sequence:
                name_to_accessions[normalize(name)].add(accession)
    return accession_to_sequence, name_to_accessions


def match_name(raw_name, name_to_accessions, available_names):
    query = normalize(raw_name)
    candidates, rule = resolve(query, available_names)
    accessions = set()
    for candidate in candidates:
        accessions |= name_to_accessions.get(candidate, set())
    return accessions, rule


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fluantigen-xlsx", type=Path, required=True)
    parser.add_argument("--sheet", default="H1N1")
    parser.add_argument("--merged-dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    accession_to_sequence, name_to_accessions = load_merged(args.merged_dataset)
    available_names = set(name_to_accessions.keys())

    frame = pd.read_excel(args.fluantigen_xlsx, sheet_name=args.sheet)
    rows_out = []
    mismatch_rows = []
    ambiguous_rows = []
    match_rule_counts = defaultdict(int)
    unmatched_virus_names = defaultdict(int)
    unmatched_serum_names = defaultdict(int)

    for record in frame.to_dict(orient="records"):
        virus_name = record.get("VirusStrain")
        serum_name = record.get("SerumStrain")
        hi_value = record.get("HI titre")
        if pd.isna(virus_name) or pd.isna(serum_name):
            mismatch_rows.append({"reason": "missing_strain_name", "virus": virus_name,
                                  "reference": serum_name, "hi_titre": hi_value})
            continue
        virus_accessions, virus_rule = match_name(virus_name, name_to_accessions, available_names)
        serum_accessions, serum_rule = match_name(serum_name, name_to_accessions, available_names)
        match_rule_counts["virus:" + virus_rule] += 1
        match_rule_counts["reference:" + serum_rule] += 1

        if not virus_accessions:
            unmatched_virus_names[str(virus_name)] += 1
        if not serum_accessions:
            unmatched_serum_names[str(serum_name)] += 1

        if not virus_accessions or not serum_accessions:
            mismatch_rows.append({"reason": "name_not_found", "virus": virus_name,
                                  "reference": serum_name, "hi_titre": hi_value,
                                  "virus_rule": virus_rule, "reference_rule": serum_rule})
            continue

        virus_sequences = {accession_to_sequence[a] for a in virus_accessions}
        serum_sequences = {accession_to_sequence[a] for a in serum_accessions}
        if len(virus_sequences) > 1 or len(serum_sequences) > 1:
            ambiguous_rows.append({"reason": "name_maps_to_multiple_sequences", "virus": virus_name,
                                   "reference": serum_name, "hi_titre": hi_value,
                                   "virus_accessions": ";".join(sorted(virus_accessions)),
                                   "reference_accessions": ";".join(sorted(serum_accessions))})
            continue

        virus_accession = sorted(virus_accessions)[0]
        serum_accession = sorted(serum_accessions)[0]
        rows_out.append({
            "virus": virus_accession,
            "reference": serum_accession,
            "hi_titre": hi_value,
            "virus_name": virus_name,
            "reference_name": serum_name,
            "virus_match_rule": virus_rule,
            "reference_match_rule": serum_rule,
        })

    fieldnames = list(rows_out[0].keys()) if rows_out else [
        "virus", "reference", "hi_titre", "virus_name", "reference_name",
        "virus_match_rule", "reference_match_rule"]
    with (output / "matched_pairs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    with (output / "unmatched_rows.csv").open("w", newline="") as handle:
        fields = ["reason", "virus", "reference", "hi_titre", "virus_rule", "reference_rule"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(mismatch_rows)

    with (output / "ambiguous_rows.csv").open("w", newline="") as handle:
        fields = ["reason", "virus", "reference", "hi_titre", "virus_accessions", "reference_accessions"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ambiguous_rows)

    with (output / "unmatched_virus_names.tsv").open("w") as handle:
        handle.write("name\tcount\n")
        for name, count in sorted(unmatched_virus_names.items(), key=lambda x: -x[1]):
            handle.write("%s\t%d\n" % (name, count))

    with (output / "unmatched_reference_names.tsv").open("w") as handle:
        handle.write("name\tcount\n")
        for name, count in sorted(unmatched_serum_names.items(), key=lambda x: -x[1]):
            handle.write("%s\t%d\n" % (name, count))

    report = {
        "fluantigen_rows": len(frame),
        "matched_rows": len(rows_out),
        "unmatched_rows": len(mismatch_rows),
        "ambiguous_rows": len(ambiguous_rows),
        "unique_unmatched_virus_names": len(unmatched_virus_names),
        "unique_unmatched_reference_names": len(unmatched_serum_names),
        "match_rule_counts": dict(match_rule_counts),
        "matched_rows_percent": round(100.0 * len(rows_out) / len(frame), 4) if len(frame) else 0.0,
    }
    (output / "CHECK_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
