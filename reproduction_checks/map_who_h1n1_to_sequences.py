"""Map collected WHO H1N1 vaccine strain recommendations to local HA records."""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def normalize(name):
    return re.sub(r"[_-]", "_", name.strip().replace(" ", "_")).lower()


def accession(description):
    ids = [x.strip() for x in description.split("|") if re.fullmatch(r"EPI\d+", x.strip())]
    return ids[-1] if ids else ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recommendations", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--merged-dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    sequence_by_accession = {}
    for record in SeqIO.parse(str(args.merged_dataset / "ha.fasta"), "fasta"):
        acc = accession(record.description)
        if acc:
            sequence_by_accession[acc] = str(record.seq)

    name_to_accessions = defaultdict(set)
    with (args.merged_dataset / "metadata.csv").open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            acc = (row.get("HA Segment_Id") or "").split("|")[0].strip()
            name = row.get("Isolate_Name") or ""
            if acc and name:
                name_to_accessions[normalize(name)].add(acc)

    recommendations = []
    for path, name_column in ((args.recommendations, "A/H1N1"), (args.archive, "A/H1N1_recommendation")):
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                name = (row.get(name_column) or "").strip()
                if not name:
                    continue
                clean_name = re.sub(r"\s*\(H1N1\).*|\s*-like$", "", name, flags=re.IGNORECASE).strip()
                accs = sorted(name_to_accessions.get(normalize(clean_name), set()))
                seqs = {sequence_by_accession[a] for a in accs if a in sequence_by_accession}
                recommendations.append({
                    "season": row.get("Season", ""),
                    "hemisphere": row.get("Hemisphere", ""),
                    "platform": row.get("Platform", "Historical archive entry"),
                    "who_strain": name,
                    "normalized_strain": clean_name,
                    "exact_name_accessions": len(accs),
                    "unique_ha_sequences": len(seqs),
                    "accessions": ";".join(accs),
                    "sequence_lengths": ";".join(map(str, sorted({len(s) for s in seqs}))),
                    "mapping_status": "unmatched" if not accs else "unique_sequence" if len(seqs) == 1 else "ambiguous_sequences",
                    "source_url": row.get("Source_URL", ""),
                })

    fields = ["season", "hemisphere", "platform", "who_strain", "normalized_strain", "exact_name_accessions", "unique_ha_sequences", "accessions", "sequence_lengths", "mapping_status", "source_url"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(recommendations)
    print("mapped recommendation rows:", len(recommendations))
    print("unmatched:", sum(r["mapping_status"] == "unmatched" for r in recommendations))
    print("ambiguous sequences:", sum(r["mapping_status"] == "ambiguous_sequences" for r in recommendations))


if __name__ == "__main__":
    main()
