"""Separate truly absent local HA names from metadata/paper-list gaps, without inferring host."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sqlite3

from prepare_hi_audit import normalize, resolve, digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--sequence-index", type=Path, required=True)
    p.add_argument("--metadata-index", type=Path, required=True)
    args = p.parse_args()
    audit = args.audit.resolve()
    target = audit / "missing_triage"
    if target.exists():
        raise SystemExit("Preserve previous triage results; use a new audit directory")
    seq = sqlite3.connect(args.sequence_index.resolve().as_uri() + "?mode=ro", uri=True)
    meta = sqlite3.connect(args.metadata_index.resolve().as_uri() + "?mode=ro", uri=True)
    meta_ids = {r[0] for r in meta.execute("SELECT DISTINCT ha_id FROM ha_metadata")}
    paper = {r[0] for r in seq.execute("SELECT accession_id FROM expected_ha")}
    raw_names = defaultdict(dict)
    for accession, name, length in seq.execute("SELECT accession_id,strain_name,seq_length FROM records WHERE subtype='AH1N1' AND segment='HA'"):
        if length >= 553:
            raw_names[normalize(name)][accession] = True
    queries = [r["name"] for r in csv.DictReader((audit / "submission_checked/missing_names.tsv").open(), delimiter="\t")]
    rows = []
    counts = Counter()
    for query in queries:
        names, method = resolve(query, raw_names)
        ids = sorted({a for n in names for a in raw_names[n]})
        with_meta = [a for a in ids if a in meta_ids]
        paper_ids = [a for a in ids if a in paper]
        paper_missing_meta = [a for a in paper_ids if a not in meta_ids]
        if not ids:
            status = "no_local_length_eligible_HA_match"
        elif len(names) > 1:
            status = "ambiguous_raw_alias"
        elif paper_missing_meta:
            status = "local_paper_HA_needs_metadata"
        elif not paper_ids:
            status = "local_HA_outside_paper_ID_list"
        else:
            status = "local_paper_HA_needs_conflict_or_metadata_review"
        counts[status] += 1
        rows.append({"query": query, "status": status, "method": method, "raw_names": "|".join(names),
                     "local_ha_ids": "|".join(ids), "metadata_ha_ids": "|".join(with_meta),
                     "paper_ha_ids": "|".join(paper_ids), "paper_ids_needing_metadata": "|".join(paper_missing_meta)})
    target.mkdir()
    with (target / "triage.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    needed_meta = sorted({a for r in rows for a in r["paper_ids_needing_metadata"].split("|") if a})
    (target / "local_paper_ha_needing_metadata.txt").write_text("".join(a + "\n" for a in needed_meta))
    markers = ",".join("?" for _ in needed_meta)
    isolates = sorted({r[0] for r in seq.execute(
        "SELECT DISTINCT isolate_id FROM records WHERE segment=? AND subtype=? AND accession_id IN (" + markers + ")",
        ["HA", "AH1N1"] + needed_meta)}) if needed_meta else []
    (target / "local_paper_isolates_needing_metadata.txt").write_text("".join(a + "\n" for a in isolates))
    absent_names = [r["query"] for r in rows if r["status"] == "no_local_length_eligible_HA_match"]
    with (target / "names_needing_sequence_search.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t"); writer.writerow(["hi_name"])
        writer.writerows([[n] for n in absent_names])
    summary = {"query_count": len(queries), "counts": dict(counts), "metadata_ha_id_count": len(needed_meta),
               "metadata_isolate_id_count": len(isolates),
               "scope": "raw H1N1 HA name matching with length>=553; host/submission/date eligibility is not inferred",
               "sequence_index": str(args.sequence_index.resolve()), "metadata_index": str(args.metadata_index.resolve()),
               "query_sha256": digest(audit / "submission_checked/missing_names.tsv"),
               "limits": ["no match does not prove GISAID absence", "no match can also reflect alias or short-sequence issues",
                          "raw matches cannot replace verified metadata or resolve sequence ambiguity"]}
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
