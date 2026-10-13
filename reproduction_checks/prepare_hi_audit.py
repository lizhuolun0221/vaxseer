"""Audit H1N1 HI-to-HA joins and reconstruct deterministic sequence-pair labels.

Outputs are diagnostic inputs, not the original test set. No fuzzy matching,
random split, training, or biological sequence design is performed.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import date, datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize(name):
    return name.replace(" ", "_").replace("-", "_").lower()


def resolve(query, available_names):
    """Follow official alias rules, but expose multiple substring candidates."""
    if query in available_names:
        return [query], "exact"
    converted = query
    if len(converted.split("/")[-1]) == 2:
        converted = "/".join(converted.split("/")[:-1]) + "/20" + converted.split("/")[-1]
    if "(" in converted:
        groups = re.findall(r"[(](.*?)[)]", converted)
        if not groups:
            return [], "malformed_parenthesis"
        short = groups[0]
        if len(short.split("/")[-1]) == 2 and short.split("/")[-1].isdigit():
            century = "20" if int(short.split("/")[-1][0]) <= 2 else "19"
            short = "/".join(short.split("/")[:-1]) + "/" + century + short.split("/")[-1]
        prefix = converted.split("(")[0]
        candidates = sorted(n for n in available_names if short in n and prefix in n)
        if candidates:
            return candidates, "parenthetical_prefix"
        return ([short] if short in available_names else []), "parenthetical_short"
    return ([converted] if converted in available_names else []), "short_year" if converted != query else "unmatched"


def valid_submission(value):
    try:
        return date.fromisoformat(value).isoformat() == value
    except ValueError:
        return False


def write_csv(path, rows, columns, delimiter=","):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter=delimiter)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", type=Path, required=True)
    p.add_argument("--cutoff", default="2018-02-01")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--paper-only", action="store_true")
    args = p.parse_args()
    cutoff = date.fromisoformat(args.cutoff).isoformat()
    if not cutoff.endswith("-01"):
        raise SystemExit("Use the first of a month, matching the official month cutoff")
    dataset, out = args.dataset.resolve(), args.output.resolve()
    if out.exists():
        raise SystemExit("Use a new output directory; existing results are preserved")
    from Bio import SeqIO
    source = ROOT / "data/antigenicity/hi_processed/a_h1n1_pairs.csv"
    paper_path = ROOT / "data/gisaid/ha_acc_ids.csv"
    paper = {r["accession_id(ha)"] for r in csv.DictReader(paper_path.open())}
    config = json.loads((dataset / "build_config.json").read_text())
    for name in ["ha.fasta", "metadata.csv"]:
        if digest(dataset / name) != config["files"][name]["sha256"]:
            raise RuntimeError("Standard input changed: %s" % name)
    metadata = {r["HA Segment_Id"]: r for r in csv.DictReader((dataset / "metadata.csv").open())}
    all_names = defaultdict(lambda: defaultdict(list))
    pools = {"collection_only": defaultdict(lambda: defaultdict(list)),
             "submission_checked": defaultdict(lambda: defaultdict(list))}
    sequences = {}
    pool_accessions = {regime: defaultdict(list) for regime in pools}
    selection = Counter()
    seen = set()
    print("Building verified HA lookup tables", flush=True)
    for record in SeqIO.parse(str(dataset / "ha.fasta"), "fasta"):
        fields = record.description.split("|")
        accession = fields[7]
        if accession in seen:
            raise RuntimeError("Duplicate accession in standard input")
        seen.add(accession)
        m = metadata[accession]
        protein = str(record.seq)
        if args.paper_only and accession not in paper:
            selection["excluded_not_in_paper"] += 1
            continue
        if len(protein) < 553:
            selection["excluded_short_sequence"] += 1
            continue
        name = normalize(m["Isolate_Name"])
        h = hashlib.sha256(protein.encode()).hexdigest()
        sequences[h] = protein
        all_names[name][h].append(accession)
        collection = m["Collection_Date"]
        # Official HI preparation accepts month precision and rejects year-only dates.
        if len(collection) not in {7, 10} or collection[:7] >= cutoff[:7]:
            selection["excluded_collection_time_or_year_precision"] += 1
            continue
        for regime in pools:
            if regime == "submission_checked":
                submission = m["Submission_Date"]
                if not valid_submission(submission) or submission >= cutoff:
                    selection["excluded_submission_time_or_unknown"] += 1
                    continue
            pools[regime][name][h].append(accession)
            pool_accessions[regime][h].append(accession)
            selection["eligible_" + regime] += 1
    if seen != set(metadata):
        raise RuntimeError("FASTA and metadata IDs do not agree")
    rows = list(csv.DictReader(source.open()))
    raw_keys = Counter((r["virus"], r["reference"]) for r in rows)
    if any(n > 1 for n in raw_keys.values()):
        raise RuntimeError("Duplicate raw HI keys require an explicit aggregation rule")
    baseline = {}
    for r in rows:
        value = float(r["hi"])
        if not math.isfinite(value) or value <= 0:
            raise RuntimeError("HI values must be positive and finite")
        if r["virus"] == r["reference"]:
            baseline[r["reference"]] = value
    raw_names = sorted({r[key] for r in rows for key in ["virus", "reference"]})
    summary = {"built_at": datetime.now(timezone.utc).isoformat(), "scope": "HI join reconstruction; not original held-out evaluation",
               "cutoff": cutoff, "paper_only": args.paper_only, "raw_hi_rows": len(rows),
               "raw_virus_names": len({r["virus"] for r in rows}),
               "raw_reference_names": len({r["reference"] for r in rows}), "selection": dict(selection),
               "raw_hi_columns": list(rows[0]), "homologous_reference_count": len(baseline),
               "references_missing_homologous": sorted({r["reference"] for r in rows} - set(baseline)),
               "rows_missing_homologous": sum(r["reference"] not in baseline for r in rows),
               "normalization": "max(log2(HI(reference, reference)) - log2(HI(virus, reference)), 0)",
               "source_sha256": {"hi_pairs": digest(source), "ha_fasta": digest(dataset / "ha.fasta"),
                                 "metadata": digest(dataset / "metadata.csv"), "paper_ids": digest(paper_path)},
               "limits": ["HI CSV has no report/measurement dates", "original training/validation/test CSVs not recovered",
                          "current metadata cannot recover historical database revisions", "local sequence coverage is incomplete",
                          "conservative labels average CSV rows equally, unlike official accession-multiplicity weighting"],
               "regimes": {}}
    out.mkdir(parents=True)
    for regime, pool in pools.items():
        target = out / regime
        target.mkdir()
        name_maps = {}
        profile = Counter()
        for query in raw_names:
            candidates, method = resolve(query, pool)
            global_candidates, _ = resolve(query, all_names)
            hashes = {h for n in candidates for h in pool[n]}
            accessions = sorted({a for n in candidates for ids in pool[n].values() for a in ids})
            if not candidates:
                status = "time_filtered" if global_candidates else "missing_from_available_sequence"
            elif len(candidates) > 1:
                status = "ambiguous_alias"
            elif len(hashes) > 1:
                status = "multiple_sequences"
            else:
                status = "unique_sequence"
            profile[status] += 1
            name_maps[query] = {"query": query, "status": status, "method": method,
                                "matched_names": candidates, "sequence_hashes": sorted(hashes), "accessions": accessions}
        write_csv(target / "name_matching.tsv",
                  [{**r, "matched_names": "|".join(r["matched_names"]),
                    "sequence_hashes": "|".join(r["sequence_hashes"]), "accessions": "|".join(r["accessions"])} for r in name_maps.values()],
                  ["query", "status", "method", "matched_names", "sequence_hashes", "accessions"], "\t")
        missing = [r for r in name_maps.values() if r["status"] == "missing_from_available_sequence"]
        write_csv(target / "missing_names.tsv", [{"name": r["query"]} for r in missing], ["name"], "\t")
        pair_compatible = defaultdict(list)
        pair_unique = defaultdict(list)
        row_provenance = []
        pair_profile = Counter()
        for index, raw in enumerate(rows, 2):
            q, r = name_maps[raw["virus"]], name_maps[raw["reference"]]
            reason = ""
            if raw["reference"] not in baseline:
                reason = "missing_homologous"
            elif q["status"] in {"missing_from_available_sequence", "time_filtered", "ambiguous_alias"}:
                reason = "virus_" + q["status"]
            elif r["status"] in {"missing_from_available_sequence", "time_filtered", "ambiguous_alias"}:
                reason = "reference_" + r["status"]
            label = max(math.log2(baseline[raw["reference"]]) - math.log2(float(raw["hi"])), 0.0) if raw["reference"] in baseline else None
            is_unique = not reason and q["status"] == r["status"] == "unique_sequence"
            pair_profile[reason or "mapped_rows"] += 1
            if is_unique:
                pair_profile["unique_sequence_rows"] += 1
            row_provenance.append({"source_csv_row": index, "virus_name": raw["virus"], "reference_name": raw["reference"],
                                   "raw_hi": raw["hi"], "homologous_hi": baseline.get(raw["reference"], ""),
                                   "normalized_hi": label if label is not None else "", "reason": reason,
                                   "unique_sequence_pair": is_unique,
                                   "virus_sequence_hashes": "|".join(q["sequence_hashes"]),
                                   "reference_sequence_hashes": "|".join(r["sequence_hashes"])})
            if reason:
                continue
            for qh in q["sequence_hashes"]:
                for rh in r["sequence_hashes"]:
                    # Exact official prepare.py expansion weights: each accession combination contributes one label.
                    copies = len(pool[q["matched_names"][0]][qh]) * len(pool[r["matched_names"][0]][rh])
                    pair_compatible[(qh, rh)].append((index, label, copies))
                    if is_unique:
                        pair_unique[(qh, rh)].append((index, label, 1))
        write_csv(target / "raw_pair_matching.tsv", row_provenance, list(row_provenance[0]), "\t")
        canon = {h: min(ids) for h, ids in pool_accessions[regime].items()}
        used_hashes = {h for pair in pair_compatible for h in pair}
        with (target / "sequences.fasta").open("w") as f:
            for h in sorted(used_hashes, key=lambda h: canon[h]):
                f.write(">%s\n%s\n" % (canon[h], sequences[h]))
        for kind, pairs in [("official_compatible", pair_compatible), ("unique_sequence", pair_unique)]:
            output_rows = []
            provenance = []
            for (qh, rh), values in sorted(pairs.items(), key=lambda item: (canon[item[0][0]], canon[item[0][1]])):
                mass = sum(v[2] for v in values)
                value = sum(v[1] * v[2] for v in values) / mass
                output_rows.append({"virus": canon[qh], "reference": canon[rh], "hi": value})
                provenance.append({"virus": canon[qh], "reference": canon[rh], "virus_sha256": qh,
                                   "reference_sha256": rh, "source_rows": "|".join(str(v[0]) for v in values),
                                   "row_weights": "|".join(str(v[2]) for v in values), "total_weight": mass})
            write_csv(target / (kind + "_hi_folds.csv"), output_rows, ["virus", "reference", "hi"])
            write_csv(target / (kind + "_provenance.tsv"), provenance,
                      ["virus", "reference", "virus_sha256", "reference_sha256", "source_rows", "row_weights", "total_weight"], "\t")
        changed = sum(abs(sum(v[1] * v[2] for v in vals) / sum(v[2] for v in vals) -
                          sum(v[1] for v in vals) / len(vals)) > 1e-10 for vals in pair_compatible.values())
        summary["regimes"][regime] = {"names": dict(profile), "rows": dict(pair_profile),
                                      "official_compatible_sequence_pairs": len(pair_compatible),
                                      "unique_sequence_pairs": len(pair_unique), "sequence_count": len(used_hashes),
                                      "pairs_affected_by_accession_multiplicity_weighting": changed}
        print(regime, json.dumps(summary["regimes"][regime]), flush=True)
    summary["output_sha256"] = {str(p.relative_to(out)): digest(p) for p in out.rglob("*") if p.is_file()}
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    # A small inspectable notebook reads the saved evidence rather than rebuilding or fetching data.
    cells = [
        {"cell_type": "markdown", "metadata": {}, "source": ["# HI–HA 匹配复核\n", "本结果是本地映射重建，尚非原论文测试集。标签归一化、编号多重性和时间限制见 summary.json。"]},
        {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
         "source": ["from pathlib import Path\n", "import json, pandas as pd\n", "root = Path('.')\n", "summary = json.loads((root/'summary.json').read_text())\n", "summary['regimes']"]},
        {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
         "source": ["names = pd.read_csv(root/'submission_checked/name_matching.tsv', sep='\\t')\n", "names.groupby('status').size()"]},
        {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
         "source": ["pairs = pd.read_csv(root/'submission_checked/raw_pair_matching.tsv', sep='\\t')\n", "pairs.groupby('reason', dropna=False).size()"]},
        {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
         "source": ["import numpy as np\n", "eligible = pairs[pairs.homologous_hi.notna()]\n", "recomputed = np.maximum(np.log2(eligible.homologous_hi)-np.log2(eligible.raw_hi), 0)\n", "assert np.allclose(recomputed, eligible.normalized_hi)\n", "folds = pd.read_csv(root/'submission_checked/unique_sequence_hi_folds.csv')\n", "assert not folds.duplicated(['virus','reference']).any()\n", "assert np.isfinite(folds.hi).all()\n", "folds.hi.describe()"]},
    ]
    notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python (vaxseer_lzl)", "language": "python", "name": "python3"}},
                "nbformat": 4, "nbformat_minor": 4}
    (out / "hi_match_review.ipynb").write_text(json.dumps(notebook, ensure_ascii=False, indent=2) + "\n")
    print("Evidence saved to", out, flush=True)


if __name__ == "__main__":
    main()
