"""Select one vaccine candidate using dominance-weighted predicted HI distance.

This is a retrospective model-selection analysis, not a clinical VE estimate.
The HI model's current checkpoint was trained on a random split including HI
records through 2024; keep that temporal-leakage caveat with all results.
"""
import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def fasta_records(path):
    result = {}
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split()[1].split("|") if len(record.description.split()) > 1 else []
        attrs = {x.split("=", 1)[0]: x.split("=", 1)[1] for x in fields if "=" in x}
        result[record.id] = {
            "sequence": str(record.seq),
            "frequency": float(attrs.get("freq", 0.0)),
        }
    return result


def normalized_dominance_weights(path):
    df = pd.read_csv(path)
    values = df.groupby("src_id").prediction.apply(lambda s: float(np.mean(np.exp(-s.to_numpy()))))
    values = values / values.sum()
    return values.to_dict()


def score_candidates(pair_predictions, weights, future_freq_by_id):
    rows = []
    for reference, group in pair_predictions.groupby("reference", sort=False):
        supported = group[group.virus.isin(weights)]
        if len(supported) == 0:
            continue
        w = supported.virus.map(weights).to_numpy(dtype=float)
        d = supported.predicted_hi_drop.to_numpy(dtype=float)
        denom = w.sum()
        observed_w = supported.virus.map(future_freq_by_id).fillna(0).to_numpy(dtype=float)
        observed_denom = observed_w.sum()
        row = {
            "candidate_id": reference,
            "forecast_support_mass": float(denom),
            "forecast_weighted_mean_log2_drop": float(np.dot(w, d) / denom),
            "forecast_coverage_le_1": float(w[d <= 1].sum() / denom),
            "forecast_coverage_le_2": float(w[d <= 2].sum() / denom),
            "forecast_coverage_le_3": float(w[d <= 3].sum() / denom),
            "future_observed_support_mass": float(observed_denom),
            "future_weighted_mean_log2_drop": float(np.dot(observed_w, d) / observed_denom) if observed_denom else None,
            "future_coverage_le_1": float(observed_w[d <= 1].sum() / observed_denom) if observed_denom else None,
            "future_coverage_le_2": float(observed_w[d <= 2].sum() / observed_denom) if observed_denom else None,
            "future_coverage_le_3": float(observed_w[d <= 3].sum() / observed_denom) if observed_denom else None,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dominance-predictions", type=Path, required=True)
    parser.add_argument("--hi-predictions", type=Path, required=True)
    parser.add_argument("--future-fasta", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--virus-candidate-fasta", type=Path, required=True)
    parser.add_argument("--who-mapping", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dominance_weights = normalized_dominance_weights(args.dominance_predictions)
    hi = pd.read_csv(args.hi_predictions).rename(columns={"label": "predicted_hi_drop"})
    manifest = pd.read_csv(args.candidate_manifest).fillna("")
    future = fasta_records(args.future_fasta)
    virus_candidates = fasta_records(args.virus_candidate_fasta)
    future_frequency_by_sequence = {row["sequence"]: row["frequency"] for row in future.values()}
    future_freq_by_id = {
        accession: future_frequency_by_sequence.get(row["sequence"], 0.0)
        for accession, row in virus_candidates.items()
    }

    scores = score_candidates(hi, dominance_weights, future_freq_by_id)
    scores = scores.merge(manifest, left_on="candidate_id", right_on="candidate_id", how="left")
    scores["candidate_type"] = scores.candidate_type.fillna("unknown")
    scores.to_csv(args.output_dir / "all_vaccine_candidate_scores.csv", index=False)

    historical = scores[scores.candidate_type == "historical_pool"]
    if historical.empty:
        raise SystemExit("No historical-pool candidate scores found")
    selected_vaxseer = historical.sort_values("forecast_weighted_mean_log2_drop").iloc[0]

    historical_frequencies = {k: v["frequency"] for k, v in virus_candidates.items()}
    total = sum(historical_frequencies.values())
    last_weights = {k: v / total for k, v in historical_frequencies.items()}
    last_scores = score_candidates(hi, last_weights, future_freq_by_id)
    last_scores = last_scores.merge(manifest, on="candidate_id", how="left")
    selected_last = last_scores[last_scores.candidate_type == "historical_pool"].sort_values("forecast_weighted_mean_log2_drop").iloc[0]

    who = scores[scores.candidate_type.str.startswith("WHO_")]
    report = {
        "analysis": "Select one vaccine candidate by minimizing dominance-weighted predicted log2 HI drop; compare to Last-weighted selected candidate and WHO local sequence variants.",
        "future_window": "2024-10 <= Collection_Date < 2025-04",
        "historical_candidates_scored": int(len(historical)),
        "who_local_variants_scored": int(len(who)),
        "vaxseer_selected_candidate": selected_vaxseer.to_dict(),
        "last_selected_candidate": selected_last.to_dict(),
        "who_variant_range": {
            "forecast_weighted_mean_log2_drop_min": float(who.forecast_weighted_mean_log2_drop.min()),
            "forecast_weighted_mean_log2_drop_median": float(who.forecast_weighted_mean_log2_drop.median()),
            "forecast_weighted_mean_log2_drop_max": float(who.forecast_weighted_mean_log2_drop.max()),
            "future_weighted_mean_log2_drop_min": float(who.future_weighted_mean_log2_drop.min()),
            "future_weighted_mean_log2_drop_median": float(who.future_weighted_mean_log2_drop.median()),
            "future_weighted_mean_log2_drop_max": float(who.future_weighted_mean_log2_drop.max()),
        },
        "caveats": [
            "Lower predicted log2 HI drop is treated as lower predicted antigenic distance; this is not clinical vaccine effectiveness.",
            "HI checkpoint was trained on a random split with HI reports through 2024, so this is exploratory rather than a leakage-free temporal validation.",
            "WHO strain-like recommendations map to multiple local HA variants; all are reported, not a cherry-picked single sequence.",
            "Last selection uses historical candidate frequencies as weights; it is a reconstructed baseline, not the missing official dominance_from_fasta.py implementation.",
        ],
    }
    (args.output_dir / "vaccine_selection_report.json").write_text(json.dumps(report, indent=2, default=float) + "\n")
    print(json.dumps(report, indent=2, default=float))


if __name__ == "__main__":
    main()
