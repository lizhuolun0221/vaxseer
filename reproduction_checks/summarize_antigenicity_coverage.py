"""Summarize frequency-weighted HI model predictions for candidate HA strains."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def future_frequency(path):
    result = {}
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split()[1].split("|")
        fields = {item.split("=", 1)[0]: item.split("=", 1)[1] for item in fields if "=" in item}
        result[record.id] = float(fields["freq"])
    return result


def summarize(scores, weights):
    total = float(weights.sum())
    if total <= 0:
        return {"future_frequency_support": 0.0, "mean_predicted_log2_drop": None,
                "coverage_le_1_log2": None, "coverage_le_2_log2": None,
                "coverage_le_3_log2": None}
    normalized = weights.to_numpy(dtype=float) / total
    values = scores.to_numpy(dtype=float)
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = normalized[order]
    return {
        "future_frequency_support": total,
        "mean_predicted_log2_drop": float(np.dot(normalized, values)),
        "coverage_le_1_log2": float(normalized[values <= 1.0].sum()),
        "coverage_le_2_log2": float(normalized[values <= 2.0].sum()),
        "coverage_le_3_log2": float(normalized[values <= 3.0].sum()),
        "weighted_median_log2_drop": float(sorted_values[np.searchsorted(np.cumsum(sorted_weights), 0.5, side="left")]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-csv", type=Path, required=True)
    parser.add_argument("--future-fasta", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    predictions = pd.read_csv(args.prediction_csv)
    if "pred_label" in predictions.columns:
        prediction_column = "pred_label"
    elif "label" in predictions.columns:
        prediction_column = "label"
    else:
        raise ValueError("HI prediction column pred_label/label is missing")
    predictions = predictions.rename(columns={"src_id1": "virus", "src_id2": "reference", prediction_column: "predicted_log2_drop"})

    frequencies = future_frequency(args.future_fasta)
    predictions["future_frequency"] = predictions["virus"].map(frequencies).fillna(0.0)
    # Combine multiple local records that share one future amino-acid sequence ID only once per candidate.
    manifest = pd.read_csv(args.candidate_manifest)
    rows = []
    variant_details = []
    for candidate_id, group in predictions.groupby("reference", sort=False):
        weights = group["future_frequency"].astype(float)
        metrics = summarize(group["predicted_log2_drop"], weights)
        rows.append({"candidate_id": candidate_id, **metrics})
        variant_details.append({"candidate_id": candidate_id, **metrics})

    result = pd.DataFrame(rows)
    result.to_csv(args.output_dir / "candidate_antigenicity_coverage.csv", index=False)
    # Method summaries preserve WHO sequence variants separately; no arbitrary variant is selected.
    manifest["candidate_id"] = manifest["candidate_id"].astype(str)
    result = result.merge(manifest, on="candidate_id", how="left")
    groups = {
        "VaxSeer_top1": result[result.candidate_id.str.startswith("VaxSeer_top1|")],
        "Last_top1": result[result.candidate_id.str.startswith("Last_top1|")],
        "WHO_egg_based_variants": result[result.candidate_id.str.startswith("WHO_egg_variant")],
        "WHO_cell_based_variants": result[result.candidate_id.str.startswith("WHO_cell_variant")],
    }
    method_rows = []
    for method, frame in groups.items():
        if frame.empty:
            continue
        for metric in ["mean_predicted_log2_drop", "coverage_le_1_log2", "coverage_le_2_log2", "coverage_le_3_log2"]:
            vals = frame[metric].dropna().to_numpy(dtype=float)
            method_rows.append({"method": method, "metric": metric, "variant_count": len(vals),
                                "min_across_variants": float(vals.min()), "median_across_variants": float(np.median(vals)),
                                "max_across_variants": float(vals.max()),
                                "interpretation": "log2 HI reduction / predicted frequency fraction within threshold; not clinical VE"})
    pd.DataFrame(method_rows).to_csv(args.output_dir / "method_antigenicity_coverage_summary.csv", index=False)

    report = {
        "prediction_rows": len(predictions),
        "unique_future_virus_ids": predictions.virus.nunique(),
        "unique_candidate_reference_ids": predictions.reference.nunique(),
        "predicted_value_column": prediction_column,
        "label_interpretation": "model-predicted log2 HI titer reduction relative to homologous reference titer; lower means closer antigenically",
        "thresholds": {"<=1": "within 2-fold", "<=2": "within 4-fold", "<=3": "within 8-fold"},
        "future_frequency_denominator": str(args.future_fasta),
        "caveat": "HI checkpoint was trained on the FluAntigen random split including records published through 2024; these future-window HI estimates are exploratory and are not an independent temporal HI validation. WHO local name mappings have multiple sequence variants, so ranges are reported instead of selecting a favorable variant.",
    }
    (args.output_dir / "ANTIGENICITY_COVERAGE_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
