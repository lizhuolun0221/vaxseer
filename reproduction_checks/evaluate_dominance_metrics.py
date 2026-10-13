"""Extended metrics for fixed-candidate dominance holdout predictions."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import pearsonr, spearmanr


def read_fasta(path):
    records = {}
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split()[1].split("|")
        attrs = {x.split("=", 1)[0]: x.split("=", 1)[1] for x in fields if "=" in x}
        records[record.id] = {
            "sequence": str(record.seq),
            "frequency": float(attrs.get("freq", 0.0)),
            "bin_size": int(float(attrs.get("bin_size", 0))),
        }
    return records


def ranked_prediction(df, id_to_sequence):
    losses = defaultdict(list)
    for row in df.itertuples():
        losses[row.src_id].append(float(row.prediction))
    ranked = []
    for sequence_id, values in losses.items():
        m = min(values)
        log_mean_prob = -m + math.log(sum(math.exp(-(v - m)) for v in values)) - math.log(len(values))
        ranked.append((id_to_sequence[sequence_id], log_mean_prob))
    return sorted(ranked, key=lambda x: x[1], reverse=True)


def ndcg(relevance, order, k):
    gains = np.asarray([relevance.get(s, 0.0) for s in order[:k]], dtype=float)
    discounts = np.log2(np.arange(2, len(gains) + 2))
    dcg = float(np.sum(gains / discounts))
    ideal = np.sort(np.asarray([relevance.get(s, 0.0) for s in order], dtype=float))[::-1][:k]
    ideal_discounts = np.log2(np.arange(2, len(ideal) + 2))
    idcg = float(np.sum(ideal / ideal_discounts))
    return dcg / idcg if idcg else 0.0


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--future-root", type=Path, required=True)
    parser.add_argument("--prediction-files", nargs="+", type=Path, required=True)
    parser.add_argument("--future-bins", nargs="+", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1005)
    args = parser.parse_args()
    if len(args.prediction_files) != len(args.future_bins):
        raise SystemExit("Each prediction file must correspond to one future bin")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidate = read_fasta(args.candidate_fasta)
    id_to_sequence = {key: row["sequence"] for key, row in candidate.items()}
    candidate_frequency = {row["sequence"]: row["frequency"] for row in candidate.values()}
    fractions = [0.05, 0.10, 0.20, 0.30, 0.50]
    summary_rows = []
    curve_rows = []
    bootstrap_rows = []
    rng = np.random.default_rng(args.seed)

    for prediction_file, future_bin in zip(args.prediction_files, args.future_bins):
        df = pd.read_csv(prediction_file)
        ranked_predictions = ranked_prediction(df, id_to_sequence)
        vaxseer_order = [s for s, _ in ranked_predictions]
        last_order = [s for s, _ in sorted(candidate_frequency.items(), key=lambda x: x[1], reverse=True)]
        if len(vaxseer_order) != len(candidate):
            raise ValueError("Prediction candidate count does not match candidate FASTA")
        times = sorted(df.src_time.unique().tolist())
        future_path = args.future_root / (str(future_bin) + ".fasta")
        future_records = read_fasta(future_path)
        future_frequency = {row["sequence"]: row["frequency"] for row in future_records.values()}
        freq_total = sum(future_frequency.values())
        candidate_support = sum(future_frequency.get(s, 0.0) for s in candidate_frequency)
        p = df.groupby("src_id").prediction.apply(lambda x: float(np.mean(np.exp(-x.to_numpy()))))
        observed = np.asarray([future_frequency.get(id_to_sequence[i], 0.0) for i in p.index], dtype=float)
        prob = p.to_numpy()
        rho = float(spearmanr(prob, observed).statistic) if np.unique(observed).size > 1 else float("nan")
        pearson = float(pearsonr(prob, observed).statistic) if np.std(observed) and np.std(prob) else float("nan")
        last_values = np.asarray([candidate_frequency[id_to_sequence[s]] for s in p.index], dtype=float)
        last_rho = float(spearmanr(last_values, observed).statistic) if np.unique(observed).size > 1 else float("nan")

        top20_count = int(len(candidate) * 0.20)
        v_top20 = vaxseer_order[:top20_count]
        l_top20 = last_order[:top20_count]
        top20_gain = sum(future_frequency.get(s, 0.0) for s in v_top20) - sum(future_frequency.get(s, 0.0) for s in l_top20)
        n_records = next(iter(future_records.values()))["bin_size"]
        observed_counts = np.asarray([round(future_frequency.get(row["sequence"], 0.0) * n_records) for row in future_records.values()], dtype=int)
        if observed_counts.sum() > 0:
            ids = list(future_records.keys())
            seqs = [future_records[i]["sequence"] for i in ids]
            counts_prob = observed_counts / observed_counts.sum()
            v_mask = np.asarray([s in set(v_top20) for s in seqs])
            l_mask = np.asarray([s in set(l_top20) for s in seqs])
            deltas = []
            for start in range(0, args.bootstrap, 100):
                size = min(100, args.bootstrap - start)
                samples = rng.multinomial(int(observed_counts.sum()), counts_prob, size=size)
                v_cov = samples[:, v_mask].sum(axis=1) / samples.sum(axis=1)
                l_cov = samples[:, l_mask].sum(axis=1) / samples.sum(axis=1)
                deltas.extend((v_cov - l_cov).tolist())
            ci = np.quantile(deltas, [0.025, 0.975]).tolist()
            positive_probability = float(np.mean(np.asarray(deltas) > 0))
        else:
            ci, positive_probability = [float("nan"), float("nan")], float("nan")

        for frac in fractions:
            k = int(len(candidate) * frac)
            for method, order in (("VaxSeer", vaxseer_order), ("Last", last_order)):
                selected = order[:k]
                coverage = sum(future_frequency.get(s, 0.0) for s in selected)
                curve_rows.append({
                    "future_bin": future_bin, "prediction_bins": ";".join(map(str, times)),
                    "fraction": frac, "method": method, "top_k": k,
                    "frequency_coverage": coverage,
                    "candidate_support_mass": candidate_support,
                    "random_expected_coverage_within_candidates": candidate_support * frac,
                    "enrichment_vs_random_within_candidates": coverage / (candidate_support * frac) if candidate_support else "",
                    "matched_sequences": sum(future_frequency.get(s, 0.0) > 0 for s in selected),
                    "ndcg": ndcg(future_frequency, order, k),
                })
        v_cov = sum(future_frequency.get(s, 0.0) for s in v_top20)
        l_cov = sum(future_frequency.get(s, 0.0) for s in l_top20)
        summary_rows.append({
            "future_bin": future_bin, "prediction_bins": ";".join(map(str, times)),
            "future_distinct_sequences": len(future_frequency), "future_frequency_sum": freq_total,
            "candidate_count": len(candidate), "candidate_future_support_mass": candidate_support,
            "vaxseer_top20_coverage": v_cov, "last_top20_coverage": l_cov,
            "top20_absolute_gain": v_cov - l_cov,
            "top20_bootstrap_95ci_low": ci[0], "top20_bootstrap_95ci_high": ci[1],
            "bootstrap_probability_vaxseer_gt_last": positive_probability,
            "spearman_predicted_probability_vs_future_frequency": rho,
            "spearman_historical_frequency_vs_future_frequency": last_rho,
            "pearson_predicted_probability_vs_future_frequency": pearson,
            "vaxseer_ndcg_top20": ndcg(future_frequency, vaxseer_order, top20_count),
            "last_ndcg_top20": ndcg(future_frequency, last_order, top20_count),
        })
        bootstrap_rows.append({
            "future_bin": future_bin, "replicates": args.bootstrap,
            "seed": args.seed, "ci_method": "multinomial bootstrap over observed future sequence counts",
            "delta_top20_coverage_ci_low": ci[0], "delta_top20_coverage_ci_high": ci[1],
            "probability_delta_gt_zero": positive_probability,
            "warning": "This captures finite count sampling only, not database coverage/submission bias or candidate uncertainty.",
        })

    write_csv(args.output_dir / "extended_metrics_by_window.csv", summary_rows)
    write_csv(args.output_dir / "topk_coverage_curves.csv", curve_rows)
    write_csv(args.output_dir / "bootstrap_uncertainty.csv", bootstrap_rows)
    print(json.dumps({"windows": summary_rows, "outputs": ["extended_metrics_by_window.csv", "topk_coverage_curves.csv", "bootstrap_uncertainty.csv"]}, indent=2))


if __name__ == "__main__":
    main()
