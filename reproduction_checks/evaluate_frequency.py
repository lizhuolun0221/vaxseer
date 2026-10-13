"""Evaluate retrospective dominance predictions by future sequence frequency."""
import argparse
import json
import math
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


def fasta_info(path):
    sequences = {}
    frequencies = {}
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split()[1].split("|")
        fields = {item.split("=", 1)[0]: item.split("=", 1)[1] for item in fields}
        sequences[record.id] = str(record.seq)
        frequencies[str(record.seq)] = float(fields["freq"])
    return sequences, frequencies


def rank_by_prediction(predictions, sequences):
    losses = defaultdict(list)
    for row in predictions.itertuples():
        losses[row.src_id].append(float(row.prediction))
    ranked = []
    for sequence_id, values in losses.items():
        log_probability = math.log(sum(math.exp(-value) for value in values))
        log_probability -= math.log(len(values))
        ranked.append((sequences[sequence_id], log_probability))
    return [sequence for sequence, _ in sorted(ranked, key=lambda item: item[1], reverse=True)]


def evaluate(ranked_sequences, future_frequencies, percentile):
    top_count = int(len(ranked_sequences) * percentile)
    return {
        "candidate_count": len(ranked_sequences),
        "top_count": top_count,
        "future_frequency_total": sum(future_frequencies.values()),
        "top_frequency_coverage": sum(
            future_frequencies.get(sequence, 0.0)
            for sequence in ranked_sequences[:top_count]
        ),
        "future_sequence_count": len(future_frequencies),
        "matched_top_sequence_count": sum(
            sequence in future_frequencies
            for sequence in ranked_sequences[:top_count]
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-fasta", required=True)
    parser.add_argument("--future-fasta", required=True)
    parser.add_argument("--prediction-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--percentile", type=float, default=0.20)
    args = parser.parse_args()

    sequences, historical_frequencies = fasta_info(args.candidate_fasta)
    _, future_frequencies = fasta_info(args.future_fasta)
    predictions = pd.read_csv(args.prediction_csv)
    vaxseer_ranked = rank_by_prediction(predictions, sequences)
    last_ranked = [sequence for sequence, _ in sorted(
        historical_frequencies.items(), key=lambda item: item[1], reverse=True
    )]
    result = {
        "percentile": args.percentile,
        "vaxseer": evaluate(vaxseer_ranked, future_frequencies, args.percentile),
        "last": evaluate(last_ranked, future_frequencies, args.percentile),
        "protocol": {
            "prediction_aggregation": "logsumexp(-prediction) - log(number_of_prediction_bins)",
            "future_denominator": "all sequences in the future FASTA bin",
            "unmatched_future_frequency": "zero for candidate sequences absent from future bin",
        },
    }
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()