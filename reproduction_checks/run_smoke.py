"""Run the official preprocessing and prediction entries on an arbitrary small HA subset."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(command, cwd, log):
    print("Running %s" % log.name, flush=True)
    with log.open("w") as handle:
        handle.write(json.dumps(command) + "\n")
        handle.flush()
        result = subprocess.run(command, cwd=str(cwd), stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Command failed (%d); inspect %s" % (result.returncode, log))


def checkpoint(pattern):
    matches = list(ROOT.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError("Expected one checkpoint: %s" % pattern)
    return str(matches[0].resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--year", type=int, default=2018)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--size", type=int, default=4)
    parser.add_argument("--resume-preprocessing", action="store_true",
                        help="Reuse a previously prepared cohort and preprocessing output only.")
    parser.add_argument("--official-exclude-self", action="store_true",
                        help="Reproduce the original omission of pairs with equal accession IDs.")
    args = parser.parse_args()
    dataset, out = args.dataset.resolve(), args.output.resolve()
    if args.size < 2:
        raise SystemExit("At least two sequences are required")
    if out.exists() and not args.resume_preprocessing:
        raise SystemExit("Output exists; use a new directory")
    if args.resume_preprocessing and (out / "pairs.csv").exists():
        raise SystemExit("Resume is restricted to preprocessing; predictions are not overwritten")
    out.mkdir(parents=True, exist_ok=True)
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    os.environ["HF_HUB_OFFLINE"] = "1"
    import torch
    from Bio import SeqIO
    import numpy as np
    import pandas as pd
    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("Requested GPU unavailable")
    cutoff = "%d-02-01" % args.year
    cohort = out / ("metadata_available_before_%d-02.csv" % args.year)
    processed = out / "processed"
    if not args.resume_preprocessing:
        counts = {"submitted_before_cutoff": 0, "collection_window_day_precision": 0}
        with (dataset / "metadata.csv").open(newline="") as source, cohort.open("w", newline="") as target:
            reader = csv.DictReader(source)
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                submission = row["Submission_Date"]
                collection = row["Collection_Date"]
                if len(submission) == 10 and submission < cutoff:
                    writer.writerow(row)
                    counts["submitted_before_cutoff"] += 1
                    if len(collection) == 10 and "%d-02-01" % (args.year - 3) <= collection < cutoff:
                        counts["collection_window_day_precision"] += 1
        (out / "cohort_selection.json").write_text(json.dumps({"cutoff": cutoff, "counts": counts,
            "scope": "current metadata with historical submission date; not historical database versions"}, indent=2) + "\n")
        run([sys.executable, str(ROOT / "data/process_fasta.py"), "--time_interval", "9999",
             "--start_date", "%d-02" % (args.year - 3), "--end_date", "%d-02" % args.year,
             "--subtype", "a_h1n1", "--host", "human", "--split_by", "month",
             "--meta_data_path", str(cohort), "--sequences_path", str(dataset / "ha.fasta"),
             "--save_dir", str(processed), "--remove_min_size", "1000", "--min_seq_length", "553",
             "--min_count", "5"], ROOT, out / "preprocess.log")
    if not cohort.is_file():
        raise RuntimeError("Prepared cohort missing")
    paths = list(processed.rglob("*.fasta"))
    if len(paths) != 1:
        raise RuntimeError("Expected one processed FASTA")
    records = list(SeqIO.parse(str(paths[0]), "fasta"))
    selected = records[:args.size]  # Source order, without ranking or optimized selection.
    if len(selected) != args.size or len({str(r.seq) for r in selected}) != args.size:
        raise RuntimeError("Insufficient distinct sequences")
    selected_path = out / "smoke_ha.fasta"
    SeqIO.write(selected, str(selected_path), "fasta")
    lm = checkpoint("runs/flu_lm/2003-10_to_%d-02_2M/a_h1n1/**/*.ckpt" % args.year)
    hi = checkpoint("runs/flu_hi_msa_regressor/before_%d-02/a_h1n1_seed=1005/**/*.ckpt" % args.year)
    config = {"started_at": datetime.now(timezone.utc).isoformat(), "scope": "execution smoke check; no accuracy or historical performance claim",
              "subtype": "a_h1n1", "year": args.year, "cutoff": cutoff, "dataset": str(dataset),
              "source_fasta_sha256": digest(dataset / "ha.fasta"), "source_metadata_sha256": digest(dataset / "metadata.csv"),
              "cohort_sha256": digest(cohort), "processed_sha256": digest(paths[0]),
              "processed_sequence_count": len(records), "selected_ids": [r.id for r in selected],
              "selection": "first N distinct-sequence records in official preprocessing output; no ranking",
              "include_self_pairs": not args.official_exclude_self, "lm_checkpoint": lm, "hi_checkpoint": hi,
              "precision": 32, "gpu": args.gpu}
    (out / "run_config.json").write_text(json.dumps(config, indent=2) + "\n")
    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise RuntimeError("MMseqs2 unavailable")
    run([mmseqs, "easy-search", str(selected_path), str(selected_path), str(out / "align.m8"),
         str(out / "mmseqs_tmp"), "--threads", "4", "--format-output",
         "query,target,qaln,taln,qstart,qend,tstart,tend,mismatch", "--max-seqs", "2000"], ROOT, out / "alignment.log")
    pairs_command = [sys.executable, "-m", "antigenicity.build_test_pairs", "--alignment_path", str(out / "align.m8"),
                     "--virus_fasta_path", str(selected_path), "--vaccine_fasta_path", str(selected_path),
                     "--save_path", str(out / "pairs.csv")]
    if not args.official_exclude_self:
        pairs_command.append("--include_self_pairs")
    run(pairs_command, ROOT / "data", out / "pairs.log")
    pairs = pd.read_csv(out / "pairs.csv")
    expected = {(q.id, t.id) for q in selected for t in selected
                if not args.official_exclude_self or q.id != t.id}
    if set(zip(pairs.virus, pairs.reference)) != expected or len(pairs) != len(expected):
        raise RuntimeError("Alignment/pair matrix incomplete or duplicated")
    common = ["--accelerator", "gpu", "--devices", "%d," % args.gpu, "--batch_size", "4",
              "--precision", "32", "--num_workers", "0", "--max_epochs", "1", "--seed", "1005"]
    run([sys.executable, "-m", "bin.train", "--default_root_dir", str(out / "hi_logs"),
         "--data_module", "hi_regression_aln", "--model", "esm_regressor", "--predict",
         "--resume_from_checkpoint", hi, "--predict_index_path", str(out / "pairs.csv"), "--category", "false",
         "--save_prediction_path", str(out / "antigenicity_predictions.csv")] + common,
        ROOT / "vaxseer", out / "antigenicity.log")
    time_min = (args.year - 2004) * 6 + 6
    run([sys.executable, "-m", "bin.train", "--default_root_dir", str(out / "lm_logs"),
         "--data_module", "lm_weighted", "--model", "gpt2_time_new", "--test",
         "--test_data_paths", str(selected_path), "--resume_from_checkpoint", lm,
         "--min_testing_time", str(time_min), "--max_testing_time", str(time_min + 2),
         "--save_prediction_path", str(out / "dominance_predictions.csv")] + common,
        ROOT / "vaxseer", out / "dominance.log")
    run([sys.executable, str(ROOT / "evaluation/pipeline/1_calc_vaccine_score.py"),
         "--prob_pred_path", str(out / "dominance_predictions.csv"),
         "--hi_pred_path", str(out / "antigenicity_predictions.csv"),
         "--all_sequences_path", str(dataset / "ha.fasta"), "--save_path", str(out / "coverage_scores.csv")],
        ROOT, out / "coverage.log")
    dom = pd.read_csv(out / "dominance_predictions.csv")
    antigen = pd.read_csv(out / "antigenicity_predictions.csv")
    scores = pd.read_csv(out / "coverage_scores.csv")
    assert len(dom) == args.size * 3 and len(antigen) == len(expected) and len(scores) == args.size
    assert set(dom.src_id) == set(r.id for r in selected)
    assert set(zip(antigen.virus, antigen.reference)) == expected
    assert not dom.duplicated(["src_id", "src_time"]).any()
    for column in [dom.prediction, antigen.label, scores.score]:
        assert np.isfinite(column.to_numpy()).all()
    # Check the official score calculation independently, matching its probability averaging and support normalization.
    weights = np.exp(-dom.prediction).groupby(dom.src_id).mean()
    weights /= weights.sum()
    comparisons = []
    for row in scores.itertuples():
        subset = antigen[antigen.reference == row.reference]
        support = subset.virus.map(weights).to_numpy()
        value = float(np.dot(support, subset.label.to_numpy()) / support.sum())
        assert np.isclose(value, row.score, rtol=1e-5, atol=1e-6)
        comparisons.append({"reference": row.reference, "covered_probability_mass": float(support.sum()),
                            "score_calculation_agrees": True})
    summary = {"ok": True, "scope": config["scope"], "processed_sequence_count": len(records),
               "selected_sequence_count": args.size, "pair_count": len(pairs), "dominance_rows": len(dom),
               "antigenicity_rows": len(antigen), "coverage_rows": len(scores), "finite_outputs": True,
               "include_self_pairs": not args.official_exclude_self, "score_checks": comparisons,
               "score_direction": "official score is a dominance-weighted predicted log2 HI titer reduction; lower is better",
               "output_sha256": {name: digest(out / name) for name in ["pairs.csv", "dominance_predictions.csv", "antigenicity_predictions.csv", "coverage_scores.csv"]}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
