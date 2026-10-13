"""Cross-check reconstructed labels against official code and prepare aligned HI inputs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(command, cwd, log):
    print("Running", log.name, flush=True)
    with log.open("w") as f:
        f.write(json.dumps(command) + "\n")
        f.flush()
        env = dict(os.environ, PYTHONHASHSEED="1005", TRANSFORMERS_OFFLINE="1", HF_HUB_OFFLINE="1")
        result = subprocess.run(command, cwd=str(cwd), stdout=f, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise RuntimeError("Command failed; inspect %s" % log)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--dataset", type=Path, required=True)
    p.add_argument("--run-prediction", action="store_true")
    p.add_argument("--gpu", type=int, default=1)
    args = p.parse_args()
    audit, dataset = args.audit.resolve(), args.dataset.resolve()
    evidence = json.loads((audit / "summary.json").read_text())
    for name, key in [("ha.fasta", "ha_fasta"), ("metadata.csv", "metadata")]:
        if digest(dataset / name) != evidence["source_sha256"][key]:
            raise RuntimeError("Dataset does not match the audit")
    raw_hi = ROOT / "data/antigenicity/hi_processed/a_h1n1_pairs.csv"
    if digest(raw_hi) != evidence["source_sha256"]["hi_pairs"]:
        raise RuntimeError("Raw HI source does not match the audit")
    if digest(ROOT / "data/gisaid/ha_acc_ids.csv") != evidence["source_sha256"]["paper_ids"]:
        raise RuntimeError("Paper ID source does not match the audit")
    for name, expected in evidence["output_sha256"].items():
        if digest(audit / name) != expected:
            raise RuntimeError("Audit output changed: %s" % name)
    work = audit / "input_checks"
    if work.exists():
        raise SystemExit("Checks already exist; preserve results and use a new audit directory")
    work.mkdir()
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    paper = {r["accession_id(ha)"] for r in csv.DictReader((ROOT / "data/gisaid/ha_acc_ids.csv").open())}
    meta = {r["HA Segment_Id"]: r for r in csv.DictReader((dataset / "metadata.csv").open())}
    cohort = work / "official_collection_cohort.fasta"
    cohort_count = 0
    with cohort.open("w") as f:
        for record in SeqIO.parse(str(dataset / "ha.fasta"), "fasta"):
            accession = record.description.split("|")[7]
            collection = meta[accession]["Collection_Date"]
            if evidence["paper_only"] and accession not in paper:
                continue
            if len(record.seq) >= 553 and len(collection) in {7, 10} and collection[:7] < evidence["cutoff"][:7]:
                SeqIO.write([record], f, "fasta")
                cohort_count += 1
    if cohort_count != evidence["selection"]["eligible_collection_only"]:
        raise RuntimeError("Cross-check cohort count does not match")
    run([sys.executable, str(ROOT / "data/antigenicity/prepare.py"), "--pairs_path", str(raw_hi),
         "--sequences_path", str(cohort), "--max_date_exclude", evidence["cutoff"][:7],
         "--output_root_dir", str(work / "official_output")], ROOT, work / "official_prepare.log")
    original_dir = work / "official_output" / ("before_" + evidence["cutoff"][:7])

    def load_semantic_pairs(table, fasta):
        seqs = {r.id: hashlib.sha256(str(r.seq).encode()).hexdigest() for r in SeqIO.parse(str(fasta), "fasta")}
        frame = pd.read_csv(table)
        result = {(seqs[r.virus], seqs[r.reference]): r.hi for r in frame.itertuples()}
        if len(result) != len(frame):
            raise RuntimeError("Duplicate sequence-pair keys")
        return result

    official = load_semantic_pairs(original_dir / "a_h1n1_hi_folds.csv", original_dir / "a_h1n1.fasta")
    reconstructed = load_semantic_pairs(audit / "collection_only/official_compatible_hi_folds.csv",
                                         audit / "collection_only/sequences.fasta")
    if set(official) != set(reconstructed):
        raise RuntimeError("Sequence-pair keys differ from official output")
    differences = [abs(official[k] - reconstructed[k]) for k in official]
    if max(differences, default=0) > 1e-10:
        raise RuntimeError("Labels differ from official output")
    checks = {"official_cohort_accessions": cohort_count, "official_sequence_pair_keys_match": True,
              "official_sequence_pair_count": len(official), "official_max_abs_label_difference": max(differences, default=0),
              "scope": "reconstruction and input validation only; no held-out accuracy claim"}
    print("Official labels agree:", len(official), "pairs", flush=True)
    target = audit / "submission_checked"
    pairs = pd.read_csv(target / "unique_sequence_hi_folds.csv")
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(str(target / "sequences.fasta"), "fasta")}
    for role in ["virus", "reference"]:
        with (work / (role + ".fasta")).open("w") as f:
            for accession in sorted(set(pairs[role])):
                f.write(">%s\n%s\n" % (accession, sequences[accession]))
    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise RuntimeError("MMseqs2 unavailable")
    run([mmseqs, "easy-search", str(work / "reference.fasta"), str(work / "virus.fasta"),
         str(work / "align.m8"), str(work / "mmseqs_tmp"), "--threads", "4", "--format-output",
         "query,target,qaln,taln,qstart,qend,tstart,tend,mismatch", "--max-seqs", "5000"],
        ROOT, work / "alignment.log")
    aligns = {}
    with (work / "align.m8").open() as f:
        for line in f:
            ref, virus, qaln, taln, *_ = line.split()
            if len(qaln) != len(taln):
                raise RuntimeError("Unequal aligned sequence lengths")
            if (virus, ref) in aligns:
                raise RuntimeError("Duplicated MMseqs2 alignment keys")
            aligns[(virus, ref)] = (taln, qaln)
    ready, missing = [], []
    for row in pairs.itertuples():
        if (row.virus, row.reference) not in aligns:
            missing.append({"virus": row.virus, "reference": row.reference})
            continue
        virus_seq, reference_seq = aligns[(row.virus, row.reference)]
        if virus_seq.replace("-", "") not in sequences[row.virus] or reference_seq.replace("-", "") not in sequences[row.reference]:
            raise RuntimeError("Alignment orientation/source verification failed")
        ready.append({"virus": row.virus, "reference": row.reference, "virus_seq": virus_seq,
                      "reference_seq": reference_seq, "hi": row.hi})
    aligned = work / "aligned_unique_sequence_pairs.csv"
    pd.DataFrame(ready).to_csv(aligned, index=False)
    pd.DataFrame(missing, columns=["virus", "reference"]).to_csv(work / "missing_alignments.csv", index=False)
    checks.update(required_pairs=len(pairs), aligned_pairs=len(ready), missing_alignments=len(missing),
                  virus_sequence_count=pairs.virus.nunique(), reference_sequence_count=pairs.reference.nunique())
    (work / "summary.json").write_text(json.dumps(checks, indent=2) + "\n")
    if missing:
        raise RuntimeError("Some required pairs have no alignment; inspect missing_alignments.csv")
    if args.run_prediction:
        import torch
        if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
            raise RuntimeError("Requested GPU unavailable")
        frame = pd.DataFrame(ready)
        # A deliberately balanced format check, not a representative evaluation sample.
        sample = pd.concat([frame[frame.hi < 3].head(32), frame[frame.hi >= 3].head(32)], ignore_index=True)
        if len(sample) != 64:
            raise RuntimeError("Insufficient rows in either label class for the 64-pair input check")
        sample_path = work / "diagnostic_64_pairs.csv"
        sample.to_csv(sample_path, index=False)
        checkpoint_paths = list(ROOT.glob("runs/flu_hi_msa_regressor/before_2018-02/a_h1n1_seed=1005/**/*.ckpt"))
        if evidence["cutoff"] != "2018-02-01" or len(checkpoint_paths) != 1:
            raise RuntimeError("Prediction check currently supports only the verified H1N1-2018 model")
        run([sys.executable, "-m", "bin.train", "--default_root_dir", str(work / "prediction_logs"),
             "--data_module", "hi_regression_aln", "--model", "esm_regressor", "--accelerator", "gpu",
             "--devices", "%d," % args.gpu, "--batch_size", "16", "--num_workers", "0", "--precision", "32",
             "--seed", "1005", "--max_epochs", "1", "--predict", "--category", "false",
             "--resume_from_checkpoint", str(checkpoint_paths[0]), "--predict_index_path", str(sample_path),
             "--save_prediction_path", str(work / "diagnostic_64_predictions.csv")],
            ROOT / "vaxseer", work / "prediction.log")
        pred = pd.read_csv(work / "diagnostic_64_predictions.csv")
        if len(pred) != len(sample) or list(zip(pred.virus, pred.reference)) != list(zip(sample.virus, sample.reference)):
            raise RuntimeError("Predictions do not preserve input pair order")
        if not np.isfinite(pred.label).all() or not np.allclose(pred.hi, sample.hi):
            raise RuntimeError("Invalid output or labels were changed")
        checks.update(real_hi_prediction_check_rows=len(pred), predictions_finite=True, labels_preserved=True,
                      diagnostic_sampling="32 source-order pairs with hi<3 and 32 with hi>=3; not representative",
                      held_out_status="unverified: reconstructed pairs may overlap checkpoint training/validation/test data")
    checks["ok"] = True
    checks["output_sha256"] = {p.name: digest(p) for p in [aligned, work / "missing_alignments.csv"]}
    (work / "summary.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2), flush=True)


if __name__ == "__main__":
    main()
