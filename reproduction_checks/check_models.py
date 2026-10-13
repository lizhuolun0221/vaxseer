"""Strictly load official checkpoints and compare CPU/GPU reference forward passes.

This is an execution check, not a historical evaluation or a candidate ranking.
Run with vaxseer_lzl; no training, sequence generation, or checkpoint mutation.
"""
import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def unique_checkpoint(pattern):
    matches = sorted(ROOT.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError("Expected one checkpoint for %s, found %d" % (pattern, len(matches)))
    return matches[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subtype", choices=["a_h1n1", "a_h3n2"], default="a_h1n1")
    parser.add_argument("--year", type=int, default=2018)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit("Use a new output directory; existing results are preserved.")
    output.mkdir(parents=True)
    # The tiny GPT-2 config must already be cached. Do not fetch weights or configs here.
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    os.environ["HF_HUB_OFFLINE"] = "1"
    os.chdir(ROOT / "vaxseer")  # Official HI checkpoint stores a relative model-args path.
    sys.path.insert(0, str(ROOT / "vaxseer"))
    import torch
    from models import build_model
    from data.datasets import PairwiseAlnClassificationDataset
    from data.utils import default_lm_collate_func

    torch.manual_seed(1005)
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    refs = {
        "a_h1n1": "protein_A_California_07_2009_H1N1_ha.fasta",
        "a_h3n2": "prortein_A_NewYork_392_2004_H3N2_ha.fasta",
    }
    reference = ROOT / "data/reference_fasta" / refs[args.subtype]
    lines = reference.read_text().splitlines()
    if sum(line.startswith(">") for line in lines) != 1:
        raise RuntimeError("Reference file must contain exactly one sequence")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    report = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "executable": sys.executable,
        "torch": torch.__version__,
        "cuda_build": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "subtype": args.subtype,
        "year": args.year,
        "reference": str(reference.relative_to(ROOT)),
        "reference_sha256": digest(reference),
        "sequence_length": len(sequence),
        "sequence_x_count": sequence.count("X"),
        "scope": "strict loading and FP32 forward smoke check; not accuracy validation",
        "models": {},
    }
    if not args.cpu_only and (not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count()):
        raise RuntimeError("Requested GPU unavailable; explicitly use --cpu-only for CPU checks")
    gpu = torch.device("cuda:%d" % args.gpu)
    if not args.cpu_only:
        report["gpu"] = {"index": args.gpu, "name": torch.cuda.get_device_name(gpu)}
    patterns = {
        "dominance": "runs/flu_lm/2003-10_to_%d-02_2M/%s/**/*.ckpt" % (args.year, args.subtype),
        "antigenicity": "runs/flu_hi_msa_regressor/before_%d-02/%s_seed=1005/**/*.ckpt" % (args.year, args.subtype),
    }
    classes = {"dominance": "gpt2_time_new", "antigenicity": "esm_regressor"}
    for task in patterns:
        started = time.monotonic()
        item = {"ok": False}
        report["models"][task] = item
        try:
            path = unique_checkpoint(patterns[task])
            item.update(checkpoint=str(path.relative_to(ROOT)), checkpoint_sha256=digest(path))
            print("Loading %s strictly on CPU" % task, flush=True)
            model = build_model(classes[task]).load_from_checkpoint(str(path), map_location="cpu", strict=True)
            model.eval()
            item.update(strict_load=True, parameter_count=sum(p.numel() for p in model.parameters()))
            if task == "dominance":
                # Same two-month time bins as the official 2018 evaluation pipeline (90..92).
                time_min = (args.year - 2004) * 6 + 6
                time_max = time_min + 2
                times = list(range(time_min, time_max + 1))
                records = [{"src_seq": sequence, "src_time": float(t)} for t in times]
                batch = default_lm_collate_func(records, model.alphabet.get_batch_converter(max_positions=1024), model.pad_idx)
                item["time_bins"] = times

                def forward(batch):
                    return model.forward(batch, None, reduce=False, mode="test")[0]
            else:
                # An identical reference pair is already aligned; use the official dataset/collator.
                dataset = PairwiseAlnClassificationDataset(
                    [("reference", "reference", sequence, sequence, 0.0)], model.alphabet,
                    category=False,
                    prepend_special_token_for_seq1=model.config.prepend_special_token_for_virus,
                    prepend_special_token_for_seq2=model.config.prepend_special_token_for_vaccine,
                )
                batch = default_lm_collate_func([dataset[0]], model.alphabet.get_batch_converter(max_positions=1024), model.pad_idx)

                def forward(batch):
                    return model.forward(batch, None, mode="test")[model.config.labels[0]]
            item["input_shape"] = list(batch["input_ids"].shape)
            with torch.no_grad():
                cpu_result = forward(batch).detach().cpu()
            item["cpu"] = {"shape": list(cpu_result.shape), "finite": bool(torch.isfinite(cpu_result).all()), "output": cpu_result.tolist()}
            if not item["cpu"]["finite"]:
                raise RuntimeError("CPU output contains NaN or Inf")
            if not args.cpu_only:
                model.to(gpu)
                gpu_batch = {k: v.to(gpu) if torch.is_tensor(v) else v for k, v in batch.items()}
                torch.cuda.reset_peak_memory_stats(gpu)
                with torch.no_grad():
                    gpu_result = forward(gpu_batch).detach().cpu()
                torch.cuda.synchronize(gpu)
                close = bool(torch.allclose(cpu_result, gpu_result, atol=1e-3, rtol=1e-4))
                item["gpu"] = {"finite": bool(torch.isfinite(gpu_result).all()), "output": gpu_result.tolist(), "max_abs_difference": float((cpu_result - gpu_result).abs().max()), "cpu_gpu_close": close, "atol": 1e-3, "rtol": 1e-4, "peak_memory_bytes": torch.cuda.max_memory_allocated(gpu)}
                if not item["gpu"]["finite"] or not close:
                    raise RuntimeError("GPU finite-output or CPU/GPU agreement check failed")
                del gpu_batch
            item["ok"] = True
            del model, batch
            gc.collect()
            if not args.cpu_only:
                torch.cuda.empty_cache()
            print("PASS %s (%.1fs)" % (task, time.monotonic() - started), flush=True)
        except Exception:
            item["error"] = traceback.format_exc()
            print(item["error"], flush=True)
        item["elapsed_seconds"] = round(time.monotonic() - started, 3)
        (output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    report["ok"] = all(item["ok"] for item in report["models"].values())
    (output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print("Result: %s" % (output / "summary.json"), flush=True)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
