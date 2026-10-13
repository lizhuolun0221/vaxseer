"""Merge H1N1 FASTA/XLS batches into a training-ready matched dataset."""
import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


EPI_RE = re.compile(r"^EPI\d+$")


def batch_number(path):
    match = re.search(r" (\d+)\.", path.name)
    return match.group(1) if match else "1"


def fasta_records(path):
    record_id = None
    description = None
    sequence = []
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, description, "".join(sequence)
                description = line[1:].strip()
                record_id = description.split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
        if record_id is not None:
            yield record_id, description, "".join(sequence)


def ha_id(description):
    values = [item.strip() for item in description.split("|")]
    matches = [item for item in values if EPI_RE.fullmatch(item)]
    return matches[-1] if matches else None


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def convert_xls(path, output_dir):
    subprocess.run(
        ["libreoffice", "--headless", "--convert-to", "csv", "--outdir", str(output_dir), str(path)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    return output_dir / (path.stem + ".csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fasta_dir = args.source / "HIN1" / "FASTA"
    meta_dir = args.source / "HIN1" / "META"
    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    fasta_by_id = {}
    fasta_sources = defaultdict(list)
    sequence_conflicts = []
    fasta_files = sorted(fasta_dir.glob("*.fasta"))
    for path in fasta_files:
        for _, description, sequence in fasta_records(path):
            accession = ha_id(description)
            if accession is None:
                continue
            fasta_sources[accession].append(path.name)
            previous = fasta_by_id.get(accession)
            if previous is None:
                fasta_by_id[accession] = (sequence, description, path.name)
            elif previous[0] != sequence:
                sequence_conflicts.append((accession, previous[2], path.name))

    with tempfile.TemporaryDirectory(prefix="vaxseer_h1n1_xls_") as temporary:
        csv_dir = Path(temporary)
        metadata_by_id = {}
        metadata_sources = defaultdict(list)
        metadata_duplicates = []
        metadata_files = sorted(meta_dir.glob("*.xls"))
        for path in metadata_files:
            converted = convert_xls(path, csv_dir)
            with converted.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle)
                fieldnames = reader.fieldnames or []
                for row in reader:
                    accession = (row.get("HA Segment_Id") or "").split("|")[0].strip()
                    if not accession:
                        continue
                    metadata_sources[accession].append(path.name)
                    if accession in metadata_by_id:
                        metadata_duplicates.append(accession)
                        continue
                    row["source_metadata_file"] = path.name
                    metadata_by_id[accession] = (row, fieldnames)

    matched = sorted(set(fasta_by_id) & set(metadata_by_id))
    fasta_only = sorted(set(fasta_by_id) - set(metadata_by_id))
    metadata_only = sorted(set(metadata_by_id) - set(fasta_by_id))
    fields = []
    for accession in matched:
        for field in metadata_by_id[accession][1] + ["source_metadata_file"]:
            if field not in fields:
                fields.append(field)

    with (output / "ha.fasta").open("w") as handle:
        for accession in matched:
            sequence, description, source = fasta_by_id[accession]
            handle.write(">%s|source_fasta_file=%s|original_header=%s\n%s\n" %
                         (description, source, description, sequence))

    with (output / "metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for accession in matched:
            row = dict(metadata_by_id[accession][0])
            for field in fields:
                row.setdefault(field, "")
            writer.writerow(row)

    with (output / "excluded.tsv").open("w") as handle:
        handle.write("reason\taccession_id\tsource_files\n")
        for accession in fasta_only:
            handle.write("fasta_without_metadata\t%s\t%s\n" %
                         (accession, ";".join(sorted(set(fasta_sources[accession])))))
        for accession in metadata_only:
            handle.write("metadata_without_fasta\t%s\t%s\n" %
                         (accession, ";".join(sorted(set(metadata_sources[accession])))))

    with (output / "source_manifest.tsv").open("w") as handle:
        handle.write("file_type\tfile\tsha256\n")
        for path in fasta_files + metadata_files:
            handle.write("%s\t%s\t%s\n" %
                         ("fasta" if path.suffix == ".fasta" else "metadata", path.name, sha256(path)))

    with (output / "sequence_conflicts.tsv").open("w") as handle:
        handle.write("accession_id\tfirst_file\tconflicting_file\n")
        for row in sorted(set(sequence_conflicts)):
            handle.write("%s\t%s\t%s\n" % row)

    report = {
        "source": str(args.source),
        "output": str(output),
        "fasta_files": len(fasta_files),
        "metadata_files": len(metadata_files),
        "fasta_records_with_ha_id": sum(len(values) for values in fasta_sources.values()),
        "unique_fasta_ha": len(fasta_by_id),
        "unique_metadata_ha": len(metadata_by_id),
        "matched_ha": len(matched),
        "fasta_without_metadata": len(fasta_only),
        "metadata_without_fasta": len(metadata_only),
        "duplicate_metadata_rows": len(metadata_duplicates),
        "sequence_conflict_rows": len(sequence_conflicts),
        "rule": "retain only HA EPI IDs with both FASTA and metadata; keep one row per HA EPI ID",
    }
    (output / "CHECK_REPORT.json").write_text(__import__("json").dumps(report, indent=2) + "\n")
    print(__import__("json").dumps(report, indent=2))


if __name__ == "__main__":
    main()