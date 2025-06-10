import os
import subprocess
import csv
from pathlib import Path
import argparse
from Bio import SeqIO

# --- Configuration ---
LIPOP_SCRIPT = "/opt/Inovirus/LipoP1.0a/LipoP"  # Replace with your LipoP.pl path

# --- Run LipoP using Perl ---
def run_lipop_on_fasta(fasta_path, output_path):
    """Run LipoP and save short output."""
    with open(output_path, "w") as out:
        subprocess.run(
            ["perl", LIPOP_SCRIPT, fasta_path, "-short"],
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=False
        )

# --- Parse LipoP short output for SpII (lipoprotein) hits ---
def extract_spii_hits(short_output_path):
    """Return (True/False, list of SpII protein IDs)."""
    names = []
    with open(short_output_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                parts = line.strip().split()
                if len(parts) > 2 and parts[2] == "SpII":
                    names.append(parts[1])
    return (len(names) > 0), names

# --- Extract matching sequences to new FASTA file ---
def extract_lipoproteins_from_fasta(input_faa, output_faa, keep_ids):
    """Write only sequences with matching IDs to a new FASTA file."""
    with open(input_faa, "r") as infile, open(output_faa, "w") as outfile:
        records = (rec for rec in SeqIO.parse(infile, "fasta") if rec.id in keep_ids)
        SeqIO.write(records, outfile, "fasta")

# --- Main batch function ---
def batch_run_lipop(root_dir, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "lipop_summary.tsv"
    summary_data = []

    refined_faa_files = list(Path(root_dir).rglob("*-refined.faa"))
    print(f"[INFO] Found {len(refined_faa_files)} -refined.faa files under {root_dir}\n")

    for i, faa_path in enumerate(refined_faa_files, start=1):
        rel_path = faa_path.relative_to(root_dir)
        base_name = faa_path.stem
        short_output_path = output_dir / f"{base_name}.short.txt"

        print(f"[{i}/{len(refined_faa_files)}] Processing: {rel_path}")
        run_lipop_on_fasta(str(faa_path), str(short_output_path))

        has_spii, spii_ids = extract_spii_hits(short_output_path)

        if has_spii:
            print(f"    → Lipoproteins found: {', '.join(spii_ids)}")
            filtered_faa_path = output_dir / f"{base_name}_lipoP.faa"
            extract_lipoproteins_from_fasta(faa_path, filtered_faa_path, spii_ids)
        else:
            print(f"    → No SpII lipoproteins detected.")

        summary_data.append({
            "File": str(rel_path),
            "Lipoprotein_Detected": "Yes" if has_spii else "No",
            "Lipoprotein_Name": ";".join(spii_ids) if spii_ids else ""
        })

    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["File", "Lipoprotein_Detected", "Lipoprotein_Name"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(summary_data)

    print(f"\n Summary written to: {summary_path}")

# --- CLI interface ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch run LipoP on -refined.faa files in subdirectories.")
    parser.add_argument("input_dir", help="Root directory containing subfolders with -refined.faa files")
    parser.add_argument("output_dir", help="Directory to store LipoP short outputs, extracted sequences, and summary")

    args = parser.parse_args()
    batch_run_lipop(args.input_dir, args.output_dir)
