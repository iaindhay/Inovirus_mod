#!/bin/bash

set -euo pipefail

# === Input arguments ===
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <hmm_dir> <genome_root_dir> <output_root_dir>"
  exit 1
fi

hmm_dir="$1"
genome_root="$2"
output_root="$3"

timestamp=$(date +"%Y-%m-%d_%H%M")
hmmer_output_dir="${output_root}/hmmer_${timestamp}"
mkdir -p "$hmmer_output_dir"

echo "=== Starting HMMER batch search ==="
echo "[+] HMM profile dir:      $hmm_dir"
echo "[+] Genome root dir:      $genome_root"
echo "[+] Output will go to:    $hmmer_output_dir"
echo "-----------------------------------------"

# === Iterate through genomes and run hmmsearch ===
for acc_dir in "$genome_root"/*; do
  if [ ! -d "$acc_dir" ]; then
    continue
  fi

  acc=$(basename "$acc_dir")

  for faa_file in "$acc_dir"/*-refined.faa; do
    [ -e "$faa_file" ] || continue

    faa_base=$(basename "$faa_file")
    fragment_id="${faa_base%.faa}"
    fragment_id="${fragment_id#${acc}_}"

    for hmm_file in "$hmm_dir"/*.hmm; do
      [ -e "$hmm_file" ] || continue
      hmm_name=$(basename "$hmm_file" .hmm)
      out_tbl="${hmmer_output_dir}/${acc}__${fragment_id}__${hmm_name}.tbl"

      echo "[+] hmmsearch: $hmm_name vs $faa_base"
      hmmsearch --tblout "$out_tbl" "$hmm_file" "$faa_file" > /dev/null
    done
  done
done

echo "-----------------------------------------"
echo "[✓] All hmmsearch jobs completed."

# === Run summary and sequence extractor ===
script_dir="$(dirname "$0")"
echo "[✓] Running summary extractor..."
python3 "$script_dir/hmm_summary_extractor_final3.py" --hmmer_dir "$hmmer_output_dir" --root_dir "$genome_root"

echo "[✓] Output directory: $hmmer_output_dir"
