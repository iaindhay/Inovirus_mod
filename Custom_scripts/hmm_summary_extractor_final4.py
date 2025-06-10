import os
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_tblout(file_path):
    hits = set()
    with open(file_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            if len(fields) >= 1:
                hits.add(fields[0])
    return hits

def extract_seqs(ids, fasta_path, output_path):
    if not ids:
        print(f"  [!] No hit IDs provided; skipping {fasta_path}")
        return False

    if not os.path.exists(fasta_path):
        print(f"  [!] FASTA file not found: {fasta_path}")
        return False

    duplicate_log_path = os.path.join(os.path.dirname(output_path), "duplicate_ids.log")
    missing_log_path = os.path.join(os.path.dirname(output_path), "missing_sequences.log")

    id_counts = defaultdict(int)
    records = {}
    duplicates = []
    missing = []

    for rec in SeqIO.parse(fasta_path, "fasta"):
        id_counts[rec.id] += 1
        if rec.id not in records:
            records[rec.id] = rec
        else:
            duplicates.append(rec.id)

    if duplicates:
        with open(duplicate_log_path, 'a') as log:
            for dup_id in sorted(set(duplicates)):
                log.write(f"{dup_id}\t{fasta_path}\t{id_counts[dup_id]} instances\n")
        print(f"  [!] {len(set(duplicates))} duplicate ID(s) found in {fasta_path} — logged to {duplicate_log_path}")

    matched_records = []
    for i in ids:
        if i in records:
            matched_records.append(records[i])
        else:
            print(f"    [!] Warning: {i} not found in {fasta_path}")
            missing.append((i, fasta_path))

    if missing:
        with open(missing_log_path, 'a') as log:
            for missing_id, path in missing:
                log.write(f"{missing_id}\t{path}\n")
        print(f"  [!] Logged {len(missing)} missing hit(s) to {missing_log_path}")
    else:
        print(f"  [+] All {len(ids)} hit IDs were found in {fasta_path}")

    if matched_records:
        with open(output_path, 'w') as out:
            SeqIO.write(matched_records, out, "fasta")
        print(f"  [+] Wrote {len(matched_records)} sequences to {output_path}")
    else:
        print(f"  [!] No hits written from {fasta_path} (none matched); skipping file write")

    return bool(missing)

def main(hmmer_dir, root_dir=None):
    if root_dir is None:
        root_dir = os.path.abspath(os.path.join(hmmer_dir, os.pardir))

    print(f"\nStarting summary extraction from: {hmmer_dir}")
    print(f"Genome root directory: {root_dir}")
    print("--------------------------------------------------")

    summary = defaultdict(lambda: defaultdict(list))
    output_seq_dir = os.path.join(hmmer_dir, "sequences")
    os.makedirs(output_seq_dir, exist_ok=True)

    any_missing = False

    for file in os.listdir(hmmer_dir):
        if not file.endswith(".tbl"):
            continue

        tbl_path = os.path.join(hmmer_dir, file)
        base = file[:-4]  # remove .tbl extension

        # Require double-underscore delimited filenames
        if base.count("__") != 2:
            print(f"[!] Skipping malformed file (expected '__' delimiter x2): {file}")
            continue

        accession, fragment_id, hmm_name = base.split("__")

        print(f"\nProcessing: {file}")
        print(f"  Accession: {accession}")
        print(f"  Fragment:  {fragment_id}")
        print(f"  HMM:       {hmm_name}")

        hits = parse_tblout(tbl_path)
        print(f"  [+] Found {len(hits)} hit(s)")

        summary[(accession, fragment_id)][hmm_name] = hits

        faa_file = os.path.join(root_dir, accession, f"{accession}_{fragment_id}.faa")
        ffn_file = os.path.join(root_dir, accession, f"{accession}_{fragment_id}.ffn")

        hmm_subdir = os.path.join(output_seq_dir, hmm_name)
        os.makedirs(hmm_subdir, exist_ok=True)

        out_faa = os.path.join(hmm_subdir, f"{hmm_name}_{accession}_{fragment_id}.faa")
        out_ffn = os.path.join(hmm_subdir, f"{hmm_name}_{accession}_{fragment_id}.ffn")

        missing_faa = extract_seqs(hits, faa_file, out_faa)
        missing_ffn = extract_seqs(hits, ffn_file, out_ffn)
        if missing_faa or missing_ffn:
            any_missing = True

    hmm_names = sorted({hmm for d in summary.values() for hmm in d})

    # === Table 1: hmm_summary_counts.tsv (number of hits) ===
    counts_path = os.path.join(hmmer_dir, "hmm_summary_counts.tsv")
    with open(counts_path, 'w') as out:
        out.write("accession\tfragment_id\t" + "\t".join(hmm_names) + "\n")
        for (accession, fragment_id), hmm_hits in summary.items():
            row = [accession, fragment_id]
            for hmm in hmm_names:
                row.append(str(len(hmm_hits.get(hmm, []))))
            out.write("\t".join(row) + "\n")

    # === Table 2: hmm_summary.tsv (hit IDs) ===
    hits_path = os.path.join(hmmer_dir, "hmm_summary.tsv")
    with open(hits_path, 'w') as out:
        out.write("accession\tfragment_id\t" + "\t".join(hmm_names) + "\n")
        for (accession, fragment_id), hmm_hits in summary.items():
            row = [accession, fragment_id]
            for hmm in hmm_names:
                hits = hmm_hits.get(hmm, [])
                row.append(",".join(sorted(hits)) if hits else "")
            out.write("\t".join(row) + "\n")

    print("\n--------------------------------------------------")
    print(f"[✓] Summary written to: {hits_path}")
    print(f"[✓] Count summary written to: {counts_path}")
    print(f"[✓] Sequences saved in: {output_seq_dir}")
    if any_missing:
        print(f"[!] Some sequences were missing — see: {os.path.join(output_seq_dir, 'missing_sequences.log')}")
    else:
        print("[✓] All hit sequences found — no missing sequences.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize HMMER tblout files and extract matching sequences")
    parser.add_argument("--hmmer_dir", required=True, help="Path to HMMER output directory")
    parser.add_argument("--root_dir", required=False, help="Path to genome directories (default: parent of hmmer_dir)")
    args = parser.parse_args()
    main(args.hmmer_dir, args.root_dir)

