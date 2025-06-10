import os
import sys
import csv

if len(sys.argv) < 2:
    print("Usage: python3 zot_summary_extractor_gffstyle.py /path/to/parent_directory")
    sys.exit(1)

parent_dir = sys.argv[1]
summary_data = []
no_zot_files = []
no_zot_accessions = set()
target_suffix = "-refined.csv"

print(f"[INFO] Searching in: {parent_dir}")

for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.endswith(target_suffix):
            path = os.path.join(root, file)
            trimmed_name = os.path.basename(file).split("_annot")[0]
            loci, pfams, pcs = [], [], []
            print(f"[FILE] Checking: {path}")
            try:
                with open(path, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith("#") or line.startswith(">"):
                            continue  # skip comments and headers
                        fields = line.strip().split(",")
                        if len(fields) < 8:
                            continue
                        metadata = fields[7]
                        annotations = dict(
                            entry.strip().split("=", 1)
                            for entry in metadata.strip().split(";")
                            if "=" in entry
                        )
                        pfam_val = annotations.get("pfam", "").strip()
                        if pfam_val.lower().startswith("zot"):
                            print(f"  [HIT] Found pfam=Zot in: {metadata}")
                            locus = annotations.get("locus_tag", "NA").split()[0]
                            pc = annotations.get("inoPC", "NA").split()[0]
                            loci.append(locus)
                            pfams.append("pfam=Zot")
                            pcs.append(pc)
                if loci:
                    summary_data.append([
                        trimmed_name,
                        ";".join(loci),
                        ";".join(pfams),
                        ";".join(pcs)
                    ])
                else:
                    no_zot_files.append(trimmed_name)
                    accession = trimmed_name.split("_frag")[0]
                    no_zot_accessions.add(accession)
                    print("  [NO ZOT] No zot entries in this file.")
            except Exception as e:
                print(f"[ERROR] Failed to process {path}: {e}")

# Write matched summary
if summary_data:
    print(f"[RESULT] Writing {len(summary_data)} entries to zot_summary.tsv")
    with open("zot_summary.tsv", "w", newline='') as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["file", "locus_tag", "pfam", "inoPC"])
        writer.writerows(summary_data)
else:
    print("[RESULT] No zot hits found in any files.")

# Write no-zot file list
if no_zot_files:
    with open("no_zot_found.txt", "w") as nf:
        for name in no_zot_files:
            nf.write(name + "\n")
    print(f"[INFO] Logged {len(no_zot_files)} files with no Zot to no_zot_found.txt")

# Write unique accessions without Zot
if no_zot_accessions:
    with open("accessions_without_zot.txt", "w") as af:
        for acc in sorted(no_zot_accessions):
            af.write(acc + "\n")
    print(f"[INFO] Logged {len(no_zot_accessions)} accessions without Zot to accessions_without_zot.txt")
