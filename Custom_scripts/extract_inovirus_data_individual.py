import argparse
import csv
import os
import sys
import re
from Bio import SeqIO
import gffutils

# -------- File Finder Functions -------- #
"""
def find_file(csv_file_path, file_extension):
    base_name = '_'.join(os.path.basename(csv_file_path).split('_')[:2])
    directory = os.path.dirname(csv_file_path) or os.getcwd()
    for file in os.listdir(directory):
        if file.startswith(base_name) and file.endswith(file_extension):
            return os.path.join(directory, file)
    return None
"""
def find_nucleotide_file(csv_file_path):
    base_name = os.path.basename(csv_file_path).split('_frag')[0]
    directory = os.path.dirname(csv_file_path) or os.getcwd()
    for file in os.listdir(directory):
        if file.startswith(base_name) and (file.endswith('_genomic.fna') or file.endswith('.fasta')):
            return os.path.join(directory, file)
    return None

def find_gff_file(csv_file_path):
    base_name, _ = os.path.basename(csv_file_path).split('_frag', 1)
    directory = os.path.dirname(csv_file_path) or os.getcwd()
    expected_gff_filename = f"{base_name}.gff"
    if expected_gff_filename in os.listdir(directory):
        return os.path.join(directory, expected_gff_filename)
    return None

def find_protein_file(csv_file_path):
    base_name = os.path.basename(csv_file_path).split('_annot')[0]
    directory = os.path.dirname(csv_file_path) or os.getcwd()
    for file in os.listdir(directory):
        if file == f"{base_name}_prots.faa":
            return os.path.join(directory, file)
    return None

def find_cds_file(csv_file_path):
    base_name = '_'.join(os.path.basename(csv_file_path).split('_')[:2])
    directory = os.path.dirname(csv_file_path) or os.getcwd()
    for file in os.listdir(directory):
        if file.startswith(base_name) and file.endswith('_cds.fna'):
            return os.path.join(directory, file)
    return None

# -------- Core Extraction Functions -------- #

def extract_prophage_details_from_csv(csv_file):
    contig_name = None
    fragment_id = None
    start = stop = None
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader, start=1):
            if i <= 2 or not row:
                continue
            if row[0].startswith('>'):
                match = re.match(r"(.+?)_(\d+)-(\d+)", row[1])
                if match:
                    fragment_id, start, stop = match.groups()
            else:
                if contig_name is None:
                    contig_name = row[0]
                elif contig_name != row[0]:
                    sys.exit(f"Contig mismatch: {row[0]} vs {contig_name}")
            if contig_name and fragment_id and start and stop:
                break
    return contig_name, fragment_id, int(start), int(stop)

def write_sequence_region_to_fasta(fna_file, contig_name, start, stop, csv_file_path):
    output_fasta_file = csv_file_path.rsplit('.', 1)[0] + '.fna'
    for record in SeqIO.parse(fna_file, "fasta"):
        if record.id.startswith(contig_name):
            region = record.seq[start - 1:stop]
            from Bio.SeqRecord import SeqRecord
            region_record = SeqRecord(region, id=f"{contig_name}:{start}-{stop}", description="")
            SeqIO.write(region_record, output_fasta_file, "fasta")
            print(f"Wrote: {output_fasta_file}")
            return
    print(f"Contig {contig_name} not found in {fna_file}")

def process_gff_and_create_output(gff_file_path, contig_name, fragment_id, start, stop, csv_file_path):
    output_gff_file = csv_file_path.rsplit('.', 1)[0] + '.gff'
    db = gffutils.create_db(gff_file_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge')
    with open(output_gff_file, 'w') as out:
        for line in open(gff_file_path):
            if line.startswith('#'):
                out.write(line)
        out.write(f"##extracted-{fragment_id}-inovirus-region {contig_name}:{start}-{stop}\n###\n")
        for feature in db.all_features():
            if re.match(f"^{re.escape(contig_name)}(\\.\\d+)?$", feature.seqid):
                if feature.start >= start and feature.end <= stop:
                    out.write(str(feature) + '\n')
    print(f"Wrote: {output_gff_file}")

def extract_and_write_protein_sequences(gff_file_path, protein_file_path, csv_file_path):
    annotations = {}
    with open(gff_file_path, 'r') as gff:
        for line in gff:
            if not line.startswith('#') and line.strip():
                parts = line.strip().split('\t')
                if parts[2] == "CDS":
                    attrs = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
                    locus_tag = attrs.get('locus_tag', 'unknown')
                    name = attrs.get('Name', attrs.get('Note', 'unknown'))
                    annotations[locus_tag] = name
    output_faa = csv_file_path.rsplit('.', 1)[0] + '.faa'
    with open(protein_file_path) as faa, open(output_faa, 'w') as out:
        for record in SeqIO.parse(faa, "fasta"):
            locus_tag = record.id
            if locus_tag in annotations:
                name = annotations[locus_tag]
                record.description = ''
                record.id = f"{name}\t{locus_tag}"
                SeqIO.write(record, out, "fasta")
    print(f"Wrote: {output_faa}")

def extract_and_write_gene_sequences_from_cds_fasta(cds_fasta_path, protein_file_path, csv_file_path):
    valid_locus_tags = set()
    for record in SeqIO.parse(protein_file_path, "fasta"):
        valid_locus_tags.add(record.id.split('\t')[-1])
    output_gene_fasta = csv_file_path.rsplit('.', 1)[0] + '.ffn'
    gene_records = []
    for record in SeqIO.parse(cds_fasta_path, "fasta"):
        desc = record.description
        locus_tag = re.search(r'\[locus_tag=(.*?)\]', desc)
        protein_id = re.search(r'\[protein_id=(.*?)\]', desc)
        if locus_tag and protein_id:
            tag = locus_tag.group(1)
            if tag in valid_locus_tags:
                record.id = f"{protein_id.group(1)}\t{tag}"
                record.description = ''
                gene_records.append(record)
    SeqIO.write(gene_records, output_gene_fasta, "fasta")
    print(f"Wrote: {output_gene_fasta}")

# -------- CLI Handling -------- #

def check_files_exist(*file_paths):
    missing = [f for f in file_paths if not os.path.exists(f)]
    if missing:
        print("Missing required files:")
        for f in missing:
            print(" -", f)
        return False
    return True

def setup_arg_parser():
    parser = argparse.ArgumentParser(description='Extract inovirus region data.')
    parser.add_argument('-csv', '--csv_file', required=True, help='Path to the CSV file.')
    parser.add_argument('-nuc', '--nucleotide')
    parser.add_argument('-ann', '--annotation')
    parser.add_argument('-prot', '--protein')
    parser.add_argument('-cds', '--cds_fasta')
    parser.add_argument('--do_fna', action='store_true')
    parser.add_argument('--do_gff', action='store_true')
    parser.add_argument('--do_faa', action='store_true')
    parser.add_argument('--do_ffn', action='store_true')
    return parser

def main(csv_file_path, fasta_file_path=None, gff_file_path=None, protein_file_path=None, cds_file_path=None,
         do_fna=False, do_gff=False, do_faa=False, do_ffn=False):

    if not any([do_fna, do_gff, do_faa, do_ffn]):
        do_fna = do_gff = do_faa = do_ffn = True

    fasta_file_path = fasta_file_path or find_nucleotide_file(csv_file_path)
    gff_file_path = gff_file_path or find_gff_file(csv_file_path)
    protein_file_path = protein_file_path or find_protein_file(csv_file_path)
    cds_file_path = cds_file_path or find_cds_file(csv_file_path)

    print(f"\nProcessing: {csv_file_path}")
    print(f"CSV: {csv_file_path}")
    print(f"FASTA: {fasta_file_path}")
    print(f"GFF: {gff_file_path}")
    print(f"Protein: {protein_file_path}")
    print(f"CDS: {cds_file_path}")
    if not check_files_exist(csv_file_path, fasta_file_path, gff_file_path, protein_file_path):
        print("Aborting due to missing files.")
        return  # <-- allow batch script to continue


    contig_name, fragment_id, start, stop = extract_prophage_details_from_csv(csv_file_path)

    if do_gff:
        process_gff_and_create_output(gff_file_path, contig_name, fragment_id, start, stop, csv_file_path)
    if do_fna:
        write_sequence_region_to_fasta(fasta_file_path, contig_name, start, stop, csv_file_path)
    if do_faa:
        extract_and_write_protein_sequences(gff_file_path, protein_file_path, csv_file_path)
    if do_ffn and cds_file_path and os.path.exists(cds_file_path):
        extract_and_write_gene_sequences_from_cds_fasta(cds_file_path, protein_file_path, csv_file_path)

if __name__ == "__main__":
    args = setup_arg_parser().parse_args()
    main(args.csv_file, args.nucleotide, args.annotation, args.protein, args.cds_fasta,
         do_fna=args.do_fna, do_gff=args.do_gff, do_faa=args.do_faa, do_ffn=args.do_ffn)
