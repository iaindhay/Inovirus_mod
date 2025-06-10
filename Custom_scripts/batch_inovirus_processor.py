import argparse
import os
import sys
from extract_inovirus_data_individual import main

def setup_arg_parser():
    parser = argparse.ArgumentParser(description='Batch process multiple *_refined.csv files for inovirus extraction.')
    parser.add_argument('-dir', '--directory', required=True, help='Root directory to recursively process.')
    parser.add_argument('--do_fna', action='store_true', help='Only generate .fna region output')
    parser.add_argument('--do_gff', action='store_true', help='Only generate GFF output')
    parser.add_argument('--do_faa', action='store_true', help='Only generate .faa protein output')
    parser.add_argument('--do_ffn', action='store_true', help='Only generate .ffn CDS output')
    return parser

def process_directory(root_dir, do_fna, do_gff, do_faa, do_ffn):
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith("refined.csv"):
                csv_path = os.path.join(dirpath, filename)
                try:
                    print("\n============================================")
                    print(f"Processing: {csv_path}")
                    main(csv_file_path=csv_path,
                         fasta_file_path=None,
                         gff_file_path=None,
                         protein_file_path=None,
                         cds_file_path=None,
                         do_fna=do_fna,
                         do_gff=do_gff,
                         do_faa=do_faa,
                         do_ffn=do_ffn)
                    print("Done.")
                except Exception as e:
                    print(f"Error processing {csv_path}: {e}")

if __name__ == "__main__":
    args = setup_arg_parser().parse_args()
    if not os.path.isdir(args.directory):
        print("Error: --directory must be a valid folder path")
        sys.exit(1)
    process_directory(args.directory, args.do_fna, args.do_gff, args.do_faa, args.do_ffn)
