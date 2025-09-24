#!/usr/bin/env python3
import argparse
from pathlib import Path
import math

def count_fastas(fasta_dir, extension=".faa"):
    fasta_dir = Path(fasta_dir)
    return len(list(fasta_dir.glob(f"*{extension}")))

def process_proteinortho(tsv_file, fasta_dir, extension=".faa", output_dir="."):
    num_fastas = count_fastas(fasta_dir, extension)
    thresholds = list(range(100, 49, -5))  # 100, 95, 90, ... 50

    print(f"✅ Found {num_fastas} input FASTA files")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read header lines once
    with open(tsv_file) as fin:
        header_lines = [line for line in fin if line.startswith("#")]

    for thresh in thresholds:
        cutoff = math.ceil(num_fastas * (thresh / 100.0))
        out_file = output_dir / (f"SCO_all.txt" if thresh == 100 else f"SCO_great{thresh}.txt")
        count = 0

        with open(out_file, "w") as fout:
            # Write header first
            for header in header_lines:
                fout.write(header)

            # Process data lines
            with open(tsv_file) as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue  # already written header
                    parts = line.strip().split("\t")
                    try:
                        n_all = int(parts[0])
                        n_have = int(parts[1])
                    except ValueError:
                        continue  # skip malformed lines

                    if thresh == 100:
                        if n_all == num_fastas and n_have == num_fastas:
                            fout.write(line)
                            count += 1
                    else:
                        if n_all == n_have and n_all >= cutoff:
                            fout.write(line)
                            count += 1

        # Print summary for this threshold
        if thresh == 100:
            print(f"{thresh}% ({num_fastas} species): {count} SCOs written to {out_file}")
        else:
            print(f"{thresh}% (≥{cutoff} species): {count} SCOs written to {out_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Process ProteinOrtho output to extract SCOs at different percentage thresholds."
    )
    parser.add_argument("tsv_file", help="ProteinOrtho .tsv file")
    parser.add_argument("fasta_dir", help="Directory with the modified input FASTA files")
    parser.add_argument("-e", "--extension", default=".faa", help="FASTA file extension (default: .faa)")
    parser.add_argument("-o", "--output_dir", default="proteinortho_scos", help="Directory for output SCO files")

    args = parser.parse_args()
    process_proteinortho(args.tsv_file, args.fasta_dir, args.extension, args.output_dir)

if __name__ == "__main__":
    main()
