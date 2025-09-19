#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import Counter

# Allowed amino acid symbols (uppercase + lowercase + ambiguity codes)
VALID_AA = set("XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy")

def clean_sequence(seq: str):
    """Remove invalid amino acid symbols from a sequence line.
    Returns cleaned sequence and a Counter of removed characters.
    """
    cleaned = []
    removed = Counter()
    for aa in seq:
        if aa in VALID_AA:
            cleaned.append(aa)
        else:
            removed[aa] += 1
    return "".join(cleaned), removed

def process_fastas(input_dir, output_dir_fastas, output_dir_summary, prefix="sample"):
    input_dir = Path(input_dir)
    output_dir_fastas = Path(output_dir_fastas)
    output_dir_summary = Path(output_dir_summary)
    output_dir_fastas.mkdir(parents=True, exist_ok=True)
    output_dir_summary.mkdir(parents=True, exist_ok=True)

    mapping_file = output_dir_summary / "name_mapping.tsv"
    cleanup_file = output_dir_summary / "cleanup_summary.tsv"
    species_file = output_dir_summary / "species_list.fa"

    cleanup_summary = {}
    renamed_species = []

    # Collect all files in input_dir
    fasta_files = sorted([f for f in input_dir.iterdir() if f.is_file()])
    total_files = len(fasta_files)

    with open(mapping_file, "w") as map_out:
        map_out.write("new_name\told_name\n")

        for idx, fasta_file in enumerate(fasta_files, start=1):
            print(f"{idx} of {total_files}: Processing {fasta_file.name}")

            new_name = f"{prefix}{idx:02d}.fasta"
            map_out.write(f"{new_name}\t{fasta_file.name}\n")
            renamed_species.append(f"{prefix}{idx:02d}")

            removed_chars = Counter()
            with open(fasta_file) as fin, open(output_dir_fastas / new_name, "w") as fout:
                seq_count = 0
                for line in fin:
                    line = line.strip()
                    if line.startswith(">"):
                        seq_count += 1
                        fout.write(f">{prefix}{idx:02d}_{seq_count}|\n")
                    elif line:
                        cleaned, removed = clean_sequence(line)
                        removed_chars.update(removed)
                        if cleaned:
                            fout.write(cleaned + "\n")

            cleanup_summary[fasta_file.name] = removed_chars
            print(f"   ✅ Done processing {fasta_file.name}")

    # Write species list
    with open(species_file, "w") as sp_out:
        for species_name in renamed_species:
            sp_out.write(f">{species_name}\n")

    # Console messages
    print(f"\n✅ Processed FASTA files saved to {output_dir_fastas}")
    print(f"✅ Mapping file saved to {mapping_file}")
    print(f"✅ Cleanup summary saved to {cleanup_file}")
    print(f"✅ Species list saved to {species_file}\n")

    print("📊 Cleanup summary (invalid characters removed):")
    for fname, removed in cleanup_summary.items():
        total = sum(removed.values())
        if total > 0:
            removed_str = ", ".join([f"{char}:{count}" for char, count in removed.items()])
            print(f"   {fname}: {total} removed ({removed_str})")
        else:
            print(f"   {fname}: 0 removed")

    # Write detailed cleanup summary
    with open(cleanup_file, "w") as cf:
        cf.write("file\tcharacter\tcount\n")
        for fname, removed in cleanup_summary.items():
            if removed:
                for char, count in removed.items():
                    cf.write(f"{fname}\t{repr(char)}\t{count}\n")
            else:
                cf.write(f"{fname}\t-\t0\n")

def main():
    parser = argparse.ArgumentParser(
        description="Standardize FASTA filenames and headers for phylogenomics projects, "
                    "generate species list, separate summary files, and show progress."
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing input FASTA files (only FASTAs should be present)"
    )
    parser.add_argument(
        "output_dir_fastas",
        help="Directory where cleaned/renamed FASTA files will be saved"
    )
    parser.add_argument(
        "output_dir_summary",
        help="Directory where summary files (mapping, cleanup, species list) will be saved"
    )
    parser.add_argument(
        "-p", "--prefix",
        default="sample",
        help="Prefix for renamed files and headers (default: 'sample')"
    )

    args = parser.parse_args()
    process_fastas(args.input_dir, args.output_dir_fastas, args.output_dir_summary, args.prefix)


if __name__ == "__main__":
    main()
