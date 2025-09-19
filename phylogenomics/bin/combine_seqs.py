#!/usr/bin/env python3
import argparse
from collections import defaultdict, Counter
from pathlib import Path

def read_species(expected_file):
    """Read expected species identifiers from file."""
    species = {}
    with open(expected_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]
            if line:
                species[line] = ""  # initialize with empty string
    return species

def concatenate_sequences(expected, ortholog_files):
    """Concatenate sequences per species across ortholog files."""
    for ortholog_file in ortholog_files:
        with open(ortholog_file) as f:
            current_species = None
            current_sequence = []

            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    # save previous
                    if current_species and current_species in expected:
                        expected[current_species] += "".join(current_sequence)

                    # header like >species_number|
                    header = line[1:]
                    current_species = header.split("_", 1)[0]
                    current_sequence = []
                else:
                    if current_species:
                        current_sequence.append(line)

            # save last
            if current_species and current_species in expected:
                expected[current_species] += "".join(current_sequence)

    return expected

def write_output(expected, output_file, line_length=0):
    """Write concatenated sequences to FASTA file."""
    with open(output_file, "w") as out:
        for species in sorted(expected.keys()):
            seq = expected[species]
            out.write(f">{species}\n")
            if line_length > 0:
                for i in range(0, len(seq), line_length):
                    out.write(seq[i:i+line_length] + "\n")
            else:
                out.write(seq + "\n")

def check_lengths(expected):
    """Check sequence length consistency across species."""
    lengths = {sp: len(seq) for sp, seq in expected.items()}
    counts = Counter(lengths.values())

    if len(counts) == 1:
        length = next(iter(counts))
        print(f"✅ All species have the same concatenated length: {length} bp")
    else:
        print("❌ Length inconsistencies detected:")
        for sp, length in sorted(lengths.items()):
            print(f"   {sp}: {length} bp")

        # majority length
        expected_len = counts.most_common(1)[0][0]
        print(f"\nExpected (majority) length: {expected_len} bp")
        print("Problematic species:")
        for sp, length in sorted(lengths.items()):
            if length != expected_len:
                print(f"   {sp} ({length} bp)")

def main():
    parser = argparse.ArgumentParser(
        description="Concatenate ortholog FASTA files across species and validate sequence lengths."
    )
    parser.add_argument("--expected_file", required=True,
                        help="File with expected species list (lines like >species)")
    parser.add_argument("--output_file", required=True,
                        help="Path to write concatenated FASTA")
    parser.add_argument("--ortholog_files", nargs="+", required=True,
                        help="Ortholog FASTA files to concatenate")
    parser.add_argument("--line_length", type=int, default=0,
                        help="Wrap FASTA sequences at this line length (0 = no wrapping)")
    
    args = parser.parse_args()

    expected = read_species(args.expected_file)
    expected = concatenate_sequences(expected, args.ortholog_files)
    write_output(expected, args.output_file, args.line_length)
    check_lengths(expected)

    print(f"\nConcatenated FASTA file created: {args.output_file}")

if __name__ == "__main__":
    main()
