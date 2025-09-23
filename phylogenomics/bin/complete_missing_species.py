#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_species_list(species_file):
    """Read expected species identifiers from file, stripping '>' if present."""
    species = []
    with open(species_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]
            if line:
                species.append(line)
    return species

def get_sequences(fasta_file):
    """Read a FASTA file into dict {header: sequence}."""
    seqs = {}
    current_header = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    seqs[current_header] = "".join(current_seq)
                current_header = line[1:]  # keep full header
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            seqs[current_header] = "".join(current_seq)
    return seqs

def max_seq_length(seqs):
    """Return maximum sequence length in a dict of sequences."""
    return max((len(s) for s in seqs.values()), default=0)

def write_fasta(output_file, seqs, missing, max_length, line_length=60):
    """Write sequences + padded missing species into FASTA format."""
    with open(output_file, "w") as out:
        # Write existing sequences
        for header, seq in seqs.items():
            out.write(f">{header}\n")
            for i in range(0, len(seq), line_length):
                out.write(seq[i:i+line_length] + "\n")

        # Write missing species with padded sequences
        for species in missing:
            out.write(f">{species}_00|\n")
            padding = "-" * max_length
            for i in range(0, len(padding), line_length):
                out.write(padding[i:i+line_length] + "\n")

def validate_output(output_dir, output_ext, expected_species_count):
    """Check that each output file has the right species count and equal sequence lengths."""
    output_dir = Path(output_dir)
    issues = []

    for f in sorted(output_dir.glob(f"*{output_ext}")):
        seqs = get_sequences(f)
        num_species = len(seqs)
        lengths = {len(s) for s in seqs.values()}

        if num_species != expected_species_count:
            issues.append(f"{f.name}: has {num_species}, expected {expected_species_count}")

        if len(lengths) != 1:
            issues.append(f"{f.name}: inconsistent sequence lengths {sorted(lengths)}")

    if issues:
        print("\n❌ Validation issues found:")
        for issue in issues:
            print("   " + issue)
    else:
        print("\n✅ All output files passed validation: correct species count and equal lengths.")

def process_files(species_file, input_dir, input_ext, output_dir, output_ext, line_length, missing_report=None):
    species_list = read_species_list(species_file)
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # open report file if requested
    report_handle = None
    if missing_report:
        report_handle = open(missing_report, "w")
        report_handle.write("orthogroup\tmissing_species\n")

    for orthogroup in sorted(input_dir.glob(f"*{input_ext}")):
        print(f"Processing {orthogroup.name}")

        # Get sequences from orthogroup file
        seqs = get_sequences(orthogroup)

        # Extract species names (strip after first "_")
        existing_species = {h.split("_", 1)[0] for h in seqs}

        # Find missing species
        missing_species = [sp for sp in species_list if sp not in existing_species]

        if missing_species:
            for sp in missing_species:
                print(f"⚠️  Missing species in {orthogroup.name}: {sp}")
            if report_handle:
                report_handle.write(f"{orthogroup.stem}\t{','.join(missing_species)}\n")

        # Determine max sequence length
        max_length = max_seq_length(seqs)

        # Output path (replace extension with output_ext)
        output_file = output_dir / (orthogroup.stem + output_ext)

        # Write updated FASTA
        write_fasta(output_file, seqs, missing_species, max_length, line_length)

    if report_handle:
        report_handle.close()
        print(f"✅ Missing species report saved to {missing_report}")

    # Final validation step
    validate_output(output_dir, output_ext, len(species_list))

def main():
    parser = argparse.ArgumentParser(
        description="Fill missing species in orthogroup FASTA files with padded sequences.",
        formatter_class=argparse.RawTextHelpFormatter,
        usage=(
            "complete_missing_species.py\n"
            "  --species_file SPECIES_FILE\n"
            "  --input_dir INPUT_DIR\n"
            "  [--input_ext INPUT_EXT]\n"
            "  --output_dir OUTPUT_DIR\n"
            "  [--output_ext OUTPUT_EXT]\n"
            "  [--line_length LINE_LENGTH]\n"
            "  [--missing_report MISSING_REPORT]"
        )
    )

    parser.add_argument("--species_file", required=True,
                        help="File with expected species list (lines like >species)")
    parser.add_argument("--input_dir", required=True,
                        help="Directory with cleaned orthogroup files")
    parser.add_argument("--input_ext", default=".fa",
                        help="Input file extension to select (default: .fa)")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write completed orthogroup files")
    parser.add_argument("--output_ext", default=".fasta",
                        help="Output file extension (default: .fasta)")
    parser.add_argument("--line_length", type=int, default=60,
                        help="Max characters per sequence line (default: 60)")
    parser.add_argument("--missing_report",
                        help="Optional path to write a table of missing species per orthogroup")

    args = parser.parse_args()

    process_files(
        args.species_file,
        args.input_dir,
        args.input_ext,
        args.output_dir,
        args.output_ext,
        args.line_length,
        args.missing_report,
    )

if __name__ == "__main__":
    main()
