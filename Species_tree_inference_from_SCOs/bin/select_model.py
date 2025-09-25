#!/usr/bin/env python3
import re
import csv
import argparse
from pathlib import Path

def extract_model_from_log(log_file: Path, criterion: str):
    """Extract the best-fit model from an IQ-TREE log file for a given criterion."""
    with open(log_file, "r", encoding="utf-8") as f:
        for line in f:
            # Typical line: Best-fit model: JTT+G4 chosen according to AIC
            m = re.search(
                r"Best-fit model:\s+([^\s]+)\s+chosen according to\s+" + re.escape(criterion),
                line
            )
            if m:
                return m.group(1).strip()
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Extract best-fit models from IQ-TREE logs and fix SCO names."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Root directory containing all gene tree folders.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file.")
    parser.add_argument("-crit", "--criterion",
                        choices=["AIC", "BIC", "AICc"],
                        default="AIC",
                        help="Criterion for model selection (default: AIC)")
    parser.add_argument("-s", "--suffix",
                        default=".fasta",
                        help="Suffix to append to each SCO name (default: .fasta)")

    args = parser.parse_args()

    root_dir = Path(args.input)
    output_tsv = Path(args.output)
    criterion = args.criterion
    suffix = args.suffix

    results = []

    # Loop over each SCO directory
    for sco_dir in root_dir.iterdir():
        if sco_dir.is_dir() and sco_dir.name.startswith("SCO"):
            log_files = list(sco_dir.glob("*.log"))
            if not log_files:
                print(f"⚠️  No log file found in {sco_dir}")
                continue

            log_file = log_files[0]  # should only be one
            model = extract_model_from_log(log_file, criterion)
            if model:
                sco_name = sco_dir.name + suffix
                results.append((sco_name, model))
            else:
                print(f"⚠️  No {criterion} model found in {log_file}")

    if not results:
        print("❌ No models extracted.")
        return

    # Write TSV
    with open(output_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["SCO", f"{criterion}_Model"])
        writer.writerows(results)

    print(f"✅ Wrote {len(results)} results to {output_tsv}")

if __name__ == "__main__":
    main()
