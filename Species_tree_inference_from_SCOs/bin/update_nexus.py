#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

def load_models(models_path: Path):
    """Load models into a dict: {SCO_name: Model}"""
    models = {}
    with open(models_path, "r") as f:
        header = True
        for line in f:
            if header:
                header = False
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                models[parts[0]] = parts[1]
    return models

def main():
    parser = argparse.ArgumentParser(
        description="Insert best-fit models into a NEXUS file with charpartition definitions.",
        epilog="Example:\n"
               "  ./makenex.py -i alignments_concat.nex.nex "
               "-m best_models_fixed.tsv "
               "-o alignments_with_models_AIC.nex",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-nexus",
        required=True,
        help="Input NEXUS alignment file (with charset definitions)."
    )
    parser.add_argument(
        "-m", "--models",
        required=True,
        help="TSV file containing SCO and best-fit models."
    )
    parser.add_argument(
        "-o", "--output-nexus",
        required=True,
        help="Output NEXUS file with models inserted."
    )

    args = parser.parse_args()

    nexus_path = Path(args.input_nexus)
    models_path = Path(args.models)
    out_path = Path(args.output_nexus)

    if not nexus_path.is_file():
        print(f"❌ Error: NEXUS file not found: {nexus_path}")
        sys.exit(1)
    if not models_path.is_file():
        print(f"❌ Error: Models file not found: {models_path}")
        sys.exit(1)

    # Load models
    models = load_models(models_path)
    print(f"✅ Loaded {len(models)} models from {models_path}")

    charset_lines = []
    charpartition_parts = []

    # Read NEXUS and collect charsets
    with open(nexus_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.lower().startswith("charset"):
                charset_lines.append(stripped)

                # Example: "charset SCOgreat90.txt.OrthoGroup0.fasta.out.trim = 1-994;"
                nexus_gene_name = stripped.split()[1]

                model = models.get(nexus_gene_name)
                if model:
                    charpartition_parts.append(f"{model}:{nexus_gene_name}")
                else:
                    print(f"⚠️  Warning: No model found for {nexus_gene_name}")

    if not charset_lines:
        print("❌ Error: No charset lines found in NEXUS file.")
        sys.exit(1)

    # Write new NEXUS with models
    with open(out_path, "w") as out:
        out.write("#nexus\n")
        out.write("begin sets;\n")
        for line in charset_lines:
            out.write(f"  {line}\n")
        out.write("\n")
        out.write(f"  charpartition mymodels = {', '.join(charpartition_parts)};\n")
        out.write("end;\n")

    print(f"✅ New NEXUS with models written to {out_path}")

if __name__ == "__main__":
    main()
