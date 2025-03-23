import pandas as pd
import numpy as np
import os
import sys

def summarize_splice_sites(profile_file, output_file):
    """Summarize splice site usage with log-transformed and normalized counts."""

    # Load input profile
    #df = pd.read_csv(profile_file, sep="\t")
    df = pd.read_csv(profile_file,sep="\t",header=None,names=["Chromosome", "Start", "End", "Strand", "Count"])

    if df.empty or "Count" not in df.columns:
        print("Input file is empty or missing 'Count' column.")
        open(output_file, 'w').close()
        return

    # Compute log2(count + 1)
    df["log_count"] = df["Count"].apply(lambda x: 0 if x == 0 else round(np.log2(x + 1), 2))

    # Normalize to per-million scale (like TPM/RPKM)
    total = df["Count"].sum()
    df["normalized_count"] = (df["Count"] / total) * 1e6 if total > 0 else 0

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write to output
    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python summarize_splice_sites.py <input_profile> <output_file>")
        sys.exit(1)

    summarize_splice_sites(sys.argv[1], sys.argv[2])
