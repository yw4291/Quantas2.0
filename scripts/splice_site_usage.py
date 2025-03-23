import pandas as pd
import pyranges as pr
import os
import sys

def count_splice_sites(intron_bed, tag_bed, profile_output):
    # Define BED column names and assign memory-efficient dtypes
    col_names = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    dtype = {"Chromosome": "category", "Start": "int32", "End": "int32", "Strand": "category"}

    # Read intron and tag files with only the first 6 BED columns (BED6)
    introns = pd.read_csv(intron_bed, sep="\t", names=col_names, usecols=range(6), dtype=dtype)
    tags = pd.read_csv(tag_bed, sep="\t", names=col_names, usecols=range(6), dtype=dtype)

    # Optional: filter out overly long reads that could be noise/artifacts
    tags = tags[(tags["End"] - tags["Start"]) <= 1000]

    # Convert both datasets into PyRanges objects for efficient interval handling
    intron_gr = pr.PyRanges(introns)
    tag_gr = pr.PyRanges(tags)

    # Perform interval containment join (tags fully inside introns)
    overlap = intron_gr.join(tag_gr, how="containment", nb_cpu=1)

    # Extract relevant columns, drop malformed rows
    df = overlap.df[["Chromosome", "Start", "End", "Strand"]].dropna()

    # Count number of reads supporting each intron (no cartesian explosion)
    grouped = (
        df.value_counts(sort=False)
          .reset_index(name="Count")
    )

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(profile_output), exist_ok=True)

    # Always write output file (even if it's empty)
    if grouped.empty:
        open(profile_output, 'w').close()
    else:
        grouped.to_csv(profile_output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    # Expect: python splice_site_usage.py intron.bed tags.bed output.txt
    count_splice_sites(sys.argv[1], sys.argv[2], sys.argv[3])
