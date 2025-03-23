import pandas as pd
import pyranges as pr
import sys

def count_splice_sites(intron_bed, tag_bed, profile_output):
    # Only load required columns as strings for performance
    col_names = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    dtype = {"Chromosome": "category", "Start": "int32", "End": "int32", "Strand": "category"}

    introns = pd.read_csv(intron_bed, sep="\t", names=col_names, usecols=range(6), dtype=dtype)
    tags = pd.read_csv(tag_bed, sep="\t", names=col_names, usecols=range(6), dtype=dtype)

    # Convert to PyRanges (very fast interval engine)
    intron_gr = pr.PyRanges(introns)
    tag_gr = pr.PyRanges(tags)

    # Fast interval-based join (chrom + strand-aware)
    overlap = intron_gr.join(tag_gr)

    # Group by intron coordinates and strand, count tags
    grouped = (
        overlap.df
        .groupby(["Chromosome", "Start", "End", "Strand"], sort=False)
        .size()
        .reset_index(name="Count")
    )

    # Write output (no header index)
    grouped.to_csv(profile_output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    count_splice_sites(sys.argv[1], sys.argv[2], sys.argv[3])
