import pysam
import sys

def sam_to_bed(input_sam, output_bed, primary_only=False, uniq_only=False, reverse_strand=0, use_rna_strand=False):
    """
    Optimized version of SAM to BED conversion for better speed and memory efficiency.
    """
    samfile = pysam.AlignmentFile(input_sam, "r")

    with open(output_bed, "w") as bed_out:
        for read in samfile.fetch(until_eof=True):  # Fetches reads more efficiently
            if read.is_unmapped:
                continue

            # Filter based on flags
            if primary_only and read.is_secondary:
                continue
            if uniq_only and (not read.has_tag("XT") or read.get_tag("XT") != "U"):
                continue

            # Extract fields
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            name = read.query_name
            score = read.mapping_quality
            strand = "-" if read.is_reverse else "+"

            if reverse_strand:
                strand = "+" if strand == "-" else "-"

            if use_rna_strand and read.has_tag("XS"):
                strand = read.get_tag("XS")

            # Write directly to file (avoids large memory consumption)
            bed_out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

if __name__ == "__main__":
    sam_to_bed(
        sys.argv[1],
        sys.argv[2],
        primary_only=bool(int(sys.argv[3])),
        uniq_only=bool(int(sys.argv[4])),
        reverse_strand=int(sys.argv[5]),
        use_rna_strand=bool(int(sys.argv[6]))
    )
