import pysam
import sys
import pandas as pd

def sam_to_bed(input_sam, output_bed, primary_only=False, uniq_only=False, reverse_strand=0, use_rna_strand=False):
    """
    解析 SAM 文件，并转换为 BED 格式
    """
    samfile = pysam.AlignmentFile(input_sam, "r")
    bed_records = []

    for read in samfile.fetch():
        if read.is_unmapped:
            continue

        # 解析 Read Flags
        is_primary = not read.is_secondary
        is_unique = "XT:A:U" in read.to_string() if read.has_tag("XT") else False

        if primary_only and not is_primary:
            continue
        if uniq_only and not is_unique:
            continue

        # 解析 BED 信息
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

        bed_records.append([chrom, start, end, name, score, strand])

    # 转换为 DataFrame 并保存
    bed_df = pd.DataFrame(bed_records, columns=["chrom", "start", "end", "name", "score", "strand"])
    bed_df.to_csv(output_bed, sep="\t", index=False, header=False)

if __name__ == "__main__":
    sam_to_bed(sys.argv[1], sys.argv[2], primary_only=bool(int(sys.argv[3])), uniq_only=bool(int(sys.argv[4])), reverse_strand=int(sys.argv[5]), use_rna_strand=bool(int(sys.argv[6])))
