import sys
import pysam
import pandas as pd

def trim(s):
    """去除字符串前后空格"""
    return s.strip()

def count_total_reads(sam_file):
    """计算总 reads 数"""
    with pysam.AlignmentFile(sam_file, "r") as sam:
        return sum(1 for _ in sam)

def count_mapped_reads(sam_file):
    """计算 mapped reads 数（排除 FLAG 4）"""
    with pysam.AlignmentFile(sam_file, "r") as sam:
        return sum(1 for read in sam if not read.is_unmapped)

def count_unique_mapped_reads(bed_file):
    """计算唯一映射 reads 数（从 BED 文件行数）"""
    df = pd.read_csv(bed_file, sep="\t", header=None)
    return df.shape[0]

def count_unique_junction_reads(bed_file):
    df = pd.read_csv(bed_file, sep="\t", header=None)
    if df.shape[1] > 10:  # Ensure column 10 exists
        return df[df[10].astype(str).str.contains(",")].shape[0]
    else:
        return 0  # No junction reads if column doesn't exist

def compute_mapping_stats(sam_file, bed_file, outfile):
    try:
        total_reads = count_total_reads(sam_file)
        mapped_reads = count_mapped_reads(sam_file)
        unique_mapped_reads = count_unique_mapped_reads(bed_file)
        unique_junction_reads = count_unique_junction_reads(bed_file)

        mapped_percent = f"{(mapped_reads / total_reads) * 100:.3f}%"
        unique_mapped_percent = f"{(unique_mapped_reads / total_reads) * 100:.3f}%"
        unique_junction_percent = f"{(unique_junction_reads / total_reads) * 100:.3f}%"

        result = f"{sam_file}\t{total_reads}\t{mapped_reads}\t{mapped_percent}\t{unique_mapped_reads}\t{unique_mapped_percent}\t{unique_junction_reads}\t{unique_junction_percent}\n"
        print(result)

        with open(outfile, "w") as out:
            out.write(result)

    except FileNotFoundError as e:
        print(f"文件缺失: {e}")

if __name__ == "__main__":
    sam_file = sys.argv[1]
    bed_file = sys.argv[2]
    outfile = sys.argv[3]
    compute_mapping_stats(sam_file, bed_file, outfile)


