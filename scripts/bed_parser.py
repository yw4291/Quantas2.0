import pandas as pd

def read_bed(file_path):
    """读取 BED 文件，返回 DataFrame"""
    cols = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    df = pd.read_csv(file_path, sep="\t", comment="#", names=cols, header=None)
    return df

def filter_exons(input_bed, output_bed):
    """筛选至少有一个外显子的 BED 行"""
    df = read_bed(input_bed)
    df = df[df["chromEnd"] - df["chromStart"] > 1]  # 仅保留长度 > 1 的片段
    df.to_csv(output_bed, sep="\t", index=False, header=False)

if __name__ == "__main__":
    import sys
    filter_exons(sys.argv[1], sys.argv[2])
