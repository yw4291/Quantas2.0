import pandas as pd
import sys

def count_splice_sites(intron_bed, tag_bed, profile_output):
    """计算剪接位点的使用情况，并生成 profile"""
    
    # 读取 intron 和 tag BED 文件
    introns = pd.read_csv(intron_bed, sep="\t", names=["chrom", "start", "end", "name", "score", "strand"])
    tags = pd.read_csv(tag_bed, sep="\t", names=["chrom", "start", "end", "name", "score", "strand"])
    
    # 计算剪接位点的标签数
    merged = pd.merge(introns, tags, on=["chrom", "strand"], suffixes=("_intron", "_tag"))
    merged = merged[(merged["start_tag"] >= merged["start_intron"]) & (merged["end_tag"] <= merged["end_intron"])]

    # 计算 splice site 计数
    site_counts = merged.groupby(["chrom", "start_intron", "end_intron", "strand"]).size().reset_index(name="count")

    # **生成 profile 数据**
    profile = site_counts.groupby(["chrom", "strand"]).apply(lambda x: x.sort_values(by="start_intron")).reset_index(drop=True)

    # 输出
    site_counts.to_csv(profile_output, sep="\t", index=False)

if __name__ == "__main__":
    count_splice_sites(sys.argv[1], sys.argv[2], sys.argv[3])
