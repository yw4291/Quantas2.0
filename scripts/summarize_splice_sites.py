import pandas as pd
import sys

def summarize_splice_sites(profile_file, output_file):
    """汇总剪接位点使用情况，整合 tag profile"""
    df = pd.read_csv(profile_file, sep="\t")

    # 计算 log2 变换
    df["log_count"] = df["count"].apply(lambda x: 0 if x == 0 else round(pd.np.log2(x + 1), 2))

    # 计算标准化值（RPKM 替代方案）
    df["normalized_count"] = df["count"] / df["count"].sum() * 1e6

    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    summarize_splice_sites(sys.argv[1], sys.argv[2])
