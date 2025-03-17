import sys
import pandas as pd
import numpy as np

def process_expression(input_file, output_file, pseudo_count=1, log2_transform=False):
    """处理单个样本的基因表达数据"""
    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ["gene_id", "gene_name", "tagNum", "RPKM"]

    if log2_transform:
        df["RPKM"] = np.log2(df["RPKM"] + pseudo_count)

    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    process_expression(sys.argv[1], sys.argv[2], float(sys.argv[3]), bool(int(sys.argv[4])))
