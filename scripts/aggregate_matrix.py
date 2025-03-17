import sys
import pandas as pd
import json

def aggregate_samples(samples_json, output_file):
    """合并多个样本的基因表达数据"""
    with open(samples_json, "r") as f:
        groups = json.load(f)

    all_data = []
    sample_names = []

    for group, samples in groups.items():
        for sample in samples:
            df = pd.read_csv(f"processed/{sample}.txt", sep="\t")
            df = df[["gene_id", "gene_name", "RPKM"]]
            df.columns = ["gene_id", "gene_name", sample]
            all_data.append(df)
            sample_names.append(sample)

    merged_df = all_data[0]
    for df in all_data[1:]:
        merged_df = pd.merge(merged_df, df, on=["gene_id", "gene_name"], how="outer")

    merged_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    aggregate_samples(sys.argv[1], sys.argv[2])
