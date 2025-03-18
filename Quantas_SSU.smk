configfile: "config.yaml"

#Input a list of fastq files:
def get_samples(sample_file):
    with open(sample_file, "r") as f:
        return [line.strip().split('.')[0] for line in f if line.strip()]

rule all:
    input:
        "results/step1_mapping/.created",
        "results/step2_bedfiles/.created",
        "results/step3_splice_usage/.created",
        "results/step4_expression_matrix/.created",
        "results/mapping_stats/.created",
        expand("{output_dir_mapping}/{sample}.Aligned.out.sam", output_dir_mapping=config["output_dir_mapping"], sample=get_samples("samples.txt")),
        expand("{output_dir_bed}/{sample}.bed", output_dir_bed=config["output_dir_bed"], sample=get_samples("samples.txt")),
        expand("{output_dir_stats}/{sample}_mapping_stats.txt",output_dir_stats=config["output_dir_stats"],sample=get_samples("samples.txt")),
        expand("{output_dir_splice}/{sample}.splice_usage.txt", output_dir_splice=config["output_dir_splice"], sample=get_samples("samples.txt")),
        "results/expression_matrix.txt"

#Output directory
rule create_dirs:
    output:
        touch("results/step1_mapping/.created"),
        touch("results/step2_bedfiles/.created"),
        touch("results/step3_splice_usage/.created"),
        touch("results/step4_expression_matrix/.created"),
        touch("results/mapping_stats/.created")
    shell:
        """
        mkdir -p results/step1_mapping results/step2_bedfiles results/step3_splice_usage results/step4_expression_matrix
        """
# 1. STAR Mapping
rule mapping_star:
    input:
        # r1=config["fastq_dir"] + "/{sample}_R1.fastq.gz",
        # r2=config["fastq_dir"] + "/{sample}_R2.fastq.gz"
        r1=lambda wildcards: f"{config['fastq_dir']}/{wildcards.sample}_R1.fastq.gz", 
        r2=lambda wildcards: f"{config['fastq_dir']}/{wildcards.sample}_R2.fastq.gz" 
    output:
        "{output_dir_mapping}/{sample}.Aligned.out.sam" 
    # log:
    #     "logs/{sample}_star.log"
    threads: config["threads"]
    shell:         #mkdir -p {config[output_dir_mapping]} directory("results")
        """
        STAR --genomeDir {config[genome_dir]} --readFilesCommand gunzip -c --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {wildcards.output_dir_mapping}/{wildcards.sample}. --runThreadN {threads} --outFilterMismatchNmax {config[star_params][outFilterMismatchNmax]} --outFilterMismatchNoverReadLmax {config[star_params][outFilterMismatchNoverReadLmax]} --outFilterIntronMotifs {config[star_params][outFilterIntronMotifs]} --alignSJoverhangMin {config[star_params][alignSJoverhangMin]} --alignSJDBoverhangMin {config[star_params][alignSJDBoverhangMin]} --chimSegmentMin {config[star_params][chimSegmentMin]} > "logs/{wildcards.sample}_star.log" 2>&1
        """

# 2. Convert SAM to BED
rule sam_to_bed:
    input:
        lambda wildcards: f"{config['output_dir_mapping']}/{wildcards.sample}.Aligned.out.sam" #"{output_dir_mapping}/{sample}.Aligned.out.sam"
    output:
        "{output_dir_bed}/{sample}.bed"
    # log:
    #     "logs/{sample}_sam2bed.log"
    threads: config["threads_sam2bed"]
    shell:
        """
        python scripts/sam2bed_speed.py {input} {output} \
        {config[sam2bed_params][primary_only]} \
        {config[sam2bed_params][uniq_only]} \
        {config[sam2bed_params][reverse_strand]} \
        {config[sam2bed_params][use_rna_strand]} \
        > "logs/{wildcards.sample}_sam2bed.log" 2>&1
        """
#> {log} 2>&1

# 3. Reads Statistics
rule mapping_stats:
    input:
        sam=lambda wildcards: f"{config['output_dir_mapping']}/{wildcards.sample}.Aligned.out.sam", #"{output_dir_mapping}/{sample}.Aligned.out.sam",
        bed=lambda wildcards: f"{config['output_dir_bed']}/{wildcards.sample}.bed" #"{output_dir_bed}/{sample}.bed"
    output:
        "{output_dir_stats}/{sample}_mapping_stats.txt"
    # log:
    #     "logs/mapping_stats.log"
    threads: config["threads"]
    shell:
        """
        python scripts/mapping_stats.py {input.sam} {input.bed} {output} > "logs/{wildcards.sample}_read_stats.log" 2>&1
        """

# 4. Summarize Splice Site Usage
# rule summarize_splice_usage:
#     input:
#         lambda wildcards: f"{config['output_dir_bed']}/{wildcards.sample}.bed" #"{output_dir_bed}/{sample}.bed"
#     output:
#         "{output_dir_splice}/{sample}.splice_usage.txt"
#     # log:
#     #     "logs/{sample}_splice_usage.log"
#     threads: config["threads"]
#     shell:
#         """
#         python scripts/summarize_splice_sites.py {input} {output} \
#         {config[ssu_params][min_anchor]} {config[ssu_params][big_file]} \
#         {config[ssu_params][weight]} {config[ssu_params][separate_strand]} \
#         > "logs/{wildcards.sample}_sumarize_ssu.log" 2>&1
#         """

rule filter_bed:
    input:
        "data/tag.bed"
    output:
        "processed/tag_filtered.bed"
    script:
        "scripts/bed_parser.py"

rule count_splice_sites:
    input:
        intron="data/intron.bed",
        tag="processed/tag_filtered.bed"
    output:
        "results/splice_site_counts.txt"
    script:
        "scripts/splice_site_usage.py"

rule summarize_usage:
    input:
        "results/splice_site_counts.txt"
    output:
        "results/summary.txt"
    script:
        "scripts/summarize_splice_sites.py"


