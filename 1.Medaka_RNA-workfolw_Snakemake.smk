#!/usr/bin/env snakemake

samples = ["10dph_XX1", "10dph_XX2", "10dph_XX3", "10dph_XY1", "10dph_XY2", "10dph_XY3", "10dph_XYD1", "10dph_XYD2", "10dph_XYD3", "30dph_XX1", "30dph_XX2", "30dph_XX3", "30dph_XY1", "30dph_XY2", "30dph_XY3", "30dph_XYD1", "30dph_XYD2", "30dph_XYD3"，"120dph_XX1", "120dph_XX2", "120dph_XX3", "120dph_XY1", "120dph_XY2", "120dph_XY3", "120dph_XYD1", "120dph_XYD2", "120dph_XYD3"]
threads = 16

rule all:
    input:
       # expand("results_ola/prepare/fastqc/{sample}_{read}_fastqc.html", sample=samples, read=[1, 2]),
       # expand("results_ola/prepare/trimmed/{sample}_{read}.fastq.gz", sample=samples, read=["R1", "R2"]),
        "results_ola/mapping/star_index",
        expand("results_ola/mapping/mapped/{sample}.Aligned.sortedByCoord.out.bam", sample=samples),
        expand("results_ola/mapping/filtered/{sample}.bam", sample=samples),
        expand("results_ola/mapping/uniq/{sample}.bam", sample=samples),
        expand("results_ola/expression/featureCounts/{sample}.txt", sample=samples),

# FastQC

#rule fastqc:
 #   input:
  #      "data/datasets_ola/{name}.fq.gz"
  #  output:
   #     html = "results_ola/prepare/fastqc/{name}_fastqc.html",
    #    data = "results_ola/prepare/fastqc/{name}_fastqc.zip"
   # log:
    #    "results_ola/prepare/fastqc/{name}.log"
   # params:
    #    outdir = "results_ola/prepare/fastqc"
   # shell:
    #    """
     #   fastqc -o {params.outdir} {input} &> {log}
      #  """

# Trimmed

#rule cutadapt:
 #   input:
  #      r1 = "data/datasets_ola/{sample}_1.fq.gz",
   #     r2 = "data/datasets_ola/{sample}_2.fq.gz"
   # output:
    #    r1 = "results_ola/prepare/trimmed/{sample}_R1.fastq.gz",
     #   r2 = "results_ola/prepare/trimmed/{sample}_R2.fastq.gz",
    #log:
     #   "results_ola/prepare/trimmed/{sample}.log"
   # threads:
    #    threads
   # shell:
    #    """
     #   cutadapt -j {threads} -a "AGATCGGAAGAGCACACGT" -A "AGATCGGAAGAGCGTCGTG" \
      #      -m 30 -q 30 --max-n 2 -o {output.r1} -p {output.r2} {input} &> {log}
      #  """

# Mapping 

rule star_index:
    input:
        fasta = "data/genome/oryzias_latipes/GCF_002234675.1_ASM223467v1_genomic.modified.fasta",
        gff = "data/genome/oryzias_latipes/GCF_002234675.1_ASM223467v1_genomic.modified.sorted1.gff"
    output:
        directory("results_ola/mapping/star_index")
    log:
        "results_ola/mapping/star_index.log"
    threads:
        threads
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --sjdbGTFfile {input.gff} \
            --sjdbGTFtagExonParentTranscript Parent \
            --genomeFastaFiles {input.fasta} &> {log}
        """

rule star_mapping:
    input:
        r1 = "data/datasets/{sample}_1.fq.gz",
        r2 = "data/datasets/{sample}_2.fq.gz",
        idx = rules.star_index.output
    output:
        bam = "results_ola/mapping/mapped/{sample}.Aligned.sortedByCoord.out.bam"
    log:
        "results_ola/mapping/mapped/{sample}.log"
    params:
        prefix = "results_ola/mapping/mapped/{sample}"
    threads:
        threads
    shell:
        """
        STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd \
            --outFilterMultimapNmax 10 --readFilesCommand zcat --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. --genomeDir {input.idx} \
            --readFilesIn {input.r1} {input.r2} \
	    --quantMode TranscriptomeSAM GeneCounts &> {log}
        """


rule bamFilter:
    input:
        bam = rules.star_mapping.output.bam
    output:
        bam = "results_ola/mapping/filtered/{sample}.bam",
    log:
        "results_ola/mapping/filtered/{sample}.log"
    shell:
        """
        bamtools filter -in {input.bam} -out {output.bam} -tag NH:1 -isProperPair true  &> {log}
        """

rule remove_duplicate:
    input:
        bam = rules.bamFilter.output.bam,
        bai = rules.bamFilter.output.bam + ".bai"
    output:
        bam = "results_ola/mapping/uniq/{sample}.bam",
        txt = "results_ola/mapping/uniq/{sample}_metrics.txt"
    log:
        "results_ola/mapping/uniq/{sample}.log"
    shell:
        """(
        set +u
        source activate picard
        picard MarkDuplicates REMOVE_DUPLICATES=true I={input.bam} O={output.bam} M={output.txt}
        conda deactivate ) &> {log}
        """

# Expression

rule featureCounts:
    input:
        bam = rules.remove_duplicate.output.bam,
        gtf = "data/genome/oryzias_latipes/GCF_002234675.1_ASM223467v1_genomic.modified.sorted1.gff"
    output:
        "results_ola/expression/featureCounts/{sample}.txt"
    log:
        "results_ola/expression/featureCounts/{sample}.log"
    threads:
        threads
    shell:
        """
	featureCounts -T {threads} -p -C -B -t exon -g gene -a {input.gtf} -o {output} {input.bam} >{log} 2>&1
        """
#	featureCounts --fracOverlap 0.5 -s 0 -p -T {threads} \
#		            -a {input.gtf} -o {output} {input.bam} &> {log}

# Common rules

rule bamIndex:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        bamtools index -in {input}
        """

rule bamStats:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """


