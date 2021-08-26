rule fastqc:
    input:
        fwd = "samples/raw/200421/{sample}_R1.fastq.gz",
        rev = "samples/raw/200421/{sample}_R2.fastq.gz"
    output:
        "samples/fastqc/200421/{sample}/{sample}_R1_fastqc.zip",
        "samples/fastqc/200421/{sample}/{sample}_R2_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """fastqc --outdir  samples/fastqc/200421/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""


rule fastqscreen:
    input:
        fwd = "samples/raw/200421/{sample}_R1.fastq.gz",
        rev = "samples/raw/200421/{sample}_R2.fastq.gz"
    output:
        "samples/fastqscreen/200421/{sample}/{sample}_R1_screen.html",
        "samples/fastqscreen/200421/{sample}/{sample}_R1_screen.png",
        "samples/fastqscreen/200421/{sample}/{sample}_R1_screen.txt",
        "samples/fastqscreen/200421/{sample}/{sample}_R2_screen.html",
        "samples/fastqscreen/200421/{sample}/{sample}_R2_screen.png",
        "samples/fastqscreen/200421/{sample}/{sample}_R2_screen.txt",

    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/200421/{wildcards.sample} {input.fwd} {input.rev}"""

rule multiqc:
    input:
        expand("samples/fastqc/200421/{sample}/{sample}_{ext}_fastqc.zip", sample = SAMPLES, ext = ext)
    output:
        "results/multiqc/{project_id}_multiqc.html".format(project_id=config["project_id"])
    params:
        project_id=config["project_id"]
    conda:
        "../envs/multiqc.yaml"
    shell:
        """multiqc --outdir results/multiqc --filename {params.project_id}_multiqc samples/fastqc/200421"""
    

rule STAR:
    input:
        fwd = "samples/raw/200421/{sample}_R1.fastq.gz",
        rev = "samples/raw/200421/{sample}_R2.fastq.gz"
    output:
        "samples/star/200421/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/200421/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/200421/{sample}_bam/Log.final.out"
    params:
        gtf=config["gtf_file"]
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN 24 --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input} \
                --outFileNamePrefix samples/star/200421/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --twopassMode Basic
                """)


rule star_statistics:
    input:
        expand("samples/star/200421/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"


rule htseq:
    input:
        "samples/star/200421/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/htseq/200421_new/{sample}_htseq_gene_count.txt"
    params:
        gtf = config["gtf_file"]
    conda:
        "../envs/htseq.yaml"
    shell:
        """htseq-count \
                -f bam \
                -r pos \
                -s reverse \
                -m intersection-strict \
                {input} \
                {params.gtf} > {output}"""
    
rule compile_counts:
    input:
        expand("samples/htseq/200421_new/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table.py"

rule read_distribution:
    input:
        "samples/star/200421/{sample}_bam/Aligned.sortedByCoord.out.bam"
    params:
        bed=config['bed_file']
    output:
        "rseqc/read_distribution/200421/{sample}/{sample}.read_distribution.txt",
    conda:
        "../envs/rseqc.yaml"
    shell:
       "read_distribution.py -i {input} -r {params.bed} > {output}"

rule compile_rd:
    input:
        expand("rseqc/read_distribution/200421/{sample}/{sample}.read_distribution.txt", sample=SAMPLES)
    output:
        "results/tables/{project_id}_read_coverage.txt".format(project_id=config["project_id"])
    script:
        "../scripts/get_rd.py"

rule compile_rd_tagcounts:
    input:
        expand("rseqc/read_distribution/200421/{sample}/{sample}.read_distribution.txt", sample=SAMPLES)
    output:
        "results/tables/{project_id}_tag_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/get_rd_tags.py"

rule map_TE:
    input:
        fwd = "samples/raw/200421/{sample}_R1.fastq.gz",
        rev = "samples/raw/200421/{sample}_R2.fastq.gz"
    output:
        "samples/star_TE/{sample}/Aligned.out.bam",
        "samples/star_TE/{sample}/Log.final.out"
    params:
        index=config["star_index"],
        gtf=config["gtf_file"]
    run:
        STAR=config["star_tool"]

        shell("""
                {STAR} --runThreadN 12 --genomeDir {params.index} --sjdbGTFfile {params.gtf} \
                       --sjdbOverhang 100 --readFilesIn {input.fwd} {input.rev} \
                       --outFileNamePrefix samples/star_TE/{wildcards.sample}/  --readFilesCommand zcat \
                       --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100""")

rule TEtranscripts:
    input:
        control = control_paths,
        treat = get_TE
    output:
        "results/TEtranscripts/{condition}.cntTable",
        "results/TEtranscripts/{condition}_sigdiff_gene_TE.txt"
    conda:
        "../envs/TE.yaml"
    params:
        gtf=config["gtf_file"],
        TE_gtf=config["TE_gtf"]
    shell:
        """TEtranscripts --format BAM --stranded reverse -t {input.treat} -c {input.control} \
                         --minread 1 -i 10 --padj 0.05 --GTF {params.gtf} --TE {params.TE_gtf} \
                         --mode multi --project results/TEtranscripts/{wildcards.condition}"""
