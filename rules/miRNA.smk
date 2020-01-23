rule bowtie2:
    input:
        fwd = "samples/raw/191101/{sample}_R1.fastq.gz",
        rev = "samples/raw/191101/{sample}_R2.fastq.gz"
    output:
        "samples/miR_bam/{sample}.bam"
    params:
        index = "/home/groups/CEDAR/anno/indices/bowtie2/miRbase/miRbase"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """bowtie2 -x {params.index} -1 {input.fwd} -2 {input.rev} | samtools view -bS - > {output}"""

rule sort_bowtie2:
    input:
        "samples/miR_bam/{sample}.bam"
    output:
        "samples/miR_bam/{sample}_sort.bam"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """samtools sort -O bam -n {input} -o {output}"""

rule count_miRNA:
    input:
        "samples/miR_bam/{sample}_sort.bam"
    output:
        "samples/miRNA_count/{sample}_featureCount_miRNA_count.txt"
    params:
        gff = "/home/groups/CEDAR/anno/gtf/miRbase_hsa.gff3"
    shell:
        """/home/exacloud/lustre1/CEDAR/tools/subread-1.6.2-Linux-x86_64/bin/featureCounts -t miRNA -g Name -O -s 1 -M -a {params.gff} -o {output} {input}"""

rule star_miRNA:
    input:
        fwd = "samples/raw/191101/{sample}_R1.fastq.gz",
        rev = "samples/raw/191101/{sample}_R2.fastq.gz"
    output:
        "samples/star_miRNA/{sample}_bam/Aligned.out.sam",
        "samples/star_miRNA/{sample}_bam/Log.final.out"
    run:
        STAR=config["star_tool"],
        pathToGenomeIndex = "/home/groups/CEDAR/anno/indices/star/mirna/"

        shell("""{STAR} --runThreadN 12 --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} --genomeSAindexNbases 6 \
                --outFileNamePrefix samples/star_miRNA/{wildcards.sample}_bam/ \
                --readFilesCommand zcat --outFilterMultimapNmax 20 \
                --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
                --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 \
                --alignIntronMax 1 
                """)

rule star_miRNA_mature:
    input:
        fwd = "samples/raw/191101/{sample}_R1.fastq.gz",
        rev = "samples/raw/191101/{sample}_R2.fastq.gz"
    output:
        "samples/star_miRNA_mature/{sample}_bam/Aligned.out.sam",
        "samples/star_miRNA_mature/{sample}_bam/Log.final.out"
    run:
        STAR=config["star_tool"],
        pathToGenomeIndex = "/home/groups/CEDAR/anno/indices/star/miRbase_mature/"

        shell("""{STAR} --runThreadN 12 --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} --genomeSAindexNbases 6 \
                --outFileNamePrefix samples/star_miRNA_mature/{wildcards.sample}_bam/ \
                --readFilesCommand zcat --outFilterMultimapNmax 20 \
                --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
                --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 \
                --alignIntronMax 1 
                """)

rule star_miRNA_hairpin:
    input:
        fwd = "samples/raw/191101/{sample}_R1.fastq.gz",
        rev = "samples/raw/191101/{sample}_R2.fastq.gz"
    output:
        "samples/star_miRNA_hairpin/{sample}_bam/Aligned.out.sam",
        "samples/star_miRNA_hairpin/{sample}_bam/Log.final.out"
    run:
        STAR=config["star_tool"],
        pathToGenomeIndex = "/home/groups/CEDAR/anno/indices/star/miRbase_hairpin/"

        shell("""{STAR} --runThreadN 12 --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} --genomeSAindexNbases 6 \
                --outFileNamePrefix samples/star_miRNA_hairpin/{wildcards.sample}_bam/ \
                --readFilesCommand zcat --outFilterMultimapNmax 20 \
                --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
                --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 \
                --alignIntronMax 1 
                """)

rule filter_star_mirna:
    input:
        "samples/star_miRNA/{sample}_bam/Aligned.out.sam"
    output:
        "samples/star_miRNA/{sample}_bam/Aligned.filtered.sam"
    shell:
        """awk '{S=0; split($6,C,/[0-9]*/); n=split($6,L,/[NMSID]/);  if (and($2,0x10)>0 && C[n]=="S") {S=L[n-1]} else if (and($2,0x10)==0 && C[2]=="S") {S=L[1]}; if (S<=1) print }' {input} > {output}"""

rule sort_star:
    input:
        "samples/star_miRNA/{sample}_bam/Aligned.out.sam"
    output:
        "samples/star_miRNA/{sample}_bam/Aligned.filtered.sorted.bam"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """samtools sort -O bam -n {input} -o {output}"""

rule count_star_mirna:
    input:
        "samples/star_miRNA/{sample}_bam/Aligned.filtered.sorted.bam"
    output:
        "samples/hsa_miRNA_count/{sample}_featureCount_miRNA_count.txt"
    params:
        gff = "/home/groups/CEDAR/anno/gtf/miRbase_hsa.gff3"
    shell:
        """/home/exacloud/lustre1/CEDAR/tools/subread-1.6.2-Linux-x86_64/bin/featureCounts -t miRNA -g Name -O -s 1 -M -a {params.gff} -o {output} {input}"""

rule star_statistics_mirna:
    input:
        expand("samples/star_miRNA/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_miRNA_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"

rule star_statistics_mature:
    input:
        expand("samples/star_miRNA_mature/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_mature_miRNA_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"

rule star_statistics_hairpin:
    input:
        expand("samples/star_miRNA_hairpin/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_hairpin_miRNA_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"
