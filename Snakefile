import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"

SAMPLES, = glob_wildcards("samples/raw/200421/{sample}_R1.fastq.gz")

md = pd.read_table(config["omic_meta_data"], index_col="SampleID",dtype=str)
condition = config["linear_model"]
baseline = config["TE_baseline"]

# TO FILTER 
control = md[md[condition]==baseline]
treat = md.loc[~md.index.isin(control.index)] 

control_paths = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in control.index]

treat['cond_paths'] = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in treat.index]
cond_paths = {key:x['cond_paths'].values.tolist() for key,x in treat.groupby(condition)}
CONDITIONS = cond_paths.keys()

# Wildcard function to grab proper condition
def get_TE(wildcards):
    return cond_paths[wildcards.condition]

ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
read_dist_ext = ['txt']

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        "results/multiqc/{project_id}_multiqc.html".format(project_id=config["project_id"]),
        "data/{project_id}_counts.txt".format(project_id=config["project_id"]),
        expand("samples/fastqscreen/200421/{sample}/{sample}_{ext}_screen.{fastqscreen_ext}", sample=SAMPLES, ext=ext, fastqscreen_ext=fastqscreen_ext),
        "results/tables/{project_id}_read_coverage.txt".format(project_id=config["project_id"]),
        "results/tables/{project_id}_tag_counts.txt".format(project_id=config["project_id"]),
        #expand("samples/miRNA_count/{sample}_featureCount_miRNA_count.txt", sample = SAMPLES),
        #expand("samples/star_miRNA/{sample}_bam/Aligned.out.sam", sample = SAMPLES),
        #expand("samples/star_miRNA_mature/{sample}_bam/Aligned.out.sam", sample = SAMPLES),
        #expand("samples/star_miRNA_hairpin/{sample}_bam/Aligned.out.sam", sample = SAMPLES),
        #expand("samples/hsa_miRNA_count/{sample}_featureCount_miRNA_count.txt", sample = SAMPLES),
        #expand("results/tables/{project_id}_miRNA_STAR_mapping_statistics.txt", project_id = config['project_id']),
        #expand("results/tables/{project_id}_mature_miRNA_STAR_mapping_statistics.txt", project_id = config['project_id']),
        #expand("results/tables/{project_id}_hairpin_miRNA_STAR_mapping_statistics.txt", project_id = config['project_id'])
        expand("samples/star_TE/{sample}/Aligned.out.bam", sample = SAMPLES),
        expand("results/TEtranscripts/{condition}.cntTable", condition = CONDITIONS)



include: "rules/align_rmdp.smk"
include: "rules/miRNA.smk"


