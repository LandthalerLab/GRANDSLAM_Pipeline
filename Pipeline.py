import pandas as pd
import os

# dataset and samples
metadata = pd.read_table(config["metadata"], dtype={"sample": str}, sep = " ").set_index(["sample"], drop=False)
sample_set = metadata["sample"].tolist()

# target rule
rule all:
	input:
		SLAM_out = config['result_dir'] + "/BulkSLAMSeq/SLAM.tsv.gz",
		alignment_summary = expand(config['result_dir'] + "STAR/{sample}Log.final.out", sample = sample_set),
		contamination_summary = expand(config['result_dir'] + "{sample}_contamination_summary.out", sample = sample_set) if config["remove_contamination"] else [],
		dedup_summary = config['result_dir'] + "all_dedup_summary.out" if config["remove_duplicates"] else [] 

### 
# Duplicate Removal with cd-hit-dup (Optional) ####
###

rule duplicateRemoval:
	input:
	  read1 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read1"],
	  read2 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read2"]
	params:
		qsub = config['qsub_duplicateRemoval'],
		result_dir = config['result_dir']
	output:
		read1_nodup= config['result_dir'] + "fastq_files/{sample}_1_nodup.fastq",
		read2_nodup= config['result_dir'] + "fastq_files/{sample}_2_nodup.fastq",
		dedup_summary = config['result_dir'] + "{sample}_dedup_summary.out"
	shell:
		'''
		#!/bin/bash
	
		# Dedup search
		mkdir -p {params.result_dir}/fastq_files
		zcat {input.read1} > {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq
		zcat {input.read2} > {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq
		cd-hit-dup -i {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq -i2 {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq -o {output.read1_nodup} -o2 {output.read2_nodup} 
	
		# Dedup summary
		mkdir -p {params.result_dir}
		total_reads=$(wc -l < {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq)
		dedup_reads=$(wc -l < {output.read1_nodup})
		echo $((total_reads/4)) $((dedup_reads/4)) >> {output.dedup_summary}

		# remove temporary fastq files
		rm {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq 
		rm {params.result_dir}/fastq_files/{wildcards.sample}_*.clstr
		'''
  
rule cd_hit_dup_summary:
	input:
		dedup_summary = expand(config['result_dir'] + "{sample}_dedup_summary.out", sample = sample_set) if config["remove_duplicates"] else [],
	params:
	  qsub = "",
		result_dir = config['result_dir'],
		scripts_path = config["scripts_path"]
	output:
		config['result_dir'] + "all_dedup_summary.out"
	script: 
		"{params.scripts_path}/dedup_summary.py"
	
### 
# Remove Contamination (Optional) ####
###

# select fastq files
def get_fastq_files_afterduplicateremoval(sample, readn):
        if config["remove_duplicates"]:
                return config['result_dir'] + "fastq_files/" + sample + "_" + readn + "_nodup.fastq"
        else:
                return config["fastq_dir"] + metadata.loc[sample, "read" + readn]

rule removeContamination:
	input:
		read1 = lambda wildcards: get_fastq_files_afterduplicateremoval(wildcards.sample, "1"),
		read2 = lambda wildcards: get_fastq_files_afterduplicateremoval(wildcards.sample, "2")
	params:
		qsub = config['qsub_removeContamination'],
		result_dir = config['result_dir'],
		contamination_reference = config['contamination_reference'],
		bowtie2_param = config["bowtie2_param"]
	output:
                read1_norrna= config['result_dir'] + "fastq_files/{sample}_1_norrna.fastq",
                read2_norrna= config['result_dir'] + "fastq_files/{sample}_2_norrna.fastq",
		contReads_summary = config['result_dir'] + "{sample}_contamination_summary.out"
	shell:
		'''
		mkdir -p {params.result_dir}/fastq_files
		/gnu/var/guix/profiles/custom/landthaler/bin/bowtie2 {params.bowtie2_param} --un-conc {params.result_dir}/fastq_files/{wildcards.sample}_%_norrna.fastq --al-conc {params.result_dir}/fastq_files/{wildcards.sample}_%_rrna.fastq -x {params.contamination_reference} -1 {input.read1} -2 {input.read2} 2>> {output.contReads_summary}
		'''
		
### 
# Align Reads ####
###

# select fastq files
def get_fastq_files_foralignment(sample, readn):
	if config["remove_contamination"]:
		return config['result_dir'] + "fastq_files/" + sample + "_" + readn + "_norrna.fastq"
	else:
		return get_fastq_files_afterduplicateremoval(sample, readn)

rule run_STAR:
	input:
		read1 = lambda wildcards: get_fastq_files_foralignment(wildcards.sample, "1"),
                read2 = lambda wildcards: get_fastq_files_foralignment(wildcards.sample, "2")
	params:
		qsub = config['qsub_STAR'],
		result_dir = config['result_dir'],
		reference = config['reference'], 
		STAR_param = config['STAR_param']
	output:
		sortedreads = config['result_dir'] + "bam_files/{sample}_sorted.bam", 
		alignment_summary = config['result_dir'] + "STAR/{sample}Log.final.out"
	shell:
		'''
		mkdir -p {params.result_dir}/bam_files
	  # STAR --genomeDir {params.reference} --readFilesIn {input.read1} {input.read2} {params.STAR_param}  --outFileNamePrefix {params.result_dir}/STAR/{wildcards.sample} --outStd SAM > {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam	
	  STAR --genomeDir {params.reference} --readFilesIn {input.read1} {input.read2} {params.STAR_param}  --outFileNamePrefix {params.result_dir}/STAR/{wildcards.sample} --outStd SAM > {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam	
		samtools sort {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam -o {output.sortedreads}
		samtools index {output.sortedreads}
		'''

### 
# GRAND-SLAM ####
###

rule bam2cit:
	input:
		bam_files = expand(config['result_dir'] + "bam_files/{sample}_sorted.bam", sample = sample_set)
	params: 
	  result_dir = config['result_dir'],
		qsub = config['qsub_bam2cit']
	output: 
		cit_file = config['result_dir'] + "cit_files/BulkSLAMSeq.cit"
	shell:
		'''
		# make cit list
		rm -f BulkSLAMSeq
		IFS=', ' read -r -a array <<< "{input.bam_files}"
		for bamfile in "${{array[@]}}"
		do
    			echo $bamfile >> BulkSLAMSeq
		done	
		
		# make cit file
		mkdir -p {params.result_dir}/cit_files
		bamlist2cit -p BulkSLAMSeq
		mv BulkSLAMSeq.cit {params.result_dir}/cit_files/.
		mv BulkSLAMSeq.cit.metadata.json {params.result_dir}/cit_files/.
		mv BulkSLAMSeq {params.result_dir}/cit_files/.
		'''

rule run_GRANDSLAM:
	input: 
		cit_file = config['result_dir'] + "cit_files/BulkSLAMSeq.cit"
	params:
		qsub = config['qsub_GRANDSLAM'],
		result_dir = config['result_dir'],
		grandslam_param = config['grandslam_param'],
		GRANDSLAM_reference = config['GRANDSLAM_reference']
	output:
		SLAM_out = config['result_dir'] + "/BulkSLAMSeq/SLAM.tsv.gz"
	shell:
		'''
		set +eu; eval "$(conda shell.bash hook)"; conda activate r-base-4.0;
		gedi -e Slam {params.grandslam_param} -genomic {params.GRANDSLAM_reference} -prefix {params.result_dir}/BulkSLAMSeq/SLAM -reads {input.cit_file}
		'''
