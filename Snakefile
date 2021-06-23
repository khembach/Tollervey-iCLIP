## Configuration file
import os
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("".join(["Make sure there is a config.yaml file in ", os.getcwd(), 
			" or specify one with the --configfile commandline parameter."]))

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'STARindex', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_STAR']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str

config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['STARindex'] = sanitizefile(config['STARindex'])
config['metatxt'] = sanitizefile(config['metatxt'])

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("".join(["Metadata file ", config["metatxt"], " does not exist."]))

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
	sys.exit("".join(["Make sure 'names' and 'type' are columns in ", config["metatxt"]]))


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

outputdir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html")

rule setup:
	input:
		os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		os.path.join(outputdir, "Rout", "softwareversions.done")

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	  	os.path.join(outputdir, "Rout", "pkginstall_state.txt")
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"],
		Rbin = Rbin
	priority:
		50
	conda:
		Renv
	log:
		os.path.join(outputdir, "Rout", "install_pkgs.Rout")
	benchmark:
	  	os.path.join(outputdir, "benchmarks", "install_pkgs.txt")
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

rule rundemultiplex:
	input:
		expand(os.path.join(outputdir, "demultiplexed", "{sample}_demux.fastq.gz"), sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext1"]), "_val_1_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", str(config["fqext2"]), "_val_2_fastqc.zip"])), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC", "{sample}_umi_removed_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()),
		expand(os.path.join(outputdir, "STARbigwig", "{sample}_Aligned.sortedByCoord.out.bw"), sample = samples.names.values.tolist())

rule runxlsites:
	input:
		expand(os.path.join(outputdir, "xlsites", "{sample}_cDNA_unique.bed"), sample = samples.names.values.tolist())

rule runannotate:
	input:
		expand(os.path.join(outputdir, "annotate", "{sample}_biotype.tab"), sample = samples.names.values.tolist())

rule run_removedup:
	input:
		expand(os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.bam.bai"), sample = samples.names.values.tolist())

rule run_filter_second_read:
	input:
		expand(os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.bam.bai"), sample = samples.names.values.tolist())

rule run_norm_bigwig:
	input:
		expand(os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.ucsc.replicatesNorm.bw"), sample = samples.names.values.tolist())

rule run_bigwig_second_read:
	input:
		expand(os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.bw"), sample = samples.names.values.tolist())


## List all the packages that were used by the R analyses
rule listpackages:
	log:
		os.path.join(outputdir, "Rout", "list_packages.Rout")
	params:
		Routdir = os.path.join(outputdir, "Rout"),
		outtxt = os.path.join(outputdir, "R_package_versions.txt"),
		script = "scripts/list_packages.R",
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(os.path.join(outputdir, "Rout", "softwareversions.done"))
	log:
		os.path.join(outputdir, "logs", "softversions.log")
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'ARMOR version ' && cat version; "
		"salmon --version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; STAR --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		os.path.join(config["STARindex"], "SA"),
		os.path.join(config["STARindex"], "chrNameLength.txt")
	log:
		os.path.join(outputdir, "logs", "STAR_index.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_index.txt")
	params:
		STARindex = lambda wildcards, output: os.path.dirname(output[0]),   ## dirname of first output
		readlength = config["readlength"],
		starextraparams = config["additional_star_index"]
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength} "
		"{params.starextraparams}"


## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

# demultiplexed reads
rule fastqcdemultiplexed:
	input:
		fastq = os.path.join(outputdir, "demultiplexed", "{sample}.fastq.gz")
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_demultiplexed_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_demultiplexed_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"



## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}.fq.gz")
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_trimmed_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_trimmed_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"



# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "FastQC", "{sample}_demux_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
	if config["run_trimming"]:
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed", "{sample}_umi_removed_trimmed.fq.gz"), sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC", "{sample}_umi_removed_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
	if config["run_STAR"]:
		input.extend(expand(os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [os.path.join(outputdir, "FastQC")]
	if config["run_trimming"]:
		param.append(os.path.join(outputdir, "FASTQtrimmed"))
	if config["run_STAR"]:
		param.append(os.path.join(outputdir, "STAR"))
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html")
	params:
		inputdirs = multiqc_params,
		MultiQCdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "multiqc.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "multiqc.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Demultiplex and extract random barcode
## ------------------------------------------------------------------------------------ ##
# the iCount files were sequenced with Illumina using the L3 linker
# adapter recognized by cuatdap: AGATCGGAAGAGC

def get_barcode(wildcards):
    return samples.barcode[samples.names == wildcards.sample].values.tolist()

rule dumultiplex_barcode:
	input:
		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "demultiplexed", "{sample}_demux.fastq.gz")
	params:
		adapter = str(config["adapter"]),
		barcode = get_barcode,
		outdir = os.path.join(outputdir, "demultiplexed"),
		tmp = lambda wildcards: os.path.join(outputdir, "demultiplexed", "".join([wildcards.sample, "_", get_barcode(wildcards)[0], ".fastq.gz"])),
		metrics_file = os.path.join(outputdir, "demultiplexed", "{sample}_demultiplex_metrics.txt")
	log:
		os.path.join(outputdir, "logs", "iCount_demux_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "iCount_demux_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"iCount demultiplex {input.fastq} {params.adapter} {params.barcode} --out_dir {params.outdir} --prefix {wildcards.sample} --file_logpath {log} --results_file {params.metrics_file}; "
		"mv {params.tmp} {output}"


## ------------------------------------------------------------------------------------ ##
## UMI extraction
## ------------------------------------------------------------------------------------ ##
# We remove the samples barcode and random barcode from the reads (first few nucleotides) and place them in the read header

# def get_barcode(wildcards):
#     return samples.barcode_pattern[samples.names == wildcards.sample].values.tolist()
#     # config["sf"][wildcards.sample] 
# # sample = samples.names[samples.type == 'SE'].values.tolist())


# rule umi_tools:
# 	input:
# 		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
# 	output:
# 		os.path.join(outputdir, "FASTQtrimmed", "{sample}_umi_removed.fq.gz")
# 	params:
# 		barcode = get_barcode
# 	log:
# 		os.path.join(outputdir, "logs", "umi_tools_{sample}.log")
# 	benchmark:
# 		os.path.join(outputdir, "benchmarks", "umi_tools_{sample}.txt")
# 	conda:
# 		"envs/environment.yaml"
# 	shell:
# 		"umi_tools extract --stdin={input.fastq} --bc-pattern={params.barcode} --log={log} --stdout {output}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = os.path.join(outputdir, "demultiplexed", "{sample}_demux.fastq.gz")
	output:
		os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz")
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 15 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"


## second round of adapter trimming to remove duplicate adapter ligation events
## we specify the TruSeq adapters taken from https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
# --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
rule trimgalorePE_2nd_round:
	input:
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1.fq.gz"])),
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2.fq.gz"]))
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2} --illumina"



## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}_umi_removed_trimmed.fq.gz") if config["run_trimming"] else os.path.join(outputdir, "demultiplexed", "{sample}_demux.fastq.gz")
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

rule starPE_2nd_round:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "STAR_round2", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align_2nd_round"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"


## Index bam files
rule bamindex:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai")
	log:
		os.path.join(outputdir, "logs", "samtools_index_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "samtools_index_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "STARbigwig", "{sample}_Aligned.sortedByCoord.out.bw")
	params:
		STARbigwigdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "bigwig_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## Normalize the two replicates to the mean library size 
## the scaling factors were computed with R (clipper_analysis.Rmd)
def get_sf(wildcards):
    return config["sf"][wildcards.sample] 

rule normalize_replicates_bigwig:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.r2.ucsc.bam"),
		chrl = "reference/hg38.chrom.sizes_GRCh38.txt"
	output:
		os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.ucsc.replicatesNorm.bw")
	params:
		bigwigdir = os.path.join(outputdir, "bigwig"),
		bedgraph = "{sample}_deduplicated.r2.ucsc.bedGraph",
		sf = get_sf
	log:
		"logs/bigwig_filtered_{sample}.log"
	benchmark:
		"benchmarks/bigwig_{sample}.txt"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -scale {params.sf} -ibam {input.bam} -bg > {params.bigwigdir}/{params.bedgraph}; "
		"LC_ALL='C' sort -k1,1 -k2,2n {params.bigwigdir}/{params.bedgraph} > {params.bigwigdir}/{params.bedgraph}_sorted; "
		"bedGraphToBigWig {params.bigwigdir}/{params.bedgraph}_sorted {input.chrl} {output}; "
		"rm -f {params.bigwigdir}/{params.bedgraph}"
		"rm -f {params.bigwigdir}/{params.bedgraph}_sorted"

## we computed the library sizes with 
## samtools view -c {input.bam}
## and computed the mean library size for each replicate pair 


## ------------------------------------------------------------------------------------ ##
## Segment genome annotation
## ------------------------------------------------------------------------------------ ##
## generate new annotation file with segmented genome

rule indexfasta:
	input:
		config["genome"]
	output:
		config["genome_index"]	
	conda:
		"envs/environment.yaml"
	shell:
		"samtools faidx {input}"
		

rule segment:
	input:
		genome = config["genome_index"],
		gtf = config["segment_gtf"]
	output:
		gtf = config["segment"]
	params:
		metrics_file = os.path.join(outputdir, "logs", "segment_metrics.txt"),
	log:
		"logs/segment.log"
	benchmark:
		"benchmarks/segment.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"iCount segment {input.gtf} {output} {input.genome} --results_file {params.metrics_file} --file_logpath {log}  --report_progress"



## ------------------------------------------------------------------------------------ ##
## Quantify cross-link sites
## ------------------------------------------------------------------------------------ ##
## generate BED file with XL position and number of cDNA molecules per XL

rule xlsites:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
		# segment = config["segment"]
	output:
		os.path.join(outputdir, "xlsites", "{sample}_cDNA_unique.bed")
	params:
		metrics_file = os.path.join(outputdir, "xlsites", "{sample}_xlsites_metrics.txt"),
		multi = os.path.join(outputdir, "xlsites", "{sample}_cDNA_multi.bed"),
		skipped = os.path.join(outputdir, "xlsites", "{sample}_cDNA_skipped.bam")
	log:
		os.path.join(outputdir, "logs", "xlsites_{sample}.log")
	benchmark:
		"benchmarks/xlsites_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"iCount xlsites {input.bam} {output} {params.multi} {params.skipped} --quant cDNA --group_by start  --results_file {params.metrics_file} --file_logpath {log}"
		# "iCount xlsites {input.bam} {output}  NNNGGCGNN_cDNA_multiple.bed NNNGGCGNN_cDNA_skipped.bam --quant cDNA --segment {input.segment} --group_by start  --results_file {params.metrics_file} --file_logpath {log}"




## ------------------------------------------------------------------------------------ ##
## Annotate cross-link sites
## ------------------------------------------------------------------------------------ ##
## annotate each XL site according to biotype

rule annotate:
	input:
		xl = os.path.join(outputdir, "xlsites", "{sample}_cDNA_unique.bed"),
		segment = config["segment"]
	output:
		os.path.join(outputdir, "annotate", "{sample}_biotype.tab")
	params:
		metrics_file = os.path.join(outputdir, "annotate", "{sample}_annotate_biotype_metrics.txt")
	log:
		os.path.join(outputdir, "logs", "annotate_biotype_{sample}.log")
	benchmark:
		"benchmarks/annotate_biotype_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"iCount annotate {input.segment} {input.xl} {output} --results_file {params.metrics_file} --file_logpath {log}"










## ------------------------------------------------------------------------------------ ##
## PCR duplicate removal
## ------------------------------------------------------------------------------------ ##

rule remove_PCR_duplicates:
	input:
		os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.bam")
	params:
		metrics_file = os.path.join(outputdir, "BAM_deduplicated/{sample}_marked_dup_metrics.txt")
	shell:
		"picard MarkDuplicates "
		"I={input} "
		"O={output} "
		"M={params.metrics_file} "
		"REMOVE_DUPLICATES=true "


## Index bam files
rule dedupidx:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated/{sample}.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}.bam.bai")
	shell:
		"samtools index {input.bam}"



## ------------------------------------------------------------------------------------ ##
## Filter second read in pair
## ------------------------------------------------------------------------------------ ##
## the read-pair orientation is F2R1, the forward read is the second in the pair
rule filter_second_read:
	input:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.bam")
	output:
		os.path.join(outputdir, "BAM_deduplicated/{sample}_deduplicated.r2.bam")
	shell:
		"samtools view -f 128 -b -o {output} {input}"

## Convert BAM files to bigWig
rule bigwig_second_read:
	input:
		bam = os.path.join(outputdir, "BAM_deduplicated", "{sample}_deduplicated.r2.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "bigwig", "{sample}_deduplicated.r2.bw")
	params:
		bigwigdir = os.path.join(outputdir, "bigwig"),
		bedgraph = "{sample}_deduplicated.r2.bedGraph"
	log:
		os.path.join(outputdir, "logs", "bigwig_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg > {params.bigwigdir}/{params.bedgraph}; "
		"LC_ALL='C' sort -k1,1 -k2,2n {params.bigwigdir}/{params.bedgraph} > {params.bigwigdir}/{params.bedgraph}_sorted; "
		"bedGraphToBigWig {params.bigwigdir}/{params.bedgraph}_sorted {input.chrl} {output}; "
		"rm -f {params.bigwigdir}/{params.bedgraph};6M1_200219_L1-2_deduplicated.r2.bam "
		"rm -f {params.bigwigdir}/{params.bedgraph}_sorted"


## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
