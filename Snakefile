# make sure all files are named as 
# Tumor = COGN415-NMYC-20171205
# control = COGN415-NMYC-20171205-Input

workdir: "/home/patelk26/KP/test-chipseq-snakemake/"

SAMPLES = ["COGN415-NMYC-20171205",
"COGN415-NMYC-20171205-Input",
"Kelly-MYCN-20150914",
"Kelly-MYCN-20150914-Input",
"LAN5-NMYC-20171205",
"LAN5-NMYC-20171205-Input",
"NB1643-NMYC-20171205",
"NB1643-NMYC-20171205-Input"]

MACS2 = ["COGN415-NMYC-20171205","Kelly-MYCN-20150914","LAN5-NMYC-20171205","NB1643-NMYC-20171205"]

#*****----------------------------------------- TOOLS & SCRIPTS ---------------------------------------------*****#

PHRED_DETECTOR='/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/apps/bioinfo/scripts/phredDetector.pl'
TRIMGALORE='/home/patelk26/miniconda3/bin/trim_galore'
CUTADAPT='/home/patelk26/miniconda3/bin/cutadapt'
BWA='/mnt/isilon/maris_lab/target_nbl_ngs/JoLynne/install/bwa-0.7.12/bwa'
HG19='/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/bwa_index_files/bwa_index_0.7.12/hg19.fa'
BWA_SAMSE='/home/patelk26/miniconda3/bin/bwa'
SAMTOOLS='/home/patelk26/miniconda3/bin/samtools'
PICARD='/cm/shared/apps_chop/picard/1.140/picard.jar'
BEDTOOLS='/home/patelk26/miniconda3/bin/bedtools'
MASC='/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/apps/masc/masc.pl'
CHROM_SIZES='/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/ucsc-download/hg19.chrom.sizes'
BLACKLISTED_REGIONS='/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/ENCODE_download/blacklisted_regions/hg19-blacklistMergedList_noheader.bed'
BEDGRAPH_BIGWIG='/home/patelk26/miniconda3/bin/bedGraphToBigWig'

#*****-------------------------------------------------------------------------------------------------------*****#

rule all:
	input: 
		expand('peakfiles/bigwig/{macs}.IP.sorted.bw', macs=MACS2),
		expand('peakfiles/bigwig/{macs}.Input.sorted.bw', macs=MACS2)


rule idx_phred:
    input:
        expand('catfiles/{sample}.fastq.gz', sample=SAMPLES)
    output:
        "qcfiles/{sample}.phred.txt"
    shell:
        "{PHRED_DETECTOR} {input} > {output}"

#####---------------------------------- QC fastq ----------------------------------------------####

rule trim_qc:
    input:
        fastq = "catfiles/{sample}.fastq.gz",
        phred = "qcfiles/{sample}.phred.txt"

    output:
        output1 ='trim_fq/{sample}_trimmed.fq.gz',
        output2 ="trim_fq/{sample}.fastq.gz_trimming_report.txt"

    shell:
        """
        {TRIMGALORE} --phred$(cat {input.phred}) --fastqc --fastqc_args '--outdir fastqc/' -path_to_cutadapt {CUTADAPT} {input.fastq}  -o trim_fq/
        """

rule unzip:
	input:
		fq_gz='trim_fq/{sample}_trimmed.fq.gz'
	output:
		'trim_fq/{sample}_trimmed.fq'
	shell:
		"gunzip {input.fq_gz}"
#####---------------------------------- BWA mapping, sorting and indexing ----------------------------------------------####

rule bwa_aln:
	input:
		fq_gunzipped = 'trim_fq/{sample}_trimmed.fq',
		phred = 'qcfiles/{sample}.phred.txt'
	output:
		'aligned/{sample}.trim.sai'
	shell:
		"""
		pscale=$(cat {input.phred})

		if [ pscale == 64 ]
		then 
			{BWA} aln {HG19} -I -q 20 -t 16 {input.fq_gunzipped} -f {output}
		else
			{BWA} aln {HG19} -q 20 -t 16 {input.fq_gunzipped} -f {output}

		fi
		"""


rule bwa_samse:
	input:
		sai = 'aligned/{sample}.trim.sai',
		trim_fq = 'trim_fq/{sample}_trimmed.fq'
	output:
		'aligned/{sample}.sam'
	shell:
		"""
		{BWA_SAMSE} samse -f {output} {HG19} {input.sai} {input.trim_fq} 
		"""


rule sam_to_bam:
	input:
		sam = 'aligned/{sample}.sam'
	output:
		'aligned/{sample}.trim.bam'
	shell:
		"{SAMTOOLS} view -Sb {input.sam} > {output}"



rule bam_sort_index:
	# sorting BAM files and indexing
	# #CleanSam will mark all unmapped reads as MAPQ = 0 (this error pops up with the remove duplicates tool unless you select LENIENT stringency, 
	# however, if not marked as 0, may cause downstream errors, so best to mark now)
	input:
		bam = 'aligned/{sample}.trim.bam'

	output:
		bam_sorted = 'aligned/{sample}.trim.sorted.bam',
		bam_clean_sorted = 'aligned/{sample}.trim.clean.sorted.bam'

	shell:
		"""
		{SAMTOOLS} sort {input.bam} -o {output.bam_sorted}
		{SAMTOOLS} index {output.bam_sorted}
		java -jar {PICARD} CleanSam INPUT={output.bam_sorted} OUTPUT={output.bam_clean_sorted}
		rm {input.bam}
		"""


rule bam_process:
	# remove duplicates -- memory requirements for removing duplicates -- set at 16, 20, 4
	# generate fragment lengths -- direct output for this to file and then get the fragment size (NOTE: Load module perl/5.26.1 for this step to work)
	input:
		in_bam = 'aligned/{sample}.trim.clean.sorted.bam'

	output:
		bam_clean_dedup = 'aligned/{sample}.trim.clean.dedup.sorted.bam',
		dupmetrics = 'qcfiles/{sample}.dupmetrics.txt', 
		map75log = 'qcfiles/{sample}.map75log.txt', 
		fraglen = 'qcfiles/{sample}.fragmentlength.txt',
		bed = 'aligned/{sample}.bed'

	shell:
		"""
		java -jar {PICARD} MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true INPUT={input.in_bam} OUTPUT={output.bam_clean_dedup} METRICS_FILE={output.dupmetrics}
		# convert bam to bed and use those reads to estimate fragment length
		{BEDTOOLS} bamtobed -i {output.bam_clean_dedup} > {output.bed}
		# generate fragment lengths -- direct output for this to file and then get the fragment size (NOTE: Load module perl/5.26.1 for this step to work)
		module load perl/5.26.1
		perl {MASC} --verbose --mappability_path=/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/ucsc-download/mappability/hg19-align75-mapfiles/binary/ --chrom_length_file={CHROM_SIZES} --input_bed={output.bed} --prefix={output.map75log}.map75 --smooth_win_size=15 --min_shift=0 --max_shift=400 > {output.map75log}
		# make new file for fragment length
		grep 'MaSC' {output.map75log} | awk '{{print $3}}' | cut -d':' -f2 > {output.fraglen}
		rm {input.in_bam}
		"""

#####---------------------------------- Creating tag directories ----------------------------------------------####

rule chip_density:
	input:
		fraglen = 'qcfiles/{sample}.fragmentlength.txt',
		bed = 'aligned/{sample}.bed'

	output:
		tag_dir = 'tag_dir/{sample}/'
	shell:
		"""
		makeTagDirectory {output.tag_dir} -fragLength $(cat {input.fraglen}) -format bed {input.bed}"
		"""

#####---------------------------------- MACS2 peak calling and post-processing ----------------------------------------------####

rule macs2:
	input:
		tumor = 'aligned/{macs}.trim.clean.dedup.sorted.bam',
		control = 'aligned/{macs}-Input.trim.clean.dedup.sorted.bam',
		fraglen = 'qcfiles/{macs}.fragmentlength.txt',
		bed = 'aligned/{macs}.bed'

	output:
		peaks = 'peakfiles/{macs}_peaks.xls',
		narrowPeaks = 'peakfiles/{macs}_peaks.narrowPeak',
		summits = 'peakfiles/{macs}_summits.bed',
		bdg_control = 'peakfiles/{macs}_control_lambda.bdg',
		bdg_tumor = 'peakfiles/{macs}_treat_pileup.bdg'

	shell:
		"""
		# running MACS2
		macs2 callpeak -t {input.tumor} -c {input.control} -f BAM -g hs -n {wildcards.macs} --outdir peakfiles/ --nomodel -B --verbose 3 --call-summits --tempdir tmp/ --extsize $(cat {input.fraglen}) --SPMR
		"""

rule macs2_postprocess:
	input:
		peakfile = 'peakfiles/{macs}_peaks.narrowPeak',
		bdg_control = 'peakfiles/{macs}_control_lambda.bdg',
		bdg_tumor = 'peakfiles/{macs}_treat_pileup.bdg'

	output:
		filtered_peak = 'peakfiles/{macs}.filtered.narrowPeak',
		bdg_sorted_tumor = 'peakfiles/{macs}.IP.sorted.bdg',
		bdg_sorted_control = 'peakfiles/{macs}.Input.sorted.bdg',
		IP_bigwig = 'peakfiles/bigwig/{macs}.IP.sorted.bw',
		Input_bigwig = 'peakfiles/bigwig/{macs}.Input.sorted.bw'


	shell:
		"""
		# remove duplicates
		bedtools intersect -v -a {input.peakfile} -b {BLACKLISTED_REGIONS} > {output.filtered_peak}
		# sort bedgraph files
		sort -k1,1 -k2,2n {input.bdg_tumor} > {output.bdg_sorted_tumor}
		sort -k1,1 -k2,2n {input.bdg_control} > {output.bdg_sorted_control}

		#convert bedgraph to bigwig
		{BEDGRAPH_BIGWIG} {output.bdg_sorted_tumor} /mnt/isilon/maris_lab/target_nbl_ngs/JoLynne/install/IGVTools/genomes/hg19.chrom.sizes {output.IP_bigwig}
		{BEDGRAPH_BIGWIG} {output.bdg_sorted_control} /mnt/isilon/maris_lab/target_nbl_ngs/JoLynne/install/IGVTools/genomes/hg19.chrom.sizes {output.Input_bigwig}

		rm {input.bdg_control}
		rm {input.bdg_tumor}

		"""

