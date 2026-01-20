#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk / g.selten@uu.nl
# Created Date  : 07/11/2022
# version       : '1.2'
# ---------------------------------------------------------------------------
# Pipeline (bash) to perform fingerprint identification from microbiome data.
# ---------------------------------------------------------------------------

function logo() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Arrays
	array_keep=("None" "BAMs" "All")
	array_verbose=("Quiet" "Sample" "All")
	array_force=("None" "All" "Skipped" "Failed")

	# Variables
	if [[ $1 != "help" && $8 != 0 ]]; then
		folder=$1
		target=$2
		lendev=$3
		cutoff=$4
		threads=$5
		mem=$6
		keepfiles=${array_keep[${7}]}
		verbose=${array_keep[${8}]}
		force=${array_force[${9}]}

		# Logo
		echo ""
		echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
		echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
		echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
		echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
		echo ""
		echo "${bold}----------------------------------------------${normal}"
		echo "${bold}Summary:${normal}"
		echo "    Folder: ${folder}"
		echo "    Target: ${target}"
		echo "    Cutoff: ${cutoff}"
		echo "    Deviation length: ${lendev}"
		echo "    Threads: ${threads}"
		echo "    Memory: ${mem} GB"
		echo "    Keep files: ${keepfiles}"
		echo "    Verbose: ${verbose}"
		echo "    Force: ${force}"
		echo ""
		printf "${bold}Start:${normal}\n"
		date "+    date: %d/%m/%Y"
		date "+    time: %H:%M:%S"
		echo ""
		echo "${bold}Progress:${normal}"
		checkProgress ${folder}
		echo "${bold}----------------------------------------------${normal}"
		echo ""
	elif [[ $8 != 0 ]]; then
		# Logo
		echo ""
		echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
		echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
		echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
		echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
		echo ""
	fi
}

function usage()
{
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Logo
	logo help

	# Usage
	echo "Usage: SyFi.sh main -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>"
	printf "\n"
	echo "${bold}REQUIRED:${normal}"
	echo "# Input"
	echo "  -i | --input-folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder-structure>."
	echo "  -s | --search-target    Genomic region of interest in fasta format, e.g., 16S."
	printf "\n"
	echo "${bold}OPTIONAL:${normal}"
	echo "# Haplotype deviation:"
	echo "  -l | --len-deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (defaut: 100 bp)."
	echo "  -c | --cutoff            Maximum ratio deviation between haplotypes per sample. This parameter defined how much can an haplotype deviate from the minimum haplotype ratio (default: 25)."
	printf "\n"
	echo "# Input extension:"
	echo "  --fasta-extension        Reference file extension (default: fasta)."
	echo "  --fastq-extension        Illumina reads file extension (default: fastq.gz)."
	printf "\n"
	echo "# Computation:"
	echo "  -t | --threads          Number of threads (default: 1)."
	echo "  -m | --memory           Memory in GBs (default: 8GB)."
	printf "\n"
	echo "# Output options:"
	echo "  -o  | --OUTPUT_FOLDER 						 Optional output folder, otherwise it will be written to current working directory"
	echo "  -k | --keep-files       Keep temporary files [0: Minimum, 1: BAM's, or 2: All] (default: 0)."
	echo "  -v | --verbose          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2)."
	echo "  -f | --force            Force re-computation of computed samples [0: None, 1: All, 2: Skipped, or 3: Failed] (default: 0)."
	printf "\n"
}


### Variables

# Colors
red=$(tput setaf 1)
yellow=$(tput setaf 220)
blue=$(tput setaf 27)
green=$(tput setaf 10)
normal=$(tput sgr0)

# Default variables
THREADS=1
FAEXT='fasta'
FQEXT='fastq.gz'
MEM=8
BPDEV=300
CUTOFF=25
FORCE=0
KEEPF=0
VERBOSE=2

# Define software path
SYFI_BASE=$(dirname $(whereis ${0} | cut -d " " -f 2))
GENOME_PHASE=$(echo ${SYFI_BASE}/src/genome_phase.sh)
FILTER_HAPLOTYPES=$(echo ${SYFI_BASE}/src/filterHaplotypes.R)
INTEGRATION=$(echo ${SYFI_BASE}/src/Integration.R)

### Parameters

# display usage if 
if [[ $# -lt 4  && $1 != "-h" && $1 != "--help" && $1 != "-c" && $1 != "--citation" && $1 != "--folder-structure" ]]; then
	printf "\n${red}ERROR:${normal} You must provided at least the assembly and alignment files.\n"
	usage
	exit 2
fi

# Default variables
OUTPUT_FOLDER="." # Default to current directory

# Get parameters
while [[ $# -gt 0 ]]; do
	case $1 in
		-o|--OUTPUT_FOLDER)
			OUTPUT_FOLDER="$2"
			shift 2
			;;
		-i|--input-folder)
			INPUT_FOLDER="$2"
			shift 2
			;;
		-s|--search-target)
			SEARCH_TARGET="$2"
			shift 2
			;;
		-l|--len-deviation)
			BPDEV="$2"
			shift 2
			;;
		-c|--cutoff)
			CUTOFF="$2"
			shift 2
			;;
		--fasta-extension)
			FAEXT="$2"
			shift 2
			;;
		--fastq-extension)
			FQEXT="$2"
			shift 2
			;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-m|--mem)
			MEM="$2"
			shift 2
			;;
		-k|--keep-files)
			KEEPF="$2"
			shift 2
			;;
		-v|--verbose)
			VERBOSE="$2"
			shift 2
			;;
		-f|--force)
			FORCE="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit
			;;
		--citation)
			citation
			exit
			;;
		--folder-structure)
			folder_structure
			exit
			;;
		*)
			echo "Unknown parameter: $1"
			usage
			exit 1
			;;
	esac
done


# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
	printf "\n${red}Execution halted by user.${normal}\n"
	if [ -f "${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt" ]; then
		printf "\n${red}Execution halted by user.${normal}\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	fi
	exit
}

### Checks

# Check: Software dependencies
software=(blastn bwa-mem2 samtools gzip spades.py seqtk seqkit kallisto Rscript gatk3 picard plot-bamstats bcftools)
for pckg in ${software[@]}
do 
	type ${pckg} 2> /dev/null 1>&2 
	if [ $? != 0 ]
	then
		printf "\n${red}ERROR:${normal} ${pckg} missing. Install or activate SyFi conda environment.\n\n"
		exit
fi
done

# Check: Input folder
if [[ ! -d ${INPUT_FOLDER} ]]; then
	printf "\n${red}ERROR:${normal} Folder ${INPUT_FOLDER} missing.\n\n"
	exit
fi

# Check: Keep file
if [[ ! "$KEEPF" =~ ^[0-9]+$ || ${KEEPF} -gt 2 || ${KEEPF} -lt 0 ]]; then
	printf "\n${red}ERROR:${normal} Keep temporary files should be a number between 0 and 2 [0: None, 1: BAM's, or 2: All].\n\n"
	exit
fi

# Check: Force value
if [[ ! "$FORCE" =~ ^[0-9]+$ || ${FORCE} -gt 3 || ${FORCE} -lt 0 ]]; then
	printf "\n${red}ERROR:${normal} Force value should be a number between 0 and 3 [0: None, 1: All, 2: Skipped or 3: Failed].\n\n"
	exit
fi

# Check: Length deviation above zero
	tl=$(grep -v "^>" ${SEARCH_TARGET} | wc -c)
	min_tl=$((tl-${BPDEV}))
	if [ ${min_tl} -le 0 ]; then
		printf "\n${red}ERROR:${normal} length deviation is equal or larger than the target length (target length: ${tl}).\n\n"
		exit
	fi

# Check: progress table
function checkProgress() {
	# Colors
	r=$(tput setaf 1)
	y=$(tput setaf 220)
	g=$(tput setaf 10)
	n=$(tput sgr0)

	# Variables
	INPUT_FOLDER=$1

	# Check file
	if [ -e progress.txt ]; then
		total=$(ls $INPUT_FOLDER | wc -l)
		success=$(grep "Success" progress.txt | wc -l)
		skipped=$(grep "Skipped" progress.txt | wc -l)
		failed=$(grep "Failed" progress.txt | wc -l)
		printf "\tSuccess: ${g} ${success} ${n}\n\tSkipped: ${y} ${skipped} ${n}\n\tFailed: ${r}  ${failed} ${n}\n\tTotal:    ${total}\n"
	else
		total=$(ls $INPUT_FOLDER | wc -l)
		printf "\tSuccess: ${g} 0 ${n}\n\tSkipped: ${y} 0 ${n}\n\tFailed: ${r}  0 ${n}\n\tTotal:    ${total}\n"
	fi

}

### Functions

## FUNCTION: Variant Calling
function variantCalling() {
	# Variables
	bam_file=$1
	reference_fasta=$2
	reference_dict=$3
	input_bam=$4
	OUTPUT_FOLDER=$5
	threads=$6

	# Intro messages
	printf "\n"
	echo "Performing Variant Calling"

	# Mark duplicates (Mapping is already done)
	printf "\n"
	echo "Mark BAM duplicates"
	markDuplicates ${bam_file} ${OUTPUT_FOLDER} ${threads}

	# Haplotype caller
	printf "\n"
	echo "Haplotype calling"
	haplotypeCaller ${reference_fasta} ${reference_dict} ${OUTPUT_FOLDER} ${input_bam}

	# Join genotyping
	printf "\n"
	echo "Join genotyping"
	jointGenotype ${reference_fasta} ${OUTPUT_FOLDER} ${threads}

}

# FUNCTION: Mark Duplicates
function markDuplicates() {
	# Arguments
	bam_file=$1
	OUTPUT_FOLDER=$2
	threads=$3

	# Variables
	sample=$(basename ${bam_file} | sed 's/.rebuild.sort.bam//g')
	output_fn="${OUTPUT_FOLDER}/mapped_filtered/${sample}_stats.txt"

	# Execution
	picard MarkDuplicates -I ${bam_file} -O ${OUTPUT_FOLDER}/mapped_filtered/${sample}.filtered.bam -M ${OUTPUT_FOLDER}/mapped_filtered/${sample}.filtered.bam-metrics.txt -AS --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR ${OUTPUT_FOLDER}/mapped_filtered

	# Plot BAM
	samtools stats -d ${bam_file} > ${output_fn}
	plot-bamstats -p ${OUTPUT_FOLDER}//mapped_filtered/ ${output_fn}
}

# FUNCTION: Haplotype Caller
function haplotypeCaller() {
	# Arguments
	reference_fasta=$1
	reference_dict=$2
	OUTPUT_FOLDER=$3
	input_bam=$4

	# Variables (-t genomic -m joint -a False -i False)
	sample=$(basename ${input_bam} | cut -d "." -f 1)
	erc_mode="GVCF"
	min_thr=30
	output_mode="EMIT_VARIANTS_ONLY"
	ploidy=2
	
	# Create dictionary
	picard CreateSequenceDictionary -R ${reference_fasta} -O ${reference_dict}

	# Add read group to BAM
	picard AddOrReplaceReadGroups -I ${input_bam} -O ${OUTPUT_FOLDER}/genotyped/${sample}.filtered.readgroup.bam -LB lib1 -PL ILLUMINA -PU unit1 -SM ${sample}
	samtools index -b ${OUTPUT_FOLDER}/genotyped/${sample}.filtered.readgroup.bam ${OUTPUT_FOLDER}/genotyped/${sample}.filtered.readgroup.bai -@ ${threads}

	# Index fasta
	samtools faidx ${reference_fasta}

	# Haplotypes
	gatk3 -T HaplotypeCaller -ERC ${erc_mode} -ploidy ${ploidy} -stand_call_conf ${min_thr} -I ${OUTPUT_FOLDER}/genotyped/${sample}.filtered.readgroup.bam -o ${OUTPUT_FOLDER}/genotyped/${sample}.g.vcf.gz -R ${reference_fasta} --output_mode ${output_mode}
}

# FUNCTION: Joint Genotyping
function jointGenotype() {
	# Arguments
	reference_fasta=$1
	OUTPUT_FOLDER=$2
	threads=$3

	# Variables
	sample=$(basename ${reference_fasta} | sed 's/.fasta//g')
	input_gvcf="${OUTPUT_FOLDER}/genotyped/${sample}.g.vcf.gz"
	intervals_fn="${OUTPUT_FOLDER}/reference/intervals.list"
	workspace_dir="${OUTPUT_FOLDER}/variants/db_workspace"
	merge_intervals="--merge-input-intervals"
	database=${workspace_dir}
	output_vcf="${OUTPUT_FOLDER}/variants/${sample}.vcf.gz"
	comp_fn="${OUTPUT_FOLDER}/variants/${sample}.vchk"
	plot_dir="${OUTPUT_FOLDER}/variants/${sample}/"

	# Create intervals
	cat ${OUTPUT_FOLDER}/20-Alignment/${sample}/${sample}.fasta.fai | awk '{print $1":1-"$2}' > ${intervals_fn}

	# Run genotype
	gatk3 -T GenotypeGVCFs -V ${input_gvcf} -R ${reference_fasta} -o ${output_vcf} -L ${intervals_fn} -G StandardAnnotation 

	# VCF Stats
	bcftools stats -F ${reference_fasta} -s- ${output_vcf} > ${comp_fn}
}

# FUNCTION: Copy Number
function copyNumber() {
	# Arguments
	INPUT_FOLDER=$1
	subf=$2

	if [ ${VERBOSE} -eq 2 ]; then printf "Copy number; "; fi

	# Log
	printf "\n\n### Target Copy Number ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

	# Header
	printf "Strain\tGenome_length\tGenome_nbases\tTarget_length\tTarget_nbases\tCopy_number\n" > ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv

	# Variables
	# Get length of assembly
	lgen=$(cat ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} | grep -v "^>" | tr -d "\n" | wc -c)
	# Get number of bases in assembly reads (Speed up?)
	if [ ${FQEXT} == "fastq.gz" ] || [ ${FQEXT} == "fq.gz" ]; then
		braw=$(zcat ${INPUT_FOLDER}/${subf}/${subf}_R[12].${FQEXT} | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)
	else
		braw=$(cat ${INPUT_FOLDER}/${subf}/${subf}_R[12].${FQEXT} | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)
	fi

	# Get length of longest recovered target without flanking regions
	l16S=$(cut -f 15 ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | awk 'BEGIN{a=0} {if ($1>0+a) a=$1} END{print a}')

	# Define start and end of target for base count (Choose first one if l16S returns multiple (equal size))
	start=$(cut -f 7,15 ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | head -n 1 | awk -v l16S=${l16S} '{if ($2 == l16S) {print $1}}')
	end=$(cut -f 8,15 ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | head -n 1 | awk -v l16S=${l16S} '{if ($2 == l16S) {print $1}}')
	
	# Number of bases in selected reads
	b16S=$(bedtools genomecov -d -ibam ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sort.bam | awk -v start=${start} '($2 > start)' | awk -v end=${end} '($2 < end)' | awk '{sum+=$3;} END{print sum;}')

	# Compute ratio (round <1 to 1)
	cnum=$(echo "(${b16S}/${l16S}) / (${braw}/${lgen})" | bc -l)
	if (( $(echo "${cnum} < 0.5" | bc -l) )); then
		cnum="1"
	elif (( $(echo "${cnum} > 25.5" | bc -l) )); then
		cnum="25"
	fi

	# File
	printf "${subf}\t${lgen}\t${braw}\t${l16S}\t${b16S}\t${cnum}\n" >> ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv
}

# FUNCTION: Fingerprint
function fingerPrint() {
	# Variables
	mode=$1 # unique/multiple
	subf=$2
	variants=$3 # Yes/No

	# Create folder
	mkdir -p ${OUTPUT_FOLDER}/70-Fingerprints/${subf}

	if [ $mode == "unique" ];then
		# Log
		if [ ${VERBOSE} -eq 2 ]; then printf "Fingerprint [unique]\n"; fi
		printf "\n\n### Fingerprint ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# Create haplotype sequence
		printf ">seq_h1\n$(grep -v "^>" ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta)\n" > ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/seq_h1.fasta

		# Copy recovered target
		cp ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta

		# Rename fasta header
		cat ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta | sed "s/^>NODE.*/>${subf}_all_haplotypes/g" >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta

	elif [ $mode == "multiple" ];then
		# Log
		if [ ${VERBOSE} -eq 2 ]; then printf "Fingerprint\n"; fi
		printf "\n\n### Fingerprint ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# Check: Variants (Yes/No)
		if [ $variants == "Yes" ]; then

			# Print statistics of variant calling
			printf "#Position - Position of the called variant in the reference sequence\n#Reference - Nucleotide in the reference sequence in that position\n#Alternative - Alternative nucleotide in that same position (variant)\n#FS - Fisher Strand Bias - Statistic test to evaluate if one DNA strand is favored over the other - high values indicate strand-specific errors (> 60 for SNPs and > 200 for InDels)\n#SOR - Strand Odds Ratio - Statistic test to evaluate if one DNA strand is favored over the other - Values > 3.0 may indicate errors\n#MQ - Mapping Quality - Average quality of mapped reads - 60 is perfect while lower values reduce confidence\n#MQRankSum - Mapping quality comparison between reference and alternative variants - Near 0 is good, otherwise it may indicate a potential bias (< -12.5 often considered bad)\n#BaseQRankSum - Base quality comparison between reference and alternative reads. Large negative or positive values suggest sequencing bias towards reference (negative) or alternative variant (positive). Values > 2 or < -2 indicate sequencing bias to one of the sequences.\n#ReadposRankSum - Comparison of location of variants (start, middle, end) of sequence between reference and alternative sequences. Values > 8 or < -8 are considered untrustworthy variants.\n" >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
			printf "\nPosition\tReference\tAlternative\tFS\tSOR\tMQ\tMQRankSum\tBaseQRankSum\tReadposRankSum\n" >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv

			number_of_variants=$(zgrep "BaseQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | grep -v "#" | cut -f2 | wc -l)

			for i in $(seq 1 $number_of_variants) ; do
    			    zgrep "BaseQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | grep -v "#" | cut -f2 | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
    			    zgrep "BaseQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | grep -v "#" | cut -f4 | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
    		    	    zgrep "BaseQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | grep -v "#" | cut -f5 | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
    		    	    zgrep "FS" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | grep "FS" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
	    		    zgrep "SOR" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | sed 's/\t/\n/g' | grep "SOR" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
  	  		    zgrep "MQ" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | grep "MQ" | grep -v "MQRankSum" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
    			    zgrep "MQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | grep "MQRankSum" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >>  ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
	    		    zgrep "BaseQRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | grep "BaseQRankSum" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\t/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
    			    zgrep "ReadPosRankSum" ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_phased.vcf.gz | sed 's/;/\n/g' | grep -v "#" | grep "ReadPosRankSum" | cut -f2 -d "=" | head -n $i | tail -n1 | tr -d '\n' | sed 's/$/\n/g' >> ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_variants_stats.tsv
			done

			# Separate fasta into multiples fasta
			seqkit split -i ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta -O ${OUTPUT_FOLDER}/70-Fingerprints/${subf} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# Rename haplotypes files
			for i in $(ls -d ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/*); do mv ${i} $(echo ${i} | sed "s/clean_${subf}_haplotypes.part_//g"); done &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		elif [ $variants == "No" ]; then
			# Separate fasta into multiples fasta
			seqkit split -i ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta -O ${OUTPUT_FOLDER}/70-Fingerprints/${subf} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# Rename haplotypes files
			for i in $(ls -d ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/*); do mv ${i} $(echo ${i} | sed "s/${subf}.part_//g"); done &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		else
			printf "Variants variable is not yes or no but ${variants}" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		fi

		# Keep haplotypes filtered in integration
		keep=($(tail -n +2 ${OUTPUT_FOLDER}/60-Integration/${subf}/integration.tsv | cut -f 1))
		for i in $(ls ${OUTPUT_FOLDER}/70-Fingerprints/${subf} | sed 's/.fasta//g'); do
			if [[ ! "${keep[*]}" =~ "${i}" ]]
			then
				rm -f ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${i}.fasta
			fi
		done

		# Concatenate haplotypes (NF-1 to obtain "per_haplotype")
		cat ${OUTPUT_FOLDER}/60-Integration/${subf}/integration.tsv | awk '{print $1 "/" $(NF-1)}' | tail -n +2 > ${OUTPUT_FOLDER}/60-Integration/${subf}/tmp.tsv
		for h in $(cat ${OUTPUT_FOLDER}/60-Integration/${subf}/tmp.tsv); do
			seq=$(echo ${h} | cut -d "/" -f 1)
			num=$(echo ${h} | cut -d "/" -f 2)
			for n in $(seq ${num}); do
				# End when haplotype ratio > 1
				if [[ ${n} == ${num} && ${n} != 1 ]]; then
					grep -v "^>" ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${seq}.fasta >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
				# End when haplotype ratio == 1 
				elif [[ ${n} == ${num} && ${n} == 1 ]]; then
					grep -v "^>" ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${seq}.fasta >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
				# Addition of haplotypes before ending
				else
					grep -v "^>" ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${seq}.fasta >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
					echo "NNNNNNNNNN" >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
				fi
			done 
		done
		printf ">${subf}_all_haplotypes\n" > ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta
		cat ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta | tr -d "\n" >> ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta
	fi
}

function CleanFiles() {
	# Variables
	SAMPLE=$1
	FOLDER=$2
	KEEPF=$3

	# Array
	keep=("./${OUTPUT_FOLDER}/01-Logsmain/log_${SAMPLE}.txt" "./${OUTPUT_FOLDER}/10-Blast/${SAMPLE}.tsv" "./${OUTPUT_FOLDER}/11-Sequences/${SAMPLE}/${SAMPLE}.fasta" "./${OUTPUT_FOLDER}/20-Alignment/${SAMPLE}/${SAMPLE}.fasta" "./${OUTPUT_FOLDER}/20-Alignment/${SAMPLE}/${SAMPLE}_R1.fastq.gz" "./${OUTPUT_FOLDER}/20-Alignment/${SAMPLE}/${SAMPLE}_R2.fastq.gz" "./${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/variants/${SAMPLE}.vcf.gz" "./${OUTPUT_FOLDER}/40-Phasing/${SAMPLE}/${SAMPLE}_phased.vcf.gz" "./${OUTPUT_FOLDER}/40-Phasing/${SAMPLE}/${SAMPLE}_variants_stats.tsv" "./${OUTPUT_FOLDER}/50-haplotypes/${SAMPLE}/clean_${SAMPLE}_haplotypes.fasta" "./${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/abundance.tsv" "./${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/copy_number.tsv" "./${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/integration.tsv" "./${OUTPUT_FOLDER}/70-Fingerprints/${SAMPLE}/${SAMPLE}_all_haplotypes.fasta")

	# Remove files
	if [ ${KEEPF} -eq 0 ]; then
		# Find and remove
		for file in $(find . -type f -name "*${SAMPLE}*"); do
			if [[ ! "${keep[*]}" =~ "${file}" ]]; then
				if [[ $(echo ${file} | grep ${FOLDER} | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep -e "seq_h" -e "assembly_h" -e ".fasta" | grep -e "40-Phasing" -e "70-Fingerprints" | wc -l) != 0 ]]; then
					continue
				else
					rm -rf ${file}
				fi
			fi
		done

		# Remove non-captured folders
		rm -rf ${OUTPUT_FOLDER}/20-Alignment/${SAMPLE}/spades ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/genotyped ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/variants/db_workspace ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/reference ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/mapped_filtered
		# Remove empty directories
		find . -type d -empty -delete
		# Remove non-captured files
		rm -rf ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/variants/tmp_*.config ${OUTPUT_FOLDER}/40-Phasing/${SAMPLE}/${SAMPLE}_phased.vcf ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/abundance.h5 ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/run_info.json ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/tmp.tsv ${OUTPUT_FOLDER}/70-Fingerprints/${SAMPLE}/${SAMPLE}_all_haplotypes.tmp.fasta

	elif [ ${KEEPF} -eq 1 ]; then
		# Find and remove
		for file in $(find . -type f -name "*${SAMPLE}*"); do
			if [[ ! "${keep[*]}" =~ "${file}" ]]; then
				if [[ $(echo ${file} | grep ${FOLDER} | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep -e "seq_h" -e "assembly_h" -e ".fasta" | grep -e "40-Phasing" -e "70-Fingerprints" | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep ".sort.bam" | wc -l) != 0 ]]; then
					continue
				else
					rm -rf ${file}
				fi
			fi
		done    

		# Remove non-captured folders
		rm -rf ${OUTPUT_FOLDER}/20-Alignment/${SAMPLE}/spades ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/genotyped ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/variants/db_workspace ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/reference ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/mapped_filtered
		# Remove empty directories
		find . -type d -empty -delete
		# Remove non-captured files
		rm -rf ${OUTPUT_FOLDER}/30-VariantCalling/${SAMPLE}/variants/tmp_*.config ${OUTPUT_FOLDER}/40-Phasing/${SAMPLE}/${SAMPLE}_phased.vcf ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/abundance.h5 ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/run_info.json ${OUTPUT_FOLDER}/60-Integration/${SAMPLE}/tmp.tsv ${OUTPUT_FOLDER}/70-Fingerprints/${SAMPLE}/${SAMPLE}_all_haplotypes.tmp.fasta
	fi
}

### Execution

#-------------- #
# I. Haplotypes #
#-------------- #

# Call logo
logo ${INPUT_FOLDER} ${SEARCH_TARGET} ${BPDEV} ${CUTOFF} ${THREADS} ${MEM} ${KEEPF} ${VERBOSE} ${FORCE}

for subf in $(ls ${INPUT_FOLDER}); do
 
	## CHECK: Redirect workflow given status and force
	if [ -f progress.txt ]; then
		status=$(grep -w ${subf} progress.txt | cut -f 2)
		if [[ ${status} == "Success" ]] && [[ ${FORCE} != 1 ]]; then
			continue
		elif [[ ${status} == "Skipped" ]] && [[ ${FORCE} != 2 ]]; then
			continue
		elif [[ ${status} == "Failed" ]] && [[ ${FORCE} != 3 ]]; then
			continue
		elif [[ ${status} == "" ]] && [[ ${FORCE} != 0 ]]; then
			continue
		fi

		# CHECK: Avoid re-computation
		if [[ ${FORCE} != 0 ]]; then
			# Remove line in progress.txt if exists already
			if [[ $(grep -w ${subf} progress.txt | wc -l) -eq 1 ]]; then
				grep -wv ${subf} progress.txt > tmp
				mv tmp progress.txt
			fi
			
			# Remove results for sample
			rm ${OUTPUT_FOLDER}/10-Blast/${subf}.tsv
			for fld in 11-Sequences 20-Alignment 30-VariantCalling 40-Phasing 50-Haplotypes 60-Integration 70-Fingerprints; do
				rm -rf ${fld}/${subf}
			done
		fi
	fi

	## CHECK: Log file
	if [[ -f ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt ]]; then
		# Remove old log file
		rm -rf ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		touch ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	else
		# Touch Log File
		mkdir -p ${OUTPUT_FOLDER}/01-Logsmain
		touch ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	fi

	## CHECK: Input in folder
	# Fasta input
	if [[ ! -f ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} ]]; then
		printf "\n${red}ERROR:${normal} File ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} missing, please use the \"-e/--extension\" argument if the extension is not \"fasta\"." | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} computation will be skipped.\n"; fi
		printf "${subf}\tSkipped\tWrong reference extension\n" >> progress.txt
		continue
	fi

	# Fastq input
	if [[ ! -f ${INPUT_FOLDER}/${subf}/${subf}_R1.${FQEXT} || ! -f ${INPUT_FOLDER}/${subf}/${subf}_R2.${FQEXT} ]]; then
		printf "\n${red}ERROR:${normal} One or both illumina reads ${INPUT_FOLDER}/${subf}/${subf}_R[1/2].${FQEXT} are missing, please use the \"-e/--extension\" argument if the extension is not \"fastq.gz\"." | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} computation will be skipped.\n"; fi
		printf "${subf}\tSkipped\tWrong reads extension\n" >> progress.txt
		continue
	fi

	if [ ${VERBOSE} -ge 1 ]; then date "+start time: %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
	if [ ${VERBOSE} -ge 1 ]; then printf "${blue}Sample:${normal} $subf\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi

	## ------------------------------------
	## Target recovery from contigs (BLAST)
	## ------------------------------------

	if [ ${VERBOSE} -eq 2 ]; then printf "Performing: Mapping; "; fi

	# Create folder
	mkdir -p ${OUTPUT_FOLDER}/10-Blast ${OUTPUT_FOLDER}/11-Sequences/${subf}
	# Blastn
	blastn -query ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} -subject ${SEARCH_TARGET} -strand both -outfmt "6 std qseq" > ${OUTPUT_FOLDER}/10-Blast/${subf}.tsv
	printf ">${subf}\n$(cat ${OUTPUT_FOLDER}/10-Blast/${subf}.tsv | sort -n -k4 | tail -n 1 | cut -f 13 | sed 's/-//g')" > ${OUTPUT_FOLDER}/11-Sequences/${subf}/${subf}.fasta

	# CHECK: Absent target match (Perhaps, remove small hits!)
	if [ $(cat ${OUTPUT_FOLDER}/11-Sequences/${subf}/${subf}.fasta | wc -l) -lt 1 ]; then
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No target was found for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
		printf "${subf}\tFailed\tNo target recovered from reference\n" >> progress.txt
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi

	## --------------------------------------------------
	## Recover pair-end reads for target (BWA & Samtools)
	## --------------------------------------------------

	# Create folder
	mkdir -p ${OUTPUT_FOLDER}/20-Alignment/${subf}

	if [ ${VERBOSE} -eq 2 ]; then printf "Alignment; "; fi

	# Alignment
	printf "### BWA mapping ###\n\n" > ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	bwa-mem2 index ${OUTPUT_FOLDER}/11-Sequences/${subf}/${subf}.fasta &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	bwa-mem2 mem ${OUTPUT_FOLDER}/11-Sequences/${subf}/${subf}.fasta ${INPUT_FOLDER}/${subf}/${subf}_R1.${FQEXT} ${INPUT_FOLDER}/${subf}/${subf}_R2.${FQEXT} -t ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.sam

	# Sam to BAM
	samtools view -b ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.sam -@ ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.bam
	# Sort BAM (Coordinate) for Variant Call
	# samtools sort -o ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.sort.bam -O bam ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.bam -@ ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	# Obtain BAM of mapped reads (properly pair)
	samtools view -b -q 30 -f 0x2 ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.sam 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.mapped.bam

	if [ ${VERBOSE} -eq 2 ]; then printf "Reads recovery; "; fi

	# Obtain Fastq's
	printf "\n\n### Reads recovery ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	samtools collate ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.mapped.bam ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.collate 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	samtools fastq -1 ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R1.fastq -2 ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R2.fastq -s ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_leftover.fastq ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.collate.bam 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	gzip -f ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_*.fastq

	## --------------------------------------------------
	## Rebuild target from reads (SPAdes, BWA & Samtools)
	## --------------------------------------------------

	# CHECK: Absent R1/R2 for de novo assembly
	if [[ $(zcat ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R1.fastq.gz | wc -l) -eq 0 || $( zcat ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R2.fastq.gz | wc -l) -eq 0 ]]; then
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No target reads were recovered for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
		printf "${subf}\tFailed\tNo reads recovered for target\n" >> progress.txt
		# Clean files/folders
		CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
		# Date
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi

	if [ ${VERBOSE} -eq 2 ]; then printf "Contig re-build; "; fi

	# SPAdes (Unambigous nucleotides assembly)
	printf "\n\n### Contig re-build ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	spades.py -1 ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R1.fastq.gz -2 ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} -o ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades -m ${MEM} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

	# Size select SPAdes recovered target
	if [ -f "${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.fasta" ]; then
		seqtk seq -L ${min_tl} ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.fasta > ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.seqtk.fasta
	else
		printf "\n${red}ERROR:${normal} No target was recovered for ${subf} by SPAdes. This could be caused by a low number of recovered reads.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		printf "${subf}\tFailed\tSPAdes did not assembly any target sequence.\n" >> progress.txt
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi
	# Check if sequences passed filtering
	if [ $(grep "^>" ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.seqtk.fasta | wc -l) -eq 0 ]; then
		printf "\n${red}ERROR:${normal} No target was recovered for ${subf} or the target recovered was too small. In the second case, make \"-l/--len-deviation\" larger. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		printf "${subf}\tSkipped\tRecovered target length below minimum length\n" >> progress.txt
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi
	
	# Define target without flanking regions
	mkdir -p ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking
	blastn -subject ${SEARCH_TARGET} -query ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.seqtk.fasta -outfmt "6 std sseq qlen" > ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.tsv

	# Size select blast hits
	while read -r line;do
		qstart=$(echo ${line} | cut -d " " -f 7)
		qend=$(echo ${line} | cut -d " " -f 8)
		total=$((qend-qstart))
		gaps=$(echo ${line} | cut -d " " -f 13 | grep -o "-" | wc -l)
		len=$((total-gaps))
		if [ $len -ge $min_tl ]; then
			printf "${line}\t${len}\n"
		fi
	done < ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.tsv > ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv
	
	# Check if hit was not recovered
	if [ $(cat ${OUTPUT_FOLDER}/20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | wc -l) == 0 ]; then
		printf "\n${red}ERROR:${normal} No target was recovered for ${subf} or the target recovered was too small. In the second case, make \"-l/--len-deviation\" larger. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		printf "${subf}\tSkipped\tRecovered target length below minimum length\n" >> progress.txt
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi

	# Select recovered sequence given minimum length
	cp ${OUTPUT_FOLDER}/20-Alignment/${subf}/spades/contigs.seqtk.fasta ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta

	# Align
	printf "\n\n### Contig re-alignment (BWA) ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	bwa-mem2 index ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	bwa-mem2 mem ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R1.fastq.gz ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sam
	
	samtools view -b ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sam -@ ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.bam
	samtools sort -o ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sort.bam -O bam ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.bam -@ ${THREADS} 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

	## -------------------------------------------
	## Variant Calling & Phasing (GATK & BCFTools)
	## -------------------------------------------

	# CHECK: Absent rebuilt BAM
	if [[ ! -f ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sort.bam || ! -f ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta ]]; then
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}}WARNING:${normal} No target reads were recovered for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
		printf "${subf}\tFailed\tNo reads recovered for target\n" >> progress.txt
		# Clean files/folders
		CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
		# Date
		date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		echo ""
		continue
	fi

	# Create folders
	mkdir -p ${OUTPUT_FOLDER}/30-VariantCalling/${subf}
	mkdir -p ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/mapped_filtered ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/genotyped ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/reference ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants

	if [ ${VERBOSE} -eq 2 ]; then printf "Variant calling; "; fi

	# Variant call
	printf "\n\n### Variant calling ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	variantCalling ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.rebuild.sort.bam ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.dict ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/mapped_filtered/${subf}.filtered.bam ${OUTPUT_FOLDER}/30-VariantCalling/${subf} ${THREADS} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

	# CHECK: No variants recovered
	if [[ -f ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) != 0 ]]; then
		# MULTIPLE HAPLOTYPES

		# Create folder
		mkdir -p ${OUTPUT_FOLDER}/40-Phasing/${subf}

		if [ ${VERBOSE} -eq 2 ]; then printf "Phasing; "; fi

		# Phasing
		printf "\n\n### Phasing ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		bash ${GENOME_PHASE} -s ${subf} -t ${THREADS} -r ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta -v ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz -i 20-Alignment -o ${OUTPUT_FOLDER}/40-Phasing/${subf} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# CHECK: Absent haplotypes
		if [[ ! -f ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_assembly_h1.fasta || ! -f ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_assembly_h2.fasta ]]; then
			if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No haplotypes were recovered for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
			printf "${subf}\tFailed\tNo haplotypes were found\n" >> progress.txt
			# Clean files/folders
			CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
			# Date
			date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			echo ""
			continue
		fi

		# Haplotype length correction: remove short matches
		mkdir -p ${OUTPUT_FOLDER}/50-haplotypes/${subf}
		seqtk seq -L ${min_tl} <(cat ${OUTPUT_FOLDER}/40-Phasing/${subf}/${subf}_assembly_h*.fasta) > ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta 

		# Rename headers
		hnum=$(cat ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta | sed 's/^> />/g' | grep "^>" | wc -l)

		for n in $(seq ${hnum})
		do
			if [ ${n} -eq "1" ]; then
				cat ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta | sed 's/_[0-9]_length_.*//g' | sed -z "s|NODE|seq_h${n}|${n}" > tmp
				mv tmp ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta
			else
				cat ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta | sed -z "s|NODE|seq_h${n}|1" > tmp
				mv tmp ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta
			fi
		done

		# Haplotype duplicate removal (SeqKit)
		seqkit rmdup -s ${OUTPUT_FOLDER}/50-haplotypes/${subf}/${subf}_haplotypes.fasta 2>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt > ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta

		## ------------------------------------
		## Haplotype abundance ratio (Kallisto)
		## ------------------------------------

		# CHECK: Absent clean haplotypes
		if [ ! -f ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta ]; then
			if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No haplotypes were recovered for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
			printf "${subf}\tSkipped\No haplotypes were recovered\n" > progress.txt
			date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			echo ""
			continue
		fi

		# Kallisto (Only applied to strains with more that one haplotype)
		if [ $(grep "^>" ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta | wc -l) -gt 1 ]
		then
			# Create folder
			mkdir -p ${OUTPUT_FOLDER}/60-Integration

			if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio; "; fi

			# Kallisto
			printf "\n\n### Kallisto ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			kallisto index -i ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			kallisto quant -i ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx -o ${OUTPUT_FOLDER}/60-Integration/${subf} ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# Filter haplotypes
			cp ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.orig.tsv
			Rscript ${FILTER_HAPLOTYPES} -i ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv -c ${CUTOFF} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# --------------- #
			# II. Copy Number #
			# --------------- #

			# Copy number
			copyNumber ${INPUT_FOLDER} ${subf}

			# ---------------- #
			# III. Integration #
			# ---------------- #

			# CHECK: Absent kallisto output
			if [ ! -f ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv ]; then
				if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing kallisto output for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
				printf "${subf}\tFailed\tKallisto could not determined the haplotype abundances\n" >> progress.txt
				# Clean files/folders
				CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
				# Date
				date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
				echo ""
				continue
			fi

			if [ ${VERBOSE} -eq 2 ]; then printf "Integration; "; fi
			# Log
			printf "\n\n### Integration ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# Integrate step I and II
			Rscript ${INTEGRATION} -r ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv -c ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv -i ${subf} -m "multiple" &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# ---------------- #
			# IV. Fingerprints #
			# ---------------- #

			# CHECK: Absent integration output
			if [ ! -f ${OUTPUT_FOLDER}/60-Integration/${subf}/integration.tsv ]; then
				if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing integration output for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
				printf "${subf}\tFailed\tIntegration could not be performed\n" >> progress.txt
				# Clean files/folders
				CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
				# Date
				date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
				echo ""
				continue
			fi

			# Fingerprint
			fingerPrint 'multiple' ${subf} 'Yes'
		# One or several variants but 1 haplotype
		elif [ $(grep "^>" ${OUTPUT_FOLDER}/50-haplotypes/${subf}/clean_${subf}_haplotypes.fasta | wc -l) == 1 ]; then
			
			if [ ${VERBOSE} -eq 2 ]; then printf "Phasing [Avoid]; "; fi
			if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio [Avoid]; "; fi

			# --------------- #
			# II. Copy Number #
			# --------------- #

			# Create folder
			mkdir -p ${OUTPUT_FOLDER}/60-Integration/${subf}

			# Copy number
			copyNumber ${INPUT_FOLDER} ${subf}

			# ---------------- #
			# III. Integration #
			# ---------------- #

			if [ ${VERBOSE} -eq 2 ]; then printf "Integration [unique]; "; fi
			# Log
			printf "\n\n### Integration ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# Create folder
			mkdir -p ${OUTPUT_FOLDER}/60-Integration/${subf}

			# Integrate only step II 
			Rscript ${INTEGRATION} -r "None" -c ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv -i ${subf} -m 'unique' &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

			# ---------------- #
			# IV. Fingerprints #
			# ---------------- #

			# Fingerprint
			fingerPrint 'unique' ${subf} 'None'

		fi
	# Check: No variants but multiple recovered targets
	elif [[ -f ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) == 0 && $(grep "^>" ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta | wc -l) -gt 1 ]]; then

		if [ ${VERBOSE} -eq 2 ]; then printf "Phasing [Avoid]; "; fi

		## ------------------------------------
		## Haplotype abundance ratio (Kallisto)
		## ------------------------------------

		# Kallisto (Applied to strains with two recovered targets - 2 haplotypes)
		# Create folder
		mkdir -p ${OUTPUT_FOLDER}/60-Integration

		if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio; "; fi

		# Kallisto
		printf "\n\n### Kallisto ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		kallisto index -i ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta.idx ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
		kallisto quant -i ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta.idx -o ${OUTPUT_FOLDER}/60-Integration/${subf} ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# Filter haplotypes
		cp ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.orig.tsv
		Rscript ${FILTER_HAPLOTYPES} -i ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv -c ${CUTOFF} &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# --------------- #
		# II. Copy Number #
		# --------------- #

		# Copy number
		copyNumber ${INPUT_FOLDER} ${subf}

		# ---------------- #
		# III. Integration #
		# ---------------- #

		# CHECK: Absent kallisto output
		if [ ! -f ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv ]; then
			if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing kallisto output for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
			printf "${subf}\tFailed\tKallisto could not determined the haplotype abundances\n" >> progress.txt
			# Clean files/folders
			CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
			# Date
			date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			echo ""
			continue
		fi

		if [ ${VERBOSE} -eq 2 ]; then printf "Integration; "; fi
		# Log
		printf "\n\n### Integration ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# Integrate step I and II
		Rscript ${INTEGRATION} -r ${OUTPUT_FOLDER}/60-Integration/${subf}/abundance.tsv -c ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv -i ${subf} -m "multiple" &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# ---------------- #
		# IV. Fingerprints #
		# ---------------- #

		# CHECK: Absent integration output
		if [ ! -f ${OUTPUT_FOLDER}/60-Integration/${subf}/integration.tsv ]; then
			if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing integration output for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
			printf "${subf}\tFailed\tIntegration could not be performed\n" >> progress.txt
			# Clean files/folders
			CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
			# Date
			date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
			echo ""
			continue
		fi

		# Fingerprint
		fingerPrint 'multiple' ${subf} "No"

	elif [[ -f ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat ${OUTPUT_FOLDER}/30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) == 0 && $(grep "^>" ${OUTPUT_FOLDER}/20-Alignment/${subf}/${subf}.fasta | wc -l) == 1 ]]; then
		# ONE HAPLOTYPE

		if [ ${VERBOSE} -eq 2 ]; then printf "Phasing [Avoid]; "; fi
		if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio [Avoid]; "; fi

		# --------------- #
		# II. Copy Number #
		# --------------- #

		# Create folder
		mkdir -p ${OUTPUT_FOLDER}/60-Integration/${subf}

		# Copy number
		copyNumber ${INPUT_FOLDER} ${subf}

		# ---------------- #
		# III. Integration #
		# ---------------- #

		if [ ${VERBOSE} -eq 2 ]; then printf "Integration [unique]; "; fi
		# Log
		printf "\n\n### Integration ###\n\n" >> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# Create folder
		mkdir -p ${OUTPUT_FOLDER}/60-Integration/${subf}

		# Integrate only step II 
		Rscript ${INTEGRATION} -r "None" -c ${OUTPUT_FOLDER}/60-Integration/${subf}/copy_number.tsv -i ${subf} -m 'unique' &>> ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt

		# ---------------- #
		# IV. Fingerprints #
		# ---------------- #

		# Fingerprint
		fingerPrint 'unique' ${subf} 'None'

	else
		# Absent variant file
		if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} VCF file missing for ${subf}. Computation will be skipped.\n" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt; fi
		continue
	fi

	if [ -f ${OUTPUT_FOLDER}/70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta ]; then
		printf "${subf}\tSuccess\tComputation finished\n" >> progress.txt
	fi

	# ----------------- #
	# ExI. Clean folder #
	# ----------------- #

	CleanFiles ${subf} ${INPUT_FOLDER} ${KEEPF}
	date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a ${OUTPUT_FOLDER}/01-Logsmain/log_${subf}.txt
	echo ""
done

function CreateSummary() {
	# Stats report
	if [[ ! -f Summary.tsv ]]; then

		# Calculations
		tl=$(grep -v "^>" ${SEARCH_TARGET} | tr -d "\n" | wc -c)

		# Header
		printf "#Isolate\tInput_folder\tTarget_file\tTarget_length\tRecovered_mean_target_length\tLength_deviation\tRecovered_reads\tNumber_SNPs\tCutoff\tNumber_haplotypes\tCopy_number\tHaplotype_ratio\tModified_output\n" > Summary.tsv

		# Loop Success
		if [[ -f progress.txt ]]; then
			for iso in $(grep "Success" progress.txt | cut -f 1); do
				# Variables
				rtl=$(grep -v "^>" ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}.fasta | awk '{ print length }' | awk '{ total += $1 } END { print total/NR }')
				recr1=$(zcat ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}_R1.fastq.gz | grep "^@" | wc -l)
				recr2=$(zcat ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}_R2.fastq.gz | grep "^@" | wc -l)
				nsnps=$(zcat ${OUTPUT_FOLDER}/30-VariantCalling/${iso}/variants/${iso}.vcf.gz | grep -v "#" | wc -l)
				nhaplo=$(grep "^>" ${OUTPUT_FOLDER}/50-haplotypes/${iso}/clean_${iso}_haplotypes.fasta | wc -l)
				cnum=$(cut -f 9 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | tail -n 1)
				rhaplo=$(cut -f 14 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | grep -v "per_haplotype" | tr "\n" "/" | sed 's/\/$//')
				mod=$(cut -f 15 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | tail -n 1)

				# Row
				if [[ -f ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv ]]; then
					printf "${iso}\t${INPUT_FOLDER}\t${SEARCH_TARGET}\t${tl}\t${rtl}\t${BPDEV}\t${recr1}/${recr2}\t${nsnps}\t${CUTOFF}\t${nhaplo}\t%.4f\t${rhaplo}\t${mod}\n" $cnum >> Summary.tsv
				fi
			done 
		fi
	else

		# Calculations
		tl=$(grep -v "^>" ${SEARCH_TARGET} | wc -c)

		# Loop Success
		if [[ -f progress.txt ]]; then
			for iso in $(grep "Success" progress.txt | cut -f 1); do
				# Variables
				rtl=$(grep -v "^>" ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}.fasta | awk '{ print length }' | awk '{ total += $1 } END { print total/NR }')
				recr1=$(zcat ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}_R1.fastq.gz | grep "^@" | wc -l)
				recr2=$(zcat ${OUTPUT_FOLDER}/20-Alignment/${iso}/${iso}_R2.fastq.gz | grep "^@" | wc -l)
				nsnps=$(zcat ${OUTPUT_FOLDER}/30-VariantCalling/${iso}/variants/${iso}.vcf.gz | grep -v "#" | wc -l)
				nhaplo=$(grep "^>" ${OUTPUT_FOLDER}/50-haplotypes/${iso}/clean_${iso}_haplotypes.fasta | wc -l)
				cnum=$(cut -f 9 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | tail -n 1)
				rhaplo=$(cut -f 14 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | grep -v "per_haplotype" | tr "\n" "/" | sed 's/\/$//')
				mod=$(cut -f 15 ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv | tail -n 1)

				# Strain not present
				if [[ -f ${OUTPUT_FOLDER}/60-Integration/${iso}/integration.tsv ]]; then
					if [[ $(grep -w -e ${iso} Summary.tsv | wc -l) -eq 0 ]]; then
						printf "${iso}\t${INPUT_FOLDER}\t${SEARCH_TARGET}\t${tl}\t${rtl}\t${BPDEV}\t${recr1}/${recr2}\t${nsnps}\t${CUTOFF}\t${nhaplo}\t%.4f\t${rhaplo}\t${mod}\n" $cnum >> Summary.tsv
					# Strain present in the same folder
					elif [[ $(grep -w -e ${iso} Summary.tsv | grep -e ${INPUT_FOLDER} | wc -l) -eq 1 ]]; then
						grep -vw ${iso} Summary.tsv > .tempfile
						printf "${iso}\t${INPUT_FOLDER}\t${SEARCH_TARGET}\t${tl}\t${rtl}\t${BPDEV}\t${recr1}/${recr2}\t${nsnps}\t${CUTOFF}\t${nhaplo}\t%.4f\t${rhaplo}\t${mod}\n" $cnum >> .tempfile
						mv .tempfile Summary.tsv
					fi
				fi
			done
		fi
	fi
}

CreateSummary &> /dev/null

# Final format

if [ ${VERBOSE} -ge 1 ]; then
	printf "Finish:\n"
	date "+    date: %d/%m/%Y" 
	date "+    time: %H:%M:%S"
	printf "\n"
fi
