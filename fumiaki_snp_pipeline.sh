#!/usr/bin/bash

Help()
{
   # Display Help
   echo "This script runs Fumiaki's SNP pipeline."
   echo "Syntax: ${0} [-r|-R|-p|-o|-m|-g|-P|-w|-h|-G]"
   echo "options:"
   echo "-r     Full path to reference fasta file."
   echo "-R     Path to directory with fastq reads."
   echo "-p     Prefix file, one prefix per line."
   echo "-o     output prefix."
   echo "-w     path to working directory [default: current dir]."
   echo "-m     perform mapping flag (default: skip mapping)."
   echo "-g     Run only gvcf generation."
   echo "-G     input Prefix (only with -g flag)."
   echo "-P     Ploidy [int=2]."
   echo "-h     show this help."
}

while getopts "gmhr:R:p:o:w:P:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        r) REFERENCE=${OPTARG};;
        R) READS=${OPTARG};;
        p) PREFIX_FILE=${OPTARG};;
        o) OUTPUT=${OPTARG};;
		m) DO_MAPPING=true;;
		g) GVCF_ONLY=true;;
		w) WORKDIR=${OPTARG};;
		G) PREFIX=${OPTARG};;
		P) PLOIDY=${OPTARG};;
        \?) # incorrect option
                echo
                echo "Error, Invalid option"
                echo
                Help
                exit;;
   esac
done

set -o errexit

#!/bin/sh

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
if [ ! -d "UGE-output" ] 
then #Create output directory in case it does NOT exist
     mkdir UGE-output
fi


####Have to update version before run

module load bwa/0.7.17
module load GATK/4.3.0.0 #GATK/3.7-0-Java-1.8.0_92  
module load picard/2.27.3 #picard/2.22.7-Java-1.8.0_92
module load samtools/1.16.1 #SAMtools/1.3-goolf-1.7.20-HTSlib-1.3

# These are all of the module paths to the pieces of GATK that I use to map and quality control the genome
# They're all variables so that I can easily update the script if the paths in Locus changes
CD="${PICARDJARPATH}/picard.jar CreateSequenceDictionary"
PC="${PICARDJARPATH}/picard.jar SortSam"
MD="${PICARDJARPATH}/picard.jar MarkDuplicates"
MR="${PICARDJARPATH}/picard.jar MergeSamFiles"
BI="${PICARDJARPATH}/picard.jar BuildBamIndex"
#GA="$EBROOTGATK/GenomeAnalysisTK.jar"

# Path to bowtie reference genome database
BWA_INDEXES="/gpfs/gsfs12/users/lorenziha/DATABASES/mm10/"
# mm10.fa.sa

if [ -z ${PLOIDY+y} ]
then
        # -g flag off
        PLOIDY=2
fi

# These are all of the variables you should need to change in a standard mapping BEFORE running the script
# Tgindex is the Toxo reference that I use to map against. I have multiple so I comment out all except the one I'm currently using
# tempdir is the folder with a lot of available space I leave open so java has enough space for all the temporary files it makes while doing calculations
# workingdir is the folder set that I want to do my work in and generate output to
# fastqdir is the folder where my fastq files live I do this so I don't always have to have them in the same folder
# strain_text_file is a text file with each file name created from your sequencing without the extension 'sam'.  The text file should have one ID per line.  This file will indicate with files to run through the processing and variant calling pipeline.

if [ -z ${WORKDIR+y} ]
then
        # -g flag off
        WORKDIR=`pwd -P`
fi

Tgindex=${REFERENCE}
db_index=`basename ${REFERENCE}`
#Tgindex_out=/hpcdata/lpd/Fumiaki/Toxo_genomes/TgME49/ToxoDB-51_TgondiiME49_Genome
workingdir=${WORKDIR}
tempdir=/lscratch/${SLURM_JOB_ID} #${workingdir}/TMP
fastqdir=${READS}
strain_text_file=${PREFIX_FILE}

# Check if workdir exists
if [ ! -d ${workingdir} ]; then
    mkdir ${workingdir}
fi        

if [ ${GVCF_ONLY} ]; then
	if [ ${PREFIX} ]; then
		echo ${PREFIX} > ${PREFIX}.prefix
		strain_text_file=${PREFIX}.prefix	
	else
		echo ERROR, you must specify option -P with flag -g; echo
		exit 1
	fi
fi


## Step 0
## create index of genome for bwa and GATK
if [ $DO_MAPPING ]; then
	if [ ! -f ${BWA_INDEXES}/${db_index}.sa ]; then
		echo; echo "Crating bwa index files for ${db_index}"; echo
		bwa index ${REFERENCE} 
		echo "Done!!"
		echo
	fi

    if [ ! -f ${BWA_INDEXES}/${REFERENCE}.fai ]; then
        echo; echo "Crating genome index for ${REFERENCE}"; echo
        samtools faidx ${REFERENCE}
        echo "Done!!"
        echo
    fi        

	## Check for GATK dict file for reference
	DICT=${BWA_INDEXES}/${db_index}.dict
	if [ ! -f ${DICT} ]; then
		echo; echo "Creating dict file for referece ${Tgindex}"; echo
		java -Xmx1G -jar $CD R=${REFERENCE} O=${BWA_INDEXES}/${db_index}.dict ##picard CCreateSequenceDictionary
		echo "Done!"
		echo
	fi
fi

## Step 1
## map to ME49
## create a text file with each file name without the extension 'sam'.  The text file should have one ID per line.  This file will indicate with files to run through the processing and variant calling pipeline.
## You will need to ensure that the fastq files that go into your mapping are named such that this code can find them using the fastqdir and line names given in your text file combined with the extension before .fastq
cd $workingdir/

while read prefix 
do
	prefix_list+=(${prefix}) # stores prefix names for later
	if [ $DO_MAPPING && ! -f ${prefix}.bam]; then
		#echo; echo "Running bowtie2 mem on ${prefix}_R1.fastq.gz and ${prefix}_R2.fastq.gz"; echo 
		#bowtie2 -x ${db_index} \
		#--end-to-end \
		#--rg-id "${prefix}" \
		#-rg "SM:${prefix}" --rg "PL:illumina" --rg "LB:lib1" --rg "PU:unit1" \
		#--threads 16 \
		#-1 $fastqdir/${prefix}_R1.fastq.gz \
		#-2 $fastqdir/${prefix}_R2.fastq.gz -b $workingdir/${prefix}.bam
		#echo ${prefix}_BOWTIE2

		echo; echo "Running bwa mem on ${prefix}.R1.fastq.gz and ${prefix}.R2.fastq.gz"; echo 
		echo bwa mem -R "@RG\tID:${prefix}\tSM:${prefix}\tPL:illumina\tLB:lib1\tPU:unit1" -t 16 -M ${REFERENCE} $fastqdir/${prefix}.R1.fastq.gz $fastqdir/${prefix}.R2.fastq.gz
		bwa mem -R "@RG\tID:${prefix}\tSM:${prefix}\tPL:illumina\tLB:lib1\tPU:unit1" \
		-t 16 \
		-M ${REFERENCE} $fastqdir/${prefix}.R1.fastq.gz $fastqdir/${prefix}.R2.fastq.gz | samtools view -hb  > ${prefix}.bam
		echo ${prefix}_BWA
	fi
done < ${strain_text_file}

## Step 2
## run pipeline to process bam file and call variants. 
## the parameters for calling variants are basically the default ones recommended by gatk for the current version 3.1.1

#####Caution!!! HAVE TO CHECK RECOMMENDED PARAMETERS FOR THE CURRRENT CERSION

## the pipeline runs in parallel in the HaplotypeCaller step, everything else is not being parallelized.  My recommendation to run this script for large number of samples would be to submit various jobs.  
## The last step includes a ploidy argument. Change it if your genome is not haploid (1)
while read prefix
do
	####################################################################################
	# Analysis-ready reads
	####################################################################################

	echo "Calling variants (gvcf) for  ${prefix} ...." 
	if [ ! -f ${prefix}.sorted.bam ]; then
		echo; echo "Running ${PC} on ${prefix}.bam"; echo 
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $PC INPUT=${prefix}.bam OUTPUT=${prefix}.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
	fi 
	if [ ! -f ${prefix}.dedup.bam ]; then
		echo; echo "Running ${MD} on ${prefix}.sorted.bam"; echo
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $MD INPUT=${prefix}.sorted.bam OUTPUT=${prefix}.dedup.bam M=${prefix}.metrics.txt
	fi
	if [ ! -f ${prefix}.dedup.bai ]; then
		echo; echo "Running ${BI} on ${prefix}.sorted.bam"; echo	
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $BI INPUT=${prefix}.dedup.bam
	fi

	##########################################
	if [ ! -f ${prefix}.BaseRecalibrator.table ]; then
		echo; echo "Running gatk -T BaseRecalibrator on ${prefix}.dedup.bam"; echo
		gatk --java-options "-Xmx30G -Djava.io.tmpdir=$tempdir" BaseRecalibrator \
			-R ${REFERENCE} \
            -I ${prefix}.dedup.bam \
            --known-sites /fdb/GATK_resource_bundle/mm10/dbsnp146_fixedNames.vcf.gz \
            -O ${prefix}.BaseRecalibrator.table
	fi

	if [ ! -f ${prefix}.recalibrator.bam ]; then
		echo; echo "Running gatk -T ApplyBQSR on ${prefix}.dedup.bam"; echo
		gatk --java-options "-Xmx30G -Djava.io.tmpdir=$tempdir" ApplyBQSR \
			-R ${REFERENCE} \
			-I ${prefix}.dedup.bam \
			--bqsr-recal-file ${prefix}.BaseRecalibrator.table \
			-O ${prefix}.recalibrator.bam
	fi
	##########################################


	####################################################################################
	# Call variants per sample in GVCF mode
	####################################################################################
	
	if [ ! -f ${prefix}.raw.snps.indels.g.vcf ]; then
		echo; echo "Running gatk -T HaplotypeCaller -ERC GVCF>> on ${prefix}.recalibrator.bam"; echo
		gatk --java-options "-Xmx100G -Djava.io.tmpdir=${tempdir}" \
			HaplotypeCaller \
			-ERC GVCF \
			-R ${REFERENCE} \
			-I ${prefix}.recalibrator.bam \
			-O ${prefix}.raw.snps.indels.g.vcf \
			--sample-ploidy ${PLOIDY} 
	fi
 	echo Done running ${prefix}!
 done < ${strain_text_file}


if [ $GVCF_ONLY ]; then
	echo; echo Done running!!; echo
	echo Leaving pipeline !!
	exit 0
fi

	####################################################################################
	# Consolidate GVCFs and joint calls
	####################################################################################

## Step 3
## Combine the GVCF files (created in the last part of the last script into chunks small enough for the Cluster to handle with the available memory
## These will all be generated by the script and named .raw.snps.indels.g.vcf
## Modify the naming so that it makes sense to you for the -o output files
# cd $workingdir

echo; echo Starting Step3; echo

# Loop  => Number of groups to be processed together = 8
min_prefix_per_group=5
num_of_prefix=${#prefix_list[@]}
min_groups=`expr ${num_of_prefix} \/ ${min_prefix_per_group} + 1`
merge_param=()

for (( loop1=0; loop1 < min_groups; loop1++ )); do
	
	for (( j=0; j < min_prefix_per_group; j++ )); do
		if [[ `expr $loop1 + $j` > 0 ]]; then
			loop2=`expr $loop1 + $j` # `expr ${loop1} \* ${min_prefix_per_group} + ${j}`
		else
			loop2=0
		fi
		
		if [ ${prefix_list[${loop2}]} ]; then
			variant_parameter="${variant_parameter} --variant ${prefix_list[${loop2}]}.raw.snps.indels.g.vcf"
		fi
	done
	if [ ! -f CombineGVCFs_${loop1}.g.vcf ]; then
		echo; echo "Running CombineGVCFs in loop number ${loop1}"; echo
		gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
			CombineGVCFs \
 			-R ${REFERENCE} \
 			-O CombineGVCFs_${loop1}.g.vcf \
			${variant_parameter}
		echo "Done with ${loop1}!"
	fi
	merge_param+="--variant CombineGVCFs_${loop1}.g.vcf "
	variant_parameter=''
done
echo GVCF consoliodation done.; echo


# ## Step 4
# # Merge all desired g.vcf files from CombineGVCFs into a single VCF file of all strains
# # You will need to change the output of this name as well as making sure the input matches what you created in Step 3
# # This can run for a VERY long time. 96 Toxo genomes ran for 12 days straight before finishing
#cd $workingdir
if [ ! -f ${OUTPUT}_raw_variants.vcf ]; then
echo; echo "Running GenotypeGVCFs on *.g.vcf files" ; echo
gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir" \
	GenotypeGVCFs \
	-R ${REFERENCE} \
	--max-alternate-alleles 6 \
	--annotate-with-num-discovered-alleles \
	-standard-min-confidence-threshold-for-calling 30 \
	-O ${OUTPUT}_raw_variants.vcf \
	${merge_param}
fi

## Step 4 - Hard filtering

## Note: For this part , Fumiaki told me to follow the site: https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/2806-how-to-apply-hard-filters-to-a-call-set

##I did not use this part but should be:

##1. Extract the SNPs from the call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_raw_snps.vcf ]; then
echo; echo "Running SelectVariants on ${OUTPUT}_raw_variants.vcf" ; echo
gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir" \
	SelectVariants \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_variants.vcf \
	--select-type-to-include SNP \
	-O ${OUTPUT}_raw_snps.vcf
fi

##This creates a VCF file called raw_snps.vcf, containing just the SNPs from the original file of raw variants.


##2. Apply the filter to the SNP call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_filtered_snps.vcf ]; then
echo; echo "Running VariantFiltration on ${OUTPUT}_raw_snps.vcf" ; echo
gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir" \
	VariantFiltration \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_snps.vcf \
	--filter-name "QDlow" \
	--filter-expression "QD < 2.0" \
	--filter-name "FShigh" \
	--filter-expression "FS > 60.0" \
	--filter-name "MQlow" \
	--filter-expression "MQ < 40.0" \
	--filter-name "MQRSlow" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "RPRSlow" \
	--filter-expression "ReadPosRankSum < -8.0" \
	-O ${OUTPUT}_filtered_snps.vcf
fi

##3.Extract the Indels from the call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_raw_indels.vcf ]; then
echo; echo "Running SelectVariants on ${OUTPUT}_raw_variants.vcf" ; echo
gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir" \
	SelectVariants \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_variants.vcf \
	--select-type-to-include INDEL \
	-O ${OUTPUT}_raw_indels.vcf
fi

##6. Apply the filter to the Indel call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_filtered_indels.vcf ]; then
echo; echo "Running VariantFiltration on ${OUTPUT}_raw_indels.vcf" ; echo
gatk --java-options "-Xmx13G -Djava.io.tmpdir=$tempdir" \
	VariantFiltration \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_indels.vcf \
	--filter-name "QDlow" \
	--filter-expression "QD < 2.0" \
	--filter-name "FShigh" \
	--filter-expression "FS > 200.0" \
	--filter-name "RPRSlow" \
	--filter-expression "ReadPosRankSum < -20.0" \
	-O ${OUTPUT}_filtered_indels.vcf
fi
