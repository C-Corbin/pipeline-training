#Start in software folder
#Print start date and time for the console record
echo Pipeline run began at: $(date)
echo ''

#########################
#MAKE NAVIGATING EASIER #
#########################
#AN EXAMPLE:
#filepath="/home/stpuser/fake_folder"
#mkdir ${filepath}/nested_folder
#My shortcuts:
software_path="/software"
out_path="/home/stpuser"
seq_path="/example_fastqs"
ref_path="/reference_files"


#######################
# ASSIGNING VARIABLES #
#######################

###FASTQS and separating their names out###

fastqs=$1

IFS=',' read -r -a array <<< "$fastqs"
count="${#array[@]}"

fastq_L0001_R1="${array[0]}"
fastq_L0001_R2="${array[1]}"
fastq_L0002_R1="${array[2]}"
fastq_L0002_R2="${array[3]}"


###EVERY OTHER PARAMETER###
file_header=$2
#file header will be the name of the clinical sample, to which you can append 'aligned' or 'vcf' or whatevs
group=$3
library=$4
full_broad_panel_bed=$5
full_small_panel_bed=$6
req_depth=30
intron_depth=18
splice_interval=25
poly_list=${7}
year=${8}
user_initials=${9}
run_no=${10}
full_trans_file=${11}
dir=${12}
rescoping=${13}
exome=`echo ${14} | tr -d '\r'`

##Get broad_bed_name alone, separate from file path
IFS='/' read -r -a b_array <<< "$full_broad_panel_bed"
broad_bed_name="${b_array[-1]}"
broad_bed_name=${broad_bed_name%.*}
echo ${broad_bed_name}

##Get small_bed_name alone, separate from file path
IFS='/' read -r -a s_array <<< "$full_small_panel_bed"
small_bed_name="${s_array[-1]}"
small_bed_name=${small_bed_name%.*}
echo ${small_bed_name}


#################
#    FASTQC     #
#################

#Start anywhere as I'm encoding absolute paths
#Run Fastqc and write to own folder, output version number, date and times to log file
#cd ${software_path}/FastQC
#./fastqc -v 
#echo FastQC analysis start time: $(date)
#./fastqc $fastq_L0001_R1 --outdir=${out_path}/fastqcs
#./fastqc $fastq_L0001_R2 --outdir=${out_path}/fastqcs
#./fastqc $fastq_L0002_R1.gz --outdir=${out_path}/fastqcs
#./fastqc $fastq_L0002_R2.gz --outdir=${out_path}/fastqcs
#echo FastQC analysis end time: $(date)
#echo ''


#########################
# ALIGN AND MERGE LANES #
#########################
#Merging lanes after alignment conforms more to 'best practice'

#Run bwa mem on R1 R2 L001 L002, write to ../../home/stpuser/aligned_seqs
#The genome reference is in ../../reference_files/ucsc.hg19.nohap.masked.fasta

#echo bwa mem and samtools, align and merge lanes
#cd ${software_path}
#echo samtools version:
#samtools --version
#cd ${software_path}/bwa-0.7.15
#echo bwa mem version:
#bwa 2>&1 |grep Version | cut -f2 -d " "
#echo bwa mem and samtools start time: $(date)

#For sample given from arguments_for_pipeline.txt
#bwa mem ${ref_path}/ucsc.hg19.nohap.masked.fasta $fastq_L0001_R1 $fastq_L0001_R2 | samtools view -bh -o ${out_path}/aligned_seqs/temp_L001.bam
#bwa mem ${ref_path}/ucsc.hg19.nohap.masked.fasta $fastq_L0002_R1 $fastq_L0002_R2 | samtools view -bh -o ${out_path}/aligned_seqs/temp_L002.bam

##NB THE ABOVE IS TESTED AND FUNCTIONING WELL
##NEED TO WORK OUT HOW TO MERGE INTO SOMETHING WITH THE CORRECT NAME - SEE REAL PIPELINE

#Merge into one alignment file
#cd ${software_path}/samtools-1.3.1
#samtools merge ${out_path}/aligned_seqs/${file_header}.bam ${out_path}/aligned_seqs/temp_L001.bam ${out_path}/aligned_seqs/temp_L002.bam
#Delete temporary files
#rm ${out_path}/aligned_seqs/temp_L001.bam
#rm ${out_path}/aligned_seqs/temp_L002.bam
#


#########################
#  SORT AND INDEX SEQS  #
#########################

#cd ${software_path}
#echo Sort and index files
#echo samtools version:
#samtools --version
#cd ${software_path}/samtools-1.3.1
#echo samtools sort start time: $(date)
#Sort the files
#samtools sort -o ${out_path}/aligned_seqs/${file_header}_sorted.bam ${out_path}/aligned_seqs/${file_header}.bam
#echo samtools sort end time: $(date)
#echo ''

#echo samtools index start time: $(date)
#Index the files for IGV viewing
#cd ${software_path}/samtools-1.3.1
#samtools index ${out_path}/aligned_seqs/${file_header}_sorted.bam
#echo samtools index end time: $(date)
#echo ''


###################
# MARK DUPLICATES #
###################

#I need to mark duplicates for each file
#This also writes a metrics file for each flagged file
#cd ${software_path}/picard-tools-2.5.0
#echo Mark duplicates
#echo Picard version 2.5.0
#echo Picard MarkDuplicates start time: $(date)
#java -jar picard.jar MarkDuplicates I=${out_path}/aligned_seqs/${file_header}_sorted.bam O=${out_path}/aligned_seqs/${file_header}_sorted_dupflag.bam M=${out_path}/aligned_seqs/${file_header}_sorted_dupflag_metrics.txt
#echo Picard MarkDuplicates end time: $(date)
#echo ''


############
# FLAGSTAT #
############

#samtools flagstat for quality metrics - to see how they match up with the picard ones!
#cd ${software_path}/samtools-1.3.1
#echo Run flagstat quality metrics
#echo samtools version:
#samtools --version
#echo samtools flagstat start time: $(date)
#samtools flagstat ${out_path}/aligned_seqs/${file_header}_sorted_dupflag.bam > ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_flagstat.txt
#echo samtools flagstat end time: $(date)
#echo ''


###################
# ADD READ GROUPS #
###################
#In a real pipeline, read group information can be used to track how particular machines/runs/e.t.c are performing; they may also allow correct calibration?
#However in our very-few-samples example we can probably put in placeholder read groups
#Adding them quite late - 'real' ones would be added early, while still at the lanes-and-reads stage

#cd ${software_path}/picard-tools-2.5.0
#echo Add read groups
#echo Picard version 2.5.0
#echo Picard AddOrReplaceReadGroups start time: $(date)

#java -jar picard.jar AddOrReplaceReadGroups \
	I=${out_path}/aligned_seqs/${file_header}_sorted_dupflag.bam \
	O=${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=${file_header}

#echo Picard AddOrReplaceReadGroups end time: $(date)
#echo ''


#####################
# INDEL REALIGNMENT #
#####################

#cd ${software_path}/samtools-1.3.1
#echo Indel realignment: index files with samtools index
#echo samtools version:
#samtools --version
#echo samtools index start time: $(date)

###INDEXING - needed for indel-realignment files
#samtools index ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.bam

#echo samtools index end time: $(date)
#echo ''

#cd ${software_path}/GenomeAnalysisTK-3.6
#echo Indel realignment: identify target regions with GATK RealignerTargetCreator
#echo GenomeAnalysisTK version:
#java -jar GenomeAnalysisTK.jar -version
#echo GenomeAnalysisTK RealignerTargetCreator start time: $(date)

### IDENTIFY TARGET REGIONS
### CAN RESTRICT TO THE BROAD BED FILE WITH -L, the GATK 'intervals' option
#
#java -jar GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${ref_path}/ucsc.hg19.nohap.masked.fasta \
	-L ${full_broad_panel_bed} \
	-I ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.bam \
	-o ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.intervals

#echo GenomeAnalysisTK RealignerTargetCreator end time: $(date)
#echo ''

#cd ${software_path}/GenomeAnalysisTK-3.6
#echo Indel realignment: perform the realignment with GATK IndelRealigner
#echo GenomeAnalysisTK version:
#java -jar GenomeAnalysisTK.jar -version
#echo GenomeAnalysisTK IndelRealigner start time: $(date)

##  INDEL REALIGNMENT Test run
#java -jar GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${ref_path}/ucsc.hg19.nohap.masked.fasta \
	-L ${full_broad_panel_bed} \
	-I ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.bam \
	-targetIntervals ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG.intervals \
	-o ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG_realigned.bam

#echo GenomeAnalysisTK IndelRealigner end time: $(date)
#echo ''


##########################################################
# EXON PLUS/MINUS 25bp COVERAGE; ON AND OFF TARGET READS #
##########################################################

## Picard BedToIntervalList
#cd ${software_path}/picard-tools-2.5.0
#java -jar picard.jar BedToIntervalList \
#	I=${full_broad_panel_bed} \
#	O=${out_path}/bed_intervals/${broad_bed_name}.intervals \
#	SD=${ref_path}/ucsc.hg19.nohap.masked.dict
#java -jar picard.jar BedToIntervalList \
#	I=${full_small_panel_bed} \
#	O=${out_path}/bed_intervals/${small_bed_name}.intervals \
#	SD=${ref_path}/ucsc.hg19.nohap.masked.dict


#Using master bed as the bait region for now. Test run.
#cd ${software_path}/picard-tools-2.5.0
##  Picard HsMetrics tools for 30X coverage
#java -jar picard.jar CollectHsMetrics \
#	I=${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG_realigned.bam \
#	o=${out_path}/aligned_seqs/${file_header}_hsmetrics.txt \
#	R=${ref_path}/ucsc.hg19.nohap.masked.fasta \
#	BAIT_INTERVALS=${out_path}/bed_intervals/${broad_bed_name}.intervals \
#	TARGET_INTERVALS=${out_path}/bed_intervals/${small_bed_name}.intervals

###THE ABOVE ALL WORKS OKAY###
	
## Use Sambamba for exonic coverage
## DOESN'T WORK YET - CAN'T FIND DOCUMENTATION

cd ${software_path}/sambamba_v0.6.3
sambamba depth region -L ${full_broad_panel_bed} -T 30 -m -F “mapping_quality >= 30” ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG_realigned.bam > ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG_realigned.sambamba_output.bed
#

###############################################
# VARIANT CALLING, DECOMPOSING AND ANNOTATING #
###############################################
#
#Call variants with GATK haplotype caller
#cd ${software_path}/GenomeAnalysisTK-3.6
#java -jar GenomeAnalysisTK.jar \
#	-R ${ref_path}/ucsc.hg19.nohap.masked.fasta \
#	-T HaplotypeCaller \
#	-I ${out_path}/aligned_seqs/${file_header}_sorted_dupflag_RG_realigned.bam \
#	-L ${full_small_panel_bed} \
#	-stand_call_conf 30 \
#	--output_mode EMIT_VARIANTS_ONLY \
#	-maxAltAlleles 10 \
#	-o ${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants.vcf

#Decompose variants with Vt
#cd ${software_path}/vt
#vt decompose ${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants.vcf -o ${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp.vcf

#Annotate the decomposed variants. Restrict to small bed.
#SnpEff, human genome 37
#Use following command to look at human databases: java -jar snpEff.jar databases | less | grep "Homo*" 
#cd ${software_path}/snpEff
#java -Xmx4g -jar snpEff.jar \
#-fi ${full_small_panel_bed} \
#GRCh37.75 \
#${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp.vcf > \
#${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp_ann.vcf
#
#SnpSift for ExAc
#cd ${software_path}/snpEff
#java -jar SnpSift.jar annotate -v ${ref_path}/exac/ExAC.r0.3.1.sites.vep.vcf.gz \
# ${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp.vcf > \
#${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp_exac.vcf

#After: filter variants above 5%
#cd ${software_path}/snpEff
#cat ${out_path}/vcf/${file_header}_sorted_dupflag_RG_realigned_variants_decomp_exac.vcf | java -jar #SnpSift.jar filter "(#filtered thing)" | > \ 
#${out_path}/vcf/filtered.vcf


#Print end date and time for the console record
echo '' 
echo Pipeline run completed at: $(date)