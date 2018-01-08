#Start in software folder
#Print start date and time for the console record
echo Pipeline run began at: $(date)
echo ''


#################
#    FASTQC     #
#################

#Start in software folder
#Run Fastqc and write to own folder, output version number, date and times to log file
#cd FastQC
#./fastqc -v 
#echo FastQC analysis start time: $(date)
#./fastqc ../../example_fastqs/*fastq.gz --outdir=/home/stpuser/fastqcs
#echo FastQC analysis end time: $(date)
#echo ''


#########################
# ALIGN AND MERGE LANES #
#########################
#Merging lanes after alignment conforms more to 'best practice'

#Run bwa mem on R1 R2 L001 L002 for the samples, write to ../../home/stpuser/aligned_seqs
#The genome reference is in ../../reference_files/ucsc.hg19.nohap.masked.fasta

#echo bwa mem and samtools, align and merge lanes
#cd ..
#echo samtools version:
#samtools --version
#cd bwa-0.7.15
#echo bwa mem version:
#bwa 2>&1 |grep Version | cut -f2 -d " "
#echo bwa mem and samtools start time: $(date)

#Sample 1504850-S1509352-02_GCTCGGTA
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1504850-S1509352-02_GCTCGGTA_L001_R1_001.fastq.gz ../../example_fastqs/1504850-S1509352-02_GCTCGGTA_L001_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L001.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1504850-S1509352-02_GCTCGGTA_L002_R1_001.fastq.gz ../../example_fastqs/1504850-S1509352-02_GCTCGGTA_L002_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L002.bam
#Merge into one alignment file
#cd ../samtools-1.3.1
#samtools merge ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam ../../home/stpuser/aligned_seqs/temp_L001.bam ../../home/stpuser/aligned_seqs/temp_L002.bam
#Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
#
#
##Sample 1606034-S1612259-02_CGAACTTA
#cd ../software/bwa-0.7.15
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1606034-S1612259-02_CGAACTTA_L001_R1_001.fastq.gz ../../example_fastqs/1606034-S1612259-02_CGAACTTA_L001_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L001.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1606034-S1612259-02_CGAACTTA_L002_R1_001.fastq.gz ../../example_fastqs/1606034-S1612259-02_CGAACTTA_L002_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L002.bam
#Merge into one alignment file
#cd ../samtools-1.3.1
#samtools merge ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam ../../home/stpuser/aligned_seqs/temp_L001.bam ../../home/stpuser/aligned_seqs/temp_L002.bam
#Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
#
#
#Sample 1607686-S1615531-02_TTCACGCA
#cd ../software/bwa-0.7.15
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1607686-S1615531-02_TTCACGCA_L001_R1_001.fastq.gz ../../example_fastqs/1607686-S1615531-02_TTCACGCA_L001_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L001.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1607686-S1615531-02_TTCACGCA_L002_R1_001.fastq.gz ../../example_fastqs/1607686-S1615531-02_TTCACGCA_L002_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L002.bam
##Merge into one alignment file
#cd ../samtools-1.3.1
#samtools merge ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam ../../home/stpuser/aligned_seqs/temp_L001.bam ../../home/stpuser/aligned_seqs/temp_L002.bam
##Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
##Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
#
##Sample 1609778-S1620040-02_CGCTGATC
#cd ../software/bwa-0.7.15
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1609778-S1620040-02_CGCTGATC_L001_R1_001.fastq.gz ../../example_fastqs/1609778-S1620040-02_CGCTGATC_L001_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L001.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1609778-S1620040-02_CGCTGATC_L002_R1_001.fastq.gz ../../example_fastqs/1609778-S1620040-02_CGCTGATC_L002_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L002.bam
##Merge into one alignment file
#cd ../samtools-1.3.1
#samtools merge ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam ../../home/stpuser/aligned_seqs/temp_L001.bam ../../home/stpuser/aligned_seqs/temp_L002.bam
##Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
#
##Sample 1703057-S1705957-02_AAGACGGA
#cd ../software/bwa-0.7.15
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1703057-S1705957-02_AAGACGGA_L001_R1_001.fastq.gz ../../example_fastqs/1703057-S1705957-02_AAGACGGA_L001_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L001.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../example_fastqs/1703057-S1705957-02_AAGACGGA_L002_R1_001.fastq.gz ../../example_fastqs/1703057-S1705957-02_AAGACGGA_L002_R2_001.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/temp_L002.bam
##Merge into one alignment file
#cd ../samtools-1.3.1
#samtools merge ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam ../../home/stpuser/aligned_seqs/temp_L001.bam ../../home/stpuser/aligned_seqs/temp_L002.bam
##Delete temporary files
#rm ../../home/stpuser/aligned_seqs/temp_L001.bam
#rm ../../home/stpuser/aligned_seqs/temp_L002.bam
#echo bwa mem and samtools end time: $(date)
#echo ''


#########################
#  SORT AND INDEX SEQS  #
#########################

#cd ..
#echo Sort and index files
#echo samtools version:
#samtools --version
#cd samtools-1.3.1
#echo samtools sort start time: $(date)
#Sort the files
#samtools sort -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam
#echo samtools sort end time: $(date)
#echo ''

#echo samtools index start time: $(date)
#Index the files for IGV viewing
#samtools index ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam
#echo samtools index end time: $(date)
#echo ''


###################
# MARK DUPLICATES #
###################

#echo Mark duplicates
#echo Picard version 2.5.0
#echo Picard MarkDuplicates start time: $(date)
#I need to mark duplicates for each file
#This also writes a metrics file for each flagged file
#cd ../picard-tools-2.5.0 
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam O=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam O=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam O=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam O=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam O=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_metrics.txt
#echo Picard MarkDuplicates end time: $(date)
#echo ''


############
# FLAGSTAT #
############

#echo Run flagstat quality metrics
#echo samtools version:
#samtools --version
#echo samtools flagstat start time: $(date)
#samtools flagstat for quality metrics - to see how they match up with the picard ones!
#cd ../samtools-1.3.1
#samtools flagstat ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_flagstat.txt
#echo samtools flagstat end time: $(date)
#echo ''


###################
# ADD READ GROUPS #
###################
#In a real pipeline, read group information can be used to track how particular machines/runs/e.t.c are performing; they may also allow correct calibration?
#However in our 4-sample example we can probably put in placeholder read groups
#Adding them quite late - 'real' ones would be added early, while still at the lanes-and-reads stage

#echo Add read groups
#echo Picard version 2.5.0
#echo Picard AddOrReplaceReadGroups start time: $(date)
#1504850-S1509352-02_GCTCGGTA
#cd ../picard-tools-2.5.0
#java -jar picard.jar AddOrReplaceReadGroups \
	I=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam \
	O=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=1504850-S1509352-02_GCTCGGTA
#
#1606034-S1612259-02_CGAACTTA
#java -jar picard.jar AddOrReplaceReadGroups \
	I=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam \
	O=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=1606034-S1612259-02_CGAACTTA
#
#1607686-S1615531-02_TTCACGCA
#java -jar picard.jar AddOrReplaceReadGroups \
	I=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam \
	O=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=1607686-S1615531-02_TTCACGCA
#
#1609778-S1620040-02_CGCTGATC
#java -jar picard.jar AddOrReplaceReadGroups \
	I=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam \
	O=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=1609778-S1620040-02_CGCTGATC
#
#1703057-S1705957-02_AAGACGGA
#java -jar picard.jar AddOrReplaceReadGroups \
	I=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam \
	O=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_RG.bam \
	RGID=1 \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=1703057-S1705957-02_AAGACGGA
	
#echo Picard AddOrReplaceReadGroups end time: $(date)
#echo ''


#####################
# INDEL REALIGNMENT #
#####################

#cd ../samtools-1.3.1
#echo Indel realignment: index files with samtools index
#echo samtools version:
#samtools --version
#echo samtools index start time: $(date)

###INDEXING - needed for indel-realignment files
#samtools index ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_RG.bam

#echo samtools index end time: $(date)
#echo ''

#cd ../GenomeAnalysisTK-3.6
#echo Indel realignment: identify target regions with GATK RealignerTargetCreator
#echo GenomeAnalysisTK version:
#java -jar GenomeAnalysisTK.jar -version
#echo GenomeAnalysisTK RealignerTargetCreator start time: $(date)

### IDENTIFY TARGET REGIONS
### CAN RESTRICT TO THE BROAD BED FILE WITH -L, the GATK 'intervals' option
#
#1504850-S1509352-02_GCTCGGTA - NGD Complete

#java -jar GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-I ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam \
#	-o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.intervals
#
#1606034-S1612259-02_CGAACTTA - CTD
#java -jar GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/CTDFinaldesignwith25bp_v3.bed \
#	-I ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_RG.bam \
#	-o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_RG.intervals
#
#1607686-S1615531-02_TTCACGCA - Motor Complete
#java -jar GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/Motor_CompletePanel_v1.bed \
#	-I ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_RG.bam \
#	-o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_RG.intervals
#
#1609778-S1620040-02_CGCTGATC - IEM All Panels
#java -jar GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/IEM_all_panels_header.bed \
#	-I ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_RG.bam \
#	-o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_RG.intervals
#
#1703057-S1705957-02_AAGACGGA - HeredCancer full
#java -jar GenomeAnalysisTK.jar \
#	-T RealignerTargetCreator \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/HeredCancer_full_panel_25bp_v1.bed \
#	-I ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_RG.bam \
#	-o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_RG.intervals

#echo GenomeAnalysisTK RealignerTargetCreator end time: $(date)
#echo ''

#cd ../GenomeAnalysisTK-3.6
#echo Indel realignment: perform the realignment with GATK IndelRealigner
#echo GenomeAnalysisTK version:
#java -jar GenomeAnalysisTK.jar -version
#echo GenomeAnalysisTK IndelRealigner start time: $(date)

##  INDEL REALIGNMENT Test run
#cd ../GenomeAnalysisTK-3.6
#java -jar GenomeAnalysisTK.jar \
#	-T IndelRealigner \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/NGD_CompletePanel_25bp_v2.bed \
#	-I ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam \
#	-targetIntervals ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.intervals \
#	-o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.bam

#echo GenomeAnalysisTK IndelRealigner end time: $(date)
#echo ''


##########################################################
# EXON PLUS/MINUS 25bp COVERAGE; ON AND OFF TARGET READS #
##########################################################

## Picard BedToIntervalList
#1504850-S1509352-02_GCTCGGTA - required bed files
#cd ../picard-tools-2.5.0
#java -jar picard.jar BedToIntervalList \
#	I=../../example_fastqs/beds/NGD_CompletePanel_25bp_v2.bed \
#	O=../../home/stpuser/bed_intervals/NGD_CompletePanel_25bp_v2.intervals \
#	SD=../../reference_files/ucsc.hg19.nohap.masked.dict
#java -jar picard.jar BedToIntervalList \
#	I=../../example_fastqs/beds/NGD_dystonia_v3_25bp.bed \
#	O=../../home/stpuser/bed_intervals/NGD_dystonia_v3_25bp.intervals \
#	SD=../../reference_files/ucsc.hg19.nohap.masked.dict
#

#Using master bed as the bait region for now. Test run.
##  Picard HsMetrics tools for 30X coverage
#1504850-S1509352-02_GCTCGGTA
#cd ../picard-tools-2.5.0
#java -jar picard.jar CollectHsMetrics \
#	I=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.bam \
#	o=../../home/stpuser/aligned_seqs/1504850-S1509352-02_hsmetrics.txt \
#	R=../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	BAIT_INTERVALS=../../home/stpuser/bed_intervals/NGD_CompletePanel_25bp_v2.intervals \
#	TARGET_INTERVALS=../../home/stpuser/bed_intervals/NGD_dystonia_v3_25bp.intervals
#
	
## Use Sambamba for exonic coverage
## DOESN'T WORK YET

#1504850-S1509352-02_GCTCGGTA
#cd ../sambamba_v0.6.3
#sambamba depth region -L ../../example_fastqs/beds/NGD_dystonia_v3_25bp.bed -T 30 -m -F “mapping_quality >= 30” ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.bam > ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.sambamba_output.bed
#

###############################################
# VARIANT CALLING, DECOMPOSING AND ANNOTATING #
###############################################
#
#Call variants with GATK haplotype caller
#cd ../GenomeAnalysisTK-3.6
#java -jar GenomeAnalysisTK.jar \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-T HaplotypeCaller \
#	-I ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.bam \
#	-L /example_fastqs/beds/NGD_dystonia_v3_25bp.bed \
#	-stand_call_conf 30 \
#	--output_mode EMIT_VARIANTS_ONLY \
#	-maxAltAlleles 10 \
#	-o ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants.vcf

#Decompose variants with Vt
#cd ../vt
#vt decompose ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants.vcf -o ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp.vcf

#Annotate the decomposed variants. Restrict to small bed.
#SnpEff, human genome 37
#Use following command to look at human databases: java -jar snpEff.jar databases | less | grep "Homo*" 
#cd ../snpEff
#java -Xmx4g -jar snpEff.jar \
#-fi ../../example_fastqs/beds/NGD_dystonia_v3_25bp.bed \
#GRCh37.75 \
#../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp.vcf > \
#../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp_ann.vcf
#
#SnpSift for ExAc
#cd ../snpEff
#java -jar SnpSift.jar annotate -v ../../reference_files/exac/ExAC.r0.3.1.sites.vep.vcf.gz \
# ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp.vcf > \
#../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp_exac.vcf

#After: filter variants above 5%
#cd ../snpEff
#cat ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_decomp_exac.vcf | java -jar #SnpSift.jar filter "(#filtered thing)" | > \ 
#../../home/stpuser/vcf/filtered.vcf


#Print end date and time for the console record
echo '' 
echo Pipeline run completed at: $(date)


















