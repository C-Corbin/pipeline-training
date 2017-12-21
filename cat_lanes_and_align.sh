#################
#    FASTQC     #
#################

#Start in software folder
#Run Fastqc and write to own folder
#cd FastQC
#./fastqc ../../example_fastqs/*fastq.gz --outdir=/home/stpuser/fastqcs


#########################
# ALIGN AND MERGE LANES #
#########################
#Merging lanes after alignment is more 'best practice'

#Run bwa mem on R1 R2 L001 L002 for the samples, write to ../../home/stpuser/aligned_seqs
#The genome reference is in ../../reference_files/ucsc.hg19.nohap.masked.fasta
#cd ../software/bwa-0.7.15
#
#
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


######################################
# DEFUNCT; CONCATENATING NOT MERGING #
######################################

#DEFUNCT: Concatenate lanes, then run bwa mem - this gets rid of read group information
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R2_both.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R2_both.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R2_both.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R2_both.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R1_both.fastq.gz | samtools view -bh -o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam


#DEFUNCT: Join lanes 1 and 2, writes to the joined_lanes folder 
#cat 1504850-S1509352-02_GCTCGGTA_L001_R1_001.fastq.gz 1504850-S1509352-02_GCTCGGTA_L002_R1_001.fastq.gz > ../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R1_both.fastq.gz
#cat 1504850-S1509352-02_GCTCGGTA_L001_R2_001.fastq.gz 1504850-S1509352-02_GCTCGGTA_L002_R2_001.fastq.gz > ../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R2_both.fastq.gz
#cat 1606034-S1612259-02_CGAACTTA_L001_R1_001.fastq.gz 1606034-S1612259-02_CGAACTTA_L002_R1_001.fastq.gz > ../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R1_both.fastq.gz
#cat 1606034-S1612259-02_CGAACTTA_L001_R2_001.fastq.gz 1606034-S1612259-02_CGAACTTA_L002_R2_001.fastq.gz > ../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R2_both.fastq.gz
#cat 1607686-S1615531-02_TTCACGCA_L001_R1_001.fastq.gz 1607686-S1615531-02_TTCACGCA_L002_R1_001.fastq.gz > ../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R1_both.fastq.gz
#cat 1607686-S1615531-02_TTCACGCA_L001_R2_001.fastq.gz 1607686-S1615531-02_TTCACGCA_L002_R2_001.fastq.gz > ../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R2_both.fastq.gz
#cat 1609778-S1620040-02_CGCTGATC_L001_R1_001.fastq.gz 1609778-S1620040-02_CGCTGATC_L002_R1_001.fastq.gz > ../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R1_both.fastq.gz
#cat 1609778-S1620040-02_CGCTGATC_L001_R2_001.fastq.gz 1609778-S1620040-02_CGCTGATC_L002_R2_001.fastq.gz > ../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R2_both.fastq.gz
#cat 1703057-S1705957-02_AAGACGGA_L001_R1_001.fastq.gz 1703057-S1705957-02_AAGACGGA_L002_R1_001.fastq.gz > ../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R1_both.fastq.gz
#cat 1703057-S1705957-02_AAGACGGA_L001_R2_001.fastq.gz 1703057-S1705957-02_AAGACGGA_L002_R2_001.fastq.gz > ../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R2_both.fastq.gz


#########################
#  SORT AND INDEX SEQS  #
#########################

#Sort the files
#cd ../samtools-1.3.1
#samtools sort -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam

#Index the files for IGV viewing
#samtools index ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam


###################
# MARK DUPLICATES #
###################

#I need to mark duplicates for each file
#This also writes a metrics file for each flagged file
#cd ../picard-tools-2.5.0 
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam O=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam O=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam O=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam O=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam O=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_metrics.txt


############
# FLAGSTAT #
############

#samtools flagstat for quality metrics - to see how they match up with the picard ones!
#cd ../samtools-1.3.1
#samtools flagstat ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_flagstat.txt
#samtools flagstat ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam > ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_flagstat.txt


###################
# ADD READ GROUPS #
###################
#In a real pipeline, read group information can be used to track how particular machines/runs/e.t.c are performing; they may also allow correct calibration?
#However in our 4-sample example we can probably put in placeholder read groups
#Adding them quite late - 'real' ones would be added early, while still at the lanes-and-reads stage

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
	

#####################
# INDEL REALIGNMENT #
#####################

###INDEXING - needed for indel-realignment files
#cd ../samtools-1.3.1
#samtools index ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_RG.bam
#samtools index ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_RG.bam


### IDENTIFY TARGET REGIONS
### CAN RESTRICT TO THE BROAD BED FILE WITH -L, the GATK 'intervals' option
#
#1504850-S1509352-02_GCTCGGTA - NGD Complete
#cd ../GenomeAnalysisTK-3.6
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


##  INDEL REALIGNMENT Test run
#cd ../GenomeAnalysisTK-3.6
#java -jar GenomeAnalysisTK.jar \
#	-T IndelRealigner \
#	-R ../../reference_files/ucsc.hg19.nohap.masked.fasta \
#	-L ../../example_fastqs/beds/NGD_CompletePanel_25bp_v2.bed \
#	-I ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.bam \
#	-targetIntervals ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG.intervals \
#	-o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned.bam


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

##################################
# VARIANT CALLING AND ANNOTATING #
##################################
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

#Annotate with SnpEff (used in our pipeline)
cd ../snpEff
#Look at human databases - GRCh37.75 available
#java -jar snpEff.jar databases | less | grep "Homo*" 
#Annotate to human genome 37 with dbSNP annotations. Specify intervals broad bed. DON'T ADD dbSNP NUMBERS YET!
java -Xmx4g -jar snpEff.jar \
-fi ../../example_fastqs/beds/NGD_dystonia_v3_25bp.bed \
GRCh37.75 \
../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants.vcf > ../../home/stpuser/vcf/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_RG_realigned_variants_ann.vcf

#Next: decompose the ALT column to put second alleles on a spare row, and then annotate with dbSNP
#Filter out those with dbSNP frequency above 5%






























