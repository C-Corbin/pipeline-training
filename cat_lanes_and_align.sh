#run from example_fastqs directory
#join lanes 1 and 2, writes to the joined_lanes folder 
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

#run bwa mem on R1 and R2 for every sample, writes to ../../home/stpuser/aligned_seqs
#the genome reference is in ../../reference_files/ucsc.hg19.nohap.masked.fasta

#cd ../software/bwa-0.7.15
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1504850-S1509352-02_GCTCGGTA_R2_both.fastq.gz | samtools view -b -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1606034-S1612259-02_CGAACTTA_R2_both.fastq.gz | samtools view -b -o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1607686-S1615531-02_TTCACGCA_R2_both.fastq.gz | samtools view -b -o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1609778-S1620040-02_CGCTGATC_R2_both.fastq.gz | samtools view -b -o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam
#bwa mem ../../reference_files/ucsc.hg19.nohap.masked.fasta ../../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R1_both.fastq.gz ../../home/stpuser/joined_lanes/1703057-S1705957-02_AAGACGGA_R1_both.fastq.gz | samtools view -b -o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam

#sort the files
#cd ../samtools-1.3.1
#samtools sort -o ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC.bam
#samtools sort -o ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA.bam

#index the files for IGV viewing
#samtools index ../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam
#samtools index ../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam

#need to mark duplicates for each file
#but first I need a metrics file
cd ../picard-tools-2.5.0 
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted.bam O=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1504850-S1509352-02_GCTCGGTA_sorted_dupflag_metrics.txt
#java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted.bam O=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1606034-S1612259-02_CGAACTTA_sorted_dupflag_metrics.txt
java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted.bam O=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1607686-S1615531-02_TTCACGCA_sorted_dupflag_metrics.txt
java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted.bam O=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1609778-S1620040-02_CGCTGATC_sorted_dupflag_metrics.txt
java -jar picard.jar MarkDuplicates I=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted.bam O=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag.bam M=../../home/stpuser/aligned_seqs/1703057-S1705957-02_AAGACGGA_sorted_dupflag_metrics.txt





