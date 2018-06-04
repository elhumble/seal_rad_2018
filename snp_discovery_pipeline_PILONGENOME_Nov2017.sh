# SNP discovery for RADseq
# Emily Humble 2017

#~~~~~
# SERVER: gaiking
# Two rounds of sequencing from CEBITEC (same libraries)
#~~~~~

/external/data/Emily/Seal_RADseq/cebitec_2016/AFS/Mar_2016
/external/data/Emily/Seal_RADseq/cebitec_2016/AFS/June_2016

# in each of these are the same six directories for the six libraries:

/AFS-RAD-Pool-1
/AFS-RAD-Pool-2
/AFS-RAD-Pool-3
/AFS-RAD-Pool-4
/AFS-RAD-Pool-5
/AFS-RAD-Pool-6


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    SAMPLE PREPARATION                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# For files in both Mar_2016 & June_2016, run the following:

# (1) FASTX-TOOLKIT
# First run fastqc to check read quality. If ends are bad, we can trim them
# In each lib directory, run the following on the READ1s and the READ2s

mkdir fastqc
~/programs/FastQC/fastqc Read1_filename -out fastqc/
~/programs/FastQC/fastqc Read2_filename -out fastqc/

# (2) STACKS
# process_radtags including trimming to 225 bp
# libraries are split into 6 folders therefore we just need the within library barcodes. Barcodes file saved in /AFS/barcodes_single:
ACACTGAC
ACGTAGCA
CACACAGT
CAGTCTCA
GTACTCGT
GTCATGTG
TGCATCGT
TGTGACTG
ACACGACA
ACGTCTAC
CAGTGTGT
CATGATCA
GTACGCTG
GTCAGTGT
TGTGTGAC
TGTGCAGT

# use default value for barcodes (--inline_null: barcode is inline with sequence, occurs only on single-end read, default)
# each library is run separately; in each lib folder, for both March and June run:

nohup process_radtags -1 AFS-RAD-Pool-1_ATCACG_L001_R1_001.fastq -2 AFS-RAD-Pool-1_ATCACG_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &
nohup process_radtags -1 AFS-RAD-Pool-2_CGATGT_L001_R1_001.fastq -2 AFS-RAD-Pool-2_CGATGT_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &
nohup process_radtags -1 AFS-RAD-Pool-3_GATCAG_L001_R1_001.fastq -2 AFS-RAD-Pool-3_GATCAG_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &
nohup process_radtags -1 AFS-RAD-Pool-4_GGCTAC_L001_R1_001.fastq -2 AFS-RAD-Pool-4_GGCTAC_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &
nohup process_radtags -1 AFS-RAD-Pool-5_TCGTTA_L001_R1_001.fastq -2 AFS-RAD-Pool-5_TCGTTA_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &
nohup process_radtags -1 AFS-RAD-Pool-6_TACCGA_L001_R1_001.fastq -2 AFS-RAD-Pool-6_TACCGA_L001_R2_001.fastq -o ./ -b ../../barcodes_single.txt -c -q -r -e sbfI -E phred33 -t 225 &

# save all the log files into one directory on my local computer, changing the filenames accordingly eg:

/home/Dropbox/phd/projects/rad/rad_seals/log_files/June_2016:
process_radtags_Mar_2016_1.log
process_radtags_Mar_2016_2.log etc

# can then run rad_tag_depth.R script to determine estimated sequencing coverage

# (3) Rename every sample_ACGTAGCA.2.fq file. At the moment sample names are the same across the 6 libraries although correspond to different individuals
# in each library folder, run for loop to change all names to contain the library name and sequencing run e.g. lib1 (March) and lib1.1 (June)
# will merge all rounds of sequencing into the same directory later

# Mar_2016 files:
for f in sample*; do mv "$f" "${f/sample_/sample_lib1_}"; done
for f in sample*; do mv "$f" "${f/sample_/sample_lib2_}"; done

# ETC...

# June_2016 files:
for f in sample*; do mv "$f" "${f/sample_/sample_lib1.1_}"; done
for f in sample*; do mv "$f" "${f/sample_/sample_lib2.1_}"; done

# ETC...

# (4) remove all *rem* files to one folder (in each directory)
mkdir rem_folder
mv folder/*rem* rem_folder

#~~~~~
# SAMPLES folder
#~~~~~

# (5) Put all individual files (i.e. read 1 and read 2 for each individual, such as: sample_lib6_ACACGACA.1.fq and sample_lib6_ACACGACA.2.fq) into one folder

mkdir SAMPLES_PILON

# Mar_2016
cp March_2016/AFS-RAD-Pool-1/sample_lib1* SAMPLES_PILON/
cp March_2016/AFS-RAD-Pool-2/sample_lib2* SAMPLES_PILON/
cp March_2016/AFS-RAD-Pool-3/sample_lib3* SAMPLES_PILON/
cp March_2016/AFS-RAD-Pool-4/sample_lib4* SAMPLES_PILON/
cp March_2016/AFS-RAD-Pool-5/sample_lib5* SAMPLES_PILON/
cp March_2016/AFS-RAD-Pool-6/sample_lib6* SAMPLES_PILON/

# June_2016
cp June_2016/AFS-RAD-Pool-1/sample_lib1.1_* SAMPLES_PILON/
cp June_2016/AFS-RAD-Pool-2/sample_lib2.1_* SAMPLES_PILON/
cp June_2016/AFS-RAD-Pool-3/sample_lib3.1_* SAMPLES_PILON/
cp June_2016/AFS-RAD-Pool-4/sample_lib4.1_* SAMPLES_PILON/
cp June_2016/AFS-RAD-Pool-5/sample_lib5.1_* SAMPLES_PILON/
cp June_2016/AFS-RAD-Pool-6/sample_lib6.1_* SAMPLES_PILON/

# count reads
cd SAMPLES/
wc -l *
 1810000 sample_lib1.1_ACACGACA.1.fq
     1810000 sample_lib1.1_ACACGACA.2.fq
     2236884 sample_lib1.1_ACACTGAC.1.fq
     2236884 sample_lib1.1_ACACTGAC.2.fq


# All individual raw reads are now in the /SAMPLES folder. There are two files per individual (LibX, LibX.1) from each sequencing run
# I will map the different sequencing runs separately in order to assign the read groups correctly before merging the different runs within individuals
# The read groups may not make a big difference to the genotyping but it is probably good practice to keep this information.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    BWA                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Working in directory with raw RAD sequencing reads (SAMPLES_PILON/)
# There are two files per individual (LibX, LibX.1) from each sequencing run


# make REFERENCE folder in top level Seal_RADseq directory

mkdir ../../../reference
bwa index submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa


#~~~~~

# 1. BWA MEM

find ./ -name ‘*.fq’

nohup find ./ -name '*.fq' | grep -v .2.fq | sed 's/.1.fq//' | parallel bwa mem ../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa {}.1.fq {}.2.fq '>' '{}'.sam &


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               PREPARE SAM FILES FOR GATK                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# 1. Remove unmapped reads:

nohup sh -c 'for i in sample*_A*.sam; do samtools view -bF 4 $i > ${i%.sam}_filtered.sam; done' &
nohup sh -c 'for i in sample*_T*.sam; do samtools view -bF 4 $i > ${i%.sam}_filtered.sam; done' &
nohup sh -c 'for i in sample*_G*.sam; do samtools view -bF 4 $i > ${i%.sam}_filtered.sam; done' &
nohup sh -c 'for i in sample*_C*.sam; do samtools view -bF 4 $i > ${i%.sam}_filtered.sam; done' &



# PICARD TOOLS

#2. Run SortSam.jar: sorts sam files by coordinate and makes a bam file
# where is SortSam.jar? /home/emily/programs/ it is executable

nohup sh -c 'for i in *lib1_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib1.1_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib2_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib2.1*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib3_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib3.1*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib4_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib4.1*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib5_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib5.1*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib6_*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &

nohup sh -c 'for i in *lib6.1*filtered.sam; do java -Xmx8g -jar /usr/local/bin/picard.jar SortSam I=$i O=${i%_filtered.sam}_filtered_sorted.bam SORT_ORDER=coordinate; done' &


for i in *.bam; do samtools flagstat $i >> filtered_sorted_bam_flagstat.txt; done


#~~~~~~~

# add read groups to file so GATK can process them
# RGSM is important: individual identifier
# where is AddOrReplaceReadGroups.jar? /home/emily/programs/ it is executable

# Read group information:
# ID: read group identifier (should be unique across read groups). picard tools will modify the ID when merging to avoid collisions (libX_BARCODE)
# PL: platform. (illumina)
# LB: library id, (libX or libX.1)
# PU: platform unit. This contains as much info as possible (PU:lib1.1_TGCATCGT)
# SM: sample. Sequencing data from the same sample, this is the name used in the vcf file (libX_BARCODE)


nohup sh -c 'for i in sample*_A*_sorted.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups INPUT=$i OUTPUT=${i%bam}RG.bam RGID=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGPL=illumina RGLB=`ls -1 "$i" | cut -f 2 -d '_'` RGPU=`ls -1 "$i" | cut -f 2 -d '_'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGSM=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` VALIDATION_STRINGENCY=SILENT; done' &


nohup sh -c 'for i in sample*_G*_sorted.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups INPUT=$i OUTPUT=${i%bam}RG.bam RGID=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGPL=illumina RGLB=`ls -1 "$i" | cut -f 2 -d '_'` RGPU=`ls -1 "$i" | cut -f 2 -d '_'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGSM=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` VALIDATION_STRINGENCY=SILENT; done' &


nohup sh -c 'for i in sample*_C*_sorted.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups INPUT=$i OUTPUT=${i%bam}RG.bam RGID=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGPL=illumina RGLB=`ls -1 "$i" | cut -f 2 -d '_'` RGPU=`ls -1 "$i" | cut -f 2 -d '_'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGSM=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` VALIDATION_STRINGENCY=SILENT; done' &


nohup sh -c 'for i in sample*_T*_sorted.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups INPUT=$i OUTPUT=${i%bam}RG.bam RGID=`ls -1 "$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGPL=illumina RGLB=`ls -1 "$i" | cut -f 2 -d '_'` RGPU=`ls -1 "$i" | cut -f 2 -d '_'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` RGSM=`ls -1
"$i" | cut -f 2 -d '_' | cut -f 1 -d '.'`"_"`ls -1 "$i" | cut -f 3 -d '_' | cut -f 1 -d '.'` VALIDATION_STRINGENCY=SILENT; done' &



# Could consider changing the RGID to include the specific library run. Although it shouldn’t matter as Picard MergeSamFiles deals with this by modifying the name.
# Check read group information:

samtools view -H sample_lib1.1_TGCATCGT_filtered_sorted.RG.bam | grep '@RG'


#~~~~~~~

# 3. Remove PCR Duplicates


nohup sh -c 'for i in sample_*_G*.RG.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT=$i OUTPUT=${i%bam}rmdup.bam METRICS_FILE=${i%bam}rmdupmetrics READ_NAME_REGEX=null MAX_RECORDS_IN_RAM=500000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 > ${i%bam}rmduplog 2>&1; done' &
[1] 7134

nohup sh -c 'for i in sample_*_A*.RG.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT=$i OUTPUT=${i%bam}rmdup.bam METRICS_FILE=${i%bam}rmdupmetrics READ_NAME_REGEX=null MAX_RECORDS_IN_RAM=500000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 > ${i%bam}rmduplog 2>&1; done' &
[2] 7337

nohup sh -c 'for i in sample_*_T*.RG.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT=$i OUTPUT=${i%bam}rmdup.bam METRICS_FILE=${i%bam}rmdupmetrics READ_NAME_REGEX=null MAX_RECORDS_IN_RAM=500000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 > ${i%bam}rmduplog 2>&1; done' &
[3] 7545

nohup sh -c 'for i in sample_*_C*.RG.bam; do java -Xmx8g -jar /usr/local/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT INPUT=$i OUTPUT=${i%bam}rmdup.bam METRICS_FILE=${i%bam}rmdupmetrics READ_NAME_REGEX=null MAX_RECORDS_IN_RAM=500000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 > ${i%bam}rmduplog 2>&1; done' &
[4] 7694


#~~~~~~~

# 4. Merge individual bam files

# first move all the current files into a directory

mkdir sorted_RG_rmdup
mv *filtered_sorted.RG.rmdup.bam sorted_RG_rmdup/
cd sorted_RG_rmdup

# get the file names for libX

for i in {1..6}; do ls *lib${i}_*.bam > lib${i} ; done

# get the file names for libX.1

for i in {1..6}; do ls *lib${i}.1*.bam > lib${i}.1 ; done

# merge same individuals in each lib (the quotes must be correct here otherwise the regex don’t work).
# merged out should be e.g. sample_lib5_TGTGTGAC_merged.bam

nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib1.1 3<lib1' &


nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib2.1 3<lib2' &


nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib3.1 3<lib3' &


nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib4.1 3<lib4' &

nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib5.1 3<lib5' &


nohup sh -c 'while IFS= read -r lineA && IFS= read -r lineB <&3; do
  java -Xmx8g -jar /usr/local/bin/picard.jar MergeSamFiles  I="$lineA" I="$lineB" O=`ls -1 ${lineB}|cut -f 1,2 -d '_'`"_"`ls -1 ${lineB}|cut -f 3 -d '_'`_merged.bam
done <lib6.1 3<lib6' &


# can check the output here e.g.

samtools view -H sample_lib5_TGTGTGAC_merged.bam | grep ‘@RG’

# RGSM should be identical

#~~~~~~~

# Move all files into one place for SNP calling

mkdir ../GATK
mv *_merged.bam ../GATK/

rm lib*

# should be 96 files in /GATK
ls | wc -l



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          GATK                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Prepare reference folder
# 1) Need 3 reference files:

# a) fasta file(s) used in BWA
# b) indexed fasta: reference.fasta.fai: run this samtools code:

samtools faidx ../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa

# c) reference.dict: run this picardtools code:

java -Xmx8g -jar /usr/local/bin/picard.jar CreateSequenceDictionary R= submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa O= submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.dict


#~~~~~~~
# in GATK folder
#~~~~~~~

# index bam files
nohup sh -c 'for i in *.bam; do samtools index $i; done' &

# file structure in GATK folder:

# *merged.bam  *merged.bam.bai
# Realigner.sh  (script)

# Now ready to start your GATK runs. First step: Realigner.sh



# Realigner.sh script:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#!/bin/bash

java -Djava.io.tmpdir=./Temp/ -Xmx64g -jar /usr/local/bin/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R ../../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa \
  -o merged_output.intervals \
  $(printf ' -I %s' *.bam)

java -Djava.io.tmpdir=./Temp/ -Xmx64g -jar /usr/local/bin/GenomeAnalysisTK.jar \
  $(printf ' -I %s' *.bam) \
  -R ../../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa \
  -T IndelRealigner \
  -targetIntervals merged_output.intervals \
  -nWayOut _realigned.bam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## run with
nohup sh Realigner.sh &




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             Haplotype Caller                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Haplotype caller run for each sample ID

nohup sh -c 'for i in *realigned.bam; do java -Djava.io.tmpdir=./Temp/ -Xmx64g -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller -nct 4 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -hets 0.01 \
-R ../../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa -I $i -o ${i%merged_realigned.bam}merged_realigned.g.vcf.gz; done' &





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             GenotypeGVCFs                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Then genotype all gvcfs together using GenotypeGVCFs

nohup java -Djava.io.tmpdir=./Temp/ -Xmx64g -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T GenotypeGVCFs -nt 8 \
-R ../../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa -o ArcGaz_pilon_genotype_gvcf.vcf.gz $(printf ' -V %s' *.g.vcf.gz) &






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                      SNP filtering                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# in GATK folder, unzip final vcf file

gzip -d ArcGaz_genotype_gvcf.vcf.gz

# All SNP filtering in /home/emily/rad_analysis_pilon_Dec2017 directory in Rproject rad_analysis_pilon_Dec2017.Rproj


# END
