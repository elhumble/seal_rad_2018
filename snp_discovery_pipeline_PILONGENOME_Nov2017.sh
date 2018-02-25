# SNP discovery for RADseq
# Emily Humble 2017

#~~~~~
# SERVER: gaiking
# Two rounds of sequencing from CEBITEC (same libraries)
#~~~~~

# Working in directory with raw RAD sequencing reads (SAMPLES_PILON/)
# There are two files per individual (LibX, LibX.1) from each sequencing run


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    BWA                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make REFERENCE folder in top level Seal_RADseq directory

mkdir ../../../reference
bwa index submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa


#~~~~~

# 1. BWA MEM

find ./ -name ‘*.fq’

nohup find ./ -name '*.fq' | grep -v .2.fq | sed 's/.1.fq.fq//' | parallel bwa mem ../../../reference/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa {}.1.fq.fq {}.2.fq '>' '{}'.sam &


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

# All SNP filtering in /home/emily/SNP_filtering directory in Rproject SNP_filtering.Rproj


# END
