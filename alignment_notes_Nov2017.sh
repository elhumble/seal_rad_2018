# Visualisation and characterisation of PBJelly genome v.1.4 to CanFam3.1
# Emily Humble
# Nov 2017


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alignment to the DOG genome
# in  /PBJelly_v1.4_CanFam3.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# (1) Download CanFam3.1 genome
# CanFam downloaded onto server into:

# /assemblies/CanFam3.1/

# (2) Run LAST alignment
# See the following for support:
# http://last.cbrc.jp/tutorial/last-tutorial-2.html
# http://last.cbrc.jp/doc/
# https://wiki.bils.se/wiki/Running_the_LAST_aligner_on_UPPMAX

# In directory /alignments/LAST/PBJelly_v1.4_CanFam3.1
# First make CanFam genome database
# -c is repeat masking scheme
# -u is for seeding scheme -- uMAM4 doesn't use too much memory

nohup nice lastdb -cR11 -uMAM4 CanFam3.1.lastdb /home/emily/assemblies/CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.fna &

# Lastal command

mkdir temp
nohup parallel-fasta --tmpdir=./temp 'lastal -e120 ../PBJelly_CanFam3.1/CanFam3.1.lastdb | last-split' < /home/emily/assemblies/ArcGaz_PBJelly_v1.4/submission.assembly.ArcGaz001_AP3.GapClosed.ArcGaz.paired.PBJelly2.utg000001c.utg.quiver.nopipe.pilon.fa > PBJelly_v1.4_CanFam3.1.maf &

# (3) Prepare files for CIRCOS plotting

# PBJelly_CanFam3.1

# Extract the alignment lines
grep '^s' PBJelly_v1.4_CanFam3.1.maf > temp1
# Get rid of sequence
awk 'NF{NF-=1};1' <temp1 >temp2
# Get rid of starting s
cut -d " " -f2-6 temp2 > temp3

# Get PBJelly coords
grep '^Contig' temp3 > PBJelly_v1.4_coords.txt
# Get CanFam coords
grep '^N' temp3 > CanFam_coords.txt

rm temp*
  
# check line numbers
wc -l PBJelly_v1.4_coords.txt
wc -l CanFam_coords.txt

# Sequence similarity with MAFfilter

# First change seq headers in maf file

sed -i 's/^s\sNC/s CanFam.NC/g’ PBJelly_v1.4_CanFam3.1.maf
sed -i 's/^s\sContig/s ArcGaz.Contig/g’ PBJelly_v1.4_CanFam3.1.maf

# run maffilter

maffilter input.file=PBJelly_v1.4_CanFam3.1.maf input.file.compression=none output.log=PBJelly_v1.4_CanFam3.1.maffilter.log maf.filter="SequenceStatistics(statistics=(PairwiseDivergence(species1=ArcGaz, species2=CanFam)), ref_species=ArcGaz, file=PBJelly_v1.4_CanFam3.1.statistics.csv)"

