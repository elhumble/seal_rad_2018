#!/bin/bash

# make extra file

awk '{print $2'} $1 > pos.txt
sed 's/:.*//' pos.txt > pos2.txt
sed -i -e 's/Contig//g' pos2.txt

# remove first column

awk '{print substr($0, index($0, $2))}' $1 > ArcGaz_genotype_biallelic.2.map

# remove chr from first column

paste pos2.txt ArcGaz_genotype_biallelic.2.map | column -s $'\t' -t > $2

# clean up
rm ArcGaz_genotype_biallelic.2.map pos.txt pos2.txt

