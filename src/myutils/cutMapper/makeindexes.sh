#!/bin/bash

GENOME="/data/downloaded/HumanGenome/hg18Ref.2bit"
SNPFILE="1kSNP.YRI.myvcf"

IndexFolder="hg18indexSNPoneComb"

for i in {0..255}; 
do 
    echo "./cutMapperMakeIndex -one -idx=$i $GENOME $SNPFILE  $IndexFolder "
    echo "./cutMapperMakeIndex -one -idx=$i $GENOME $SNPFILE  $IndexFolder " | 
         qsub -l h_vmem=3g,bigio=1 -wd `pwd` -N "indexing"
done
