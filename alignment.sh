#!/bin/bash
#Make directories to hold outputs for different MSA methods
mkdir alignment alignment/clustalo alignment/mafft alignment/einsi
#Create alignment using clustal omega on default settings.
clustalo -i Norovirus_PP.fasta -o alignment/clustalo/Norovirus_clustalo.fst
#Alignment by MAFFT in sequential mode
mafft --maxiterate 1000 Norovirus_PP.fasta > alignment/mafft/Norovirus_mafft.fst
#Alignment by MAFFT in E-INS-i mode
mafft --maxiterate 1000 --genafpair Norovirus_PP.fasta > alignment/einsi/Norovirus_einsi.fst

#Survey each MSA using AliStat in amino acid mode (6) and store the output in the relevant directory. We are interested in Cc histograms and summary (-t 2)
#Clustalo
for i in clustalo mafft einsi
do
    mkdir alignment/$i/alistat 
    alistat alignment/$i/Norovirus_$i.fst 6 -t 2 -o $i
    R CMD BATCH $i.Histogram_Cc.R
    mv $i* alignment/$i/alistat
done

#Masking using Alistat by outputting the alignment where sites have a Cc > 0.3 based off Cc distribution. Masked fasta file and summary.txt stored in relevant directory
for i in clustalo mafft einsi
do 
    mkdir alignment/$i/alistat_masked
    alistat alignment/$i/Norovirus_$i.fst 6 -m 0.3 -o $i
    #generated summary file is of the original alignment
    rm $i.Summary.txt
    #rerun alistat on the masked alignment in brief mode
    alistat $i.Mask.fst 6 -o $i
    mv $i* alignment/$i/alistat_masked
done

#Run alignments through SatuRation and SatuRationHeatMapper, moves output to relevant directory
mkdir alignment/mafft/saturation
saturation alignment/mafft/Norovirus_mafft.fst v f 30
SatuRationHeatMapper.pl -i Norovirus_mafft_lambda.csv -f
mv Norovirus_mafft_* alignment/mafft/saturation