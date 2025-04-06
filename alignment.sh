#!/bin/bash
#Make directories to hold outputs for different MSA methods
mkdir alignment alignment/clustalo alignment/mafft alignment/einsi
#Create alignment using clustal omega on default settings, output to 'clustalo' directory
clustalo -i Norovirus_PP.fasta -o alignment/clustalo/Norovirus_clustalo.fst
#Alignment by MAFFT in sequential mode, output to 'mafft' directory
mafft --maxiterate 1000 Norovirus_PP.fasta > alignment/mafft/Norovirus_mafft.fst
#Alignment by MAFFT in E-INS-i mode, output to 'einsi' directory
mafft --maxiterate 1000 --genafpair Norovirus_PP.fasta > alignment/einsi/Norovirus_einsi.fst

#Survey each MSA using AliStat in amino acid mode (6) and store the output in the relevant directory (alistat). We are interested in Cc histograms and summary (-t 2)
#Clustalo
for i in clustalo mafft einsi
do
    mkdir alignment/$i/alistat 
    alistat alignment/$i/Norovirus_$i.fst 6 -t 2 -o $i
    R CMD BATCH $i.Histogram_Cc.R
    mv $i* alignment/$i/alistat
done

#Masking using Alistat by outputting the alignment where sites have a Cc > 0.3 based off Cc distribution. Masked fasta file and summary.txt stored in relevant directory (alistat_masked)
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

#Run unmasked alignment through SatuRation and SatuRationHeatMapper in amino acid mode (30), moves output to relevant directory (saturation)
mkdir alignment/mafft/saturation
saturation alignment/mafft/Norovirus_mafft.fst v f 30
SatuRationHeatMapper.pl -i alignment/mafft/Norovirus_mafft_lambda.csv -f
mv alignment/mafft/Norovirus_mafft_* alignment/mafft/saturation

#Run alignment through Homo to obtain p-values for homogenous evolutionary conditions in amino acid mode (30) and error rate 0.05, moves output to relevant directory (homo)
mkdir alignment/mafft/homo
homo alignment/mafft/alistat_masked/mafft.Mask.fst f 30 0.05
mv alignment/mafft/alistat_masked/mafft_* alignment/mafft/homo

#Make directory to store our tree data and move required alignment file to it
mkdir trees 
cp alignment/mafft/alistat_masked/mafft.Mask.fst trees
#Model selection using model finder in iqtree2 (-m MFP), incorporating Lie Markov Models (+LM) tailored for viral amino acid sequences (--seqtype AA) and inference of phylogenetic tree with 1000 bootstrap replicates (-B 1000)
iqtree2 -s trees/mafft.Mask.fst --seqtype AA -v -m MFP+LM --msub viral -B 1000 -nt AUTO
#Note that the sequences for Marmot_norovirus were removed from the alignment file manually and passed through iqtree.
