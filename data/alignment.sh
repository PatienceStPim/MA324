#!/bin/bash
mkdir alignment alignment/clustalo alignment/mafft alignment/einsi
#Create alignment using clustal omega on default settings.
clustalo -i Norovirus_PP.fasta -o alignment/clustalo/Norovirus_clustalo_.fst
#Alignment by MAFFT in sequential mode
mafft --maxiterate 1000 Norovirus_PP.fasta > alignment/mafft/Norovirus_mafft.fst
#Alignment by MAFFT in E-INS-i mode
mafft --maxiterate 1000 --genafpair Norovirus_PP.fasta > alignment/einsi/Norovirus_einsi.fst
#Survey each MSA using AliStat in amino acid mode (1)
#Clustalo
alistat alignment/clustalo/Norovirus_clustalo_.fst 6 -i -o alignment/clustalo/clustalo
#MAFFT in sequential mode 
alistat alignment/mafft/Norovirus_mafft.fst 6 -i -o alignment/mafft/mafft
#MAFFT in E-INS-i mode
alistat alignment/einsi/Norovirus_einsi.fst 6 -i -o alignment/einsi/einsi