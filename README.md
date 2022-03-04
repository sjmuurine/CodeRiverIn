# CodeRiverIn

Data analysis of 16S rRNA gene sequences and ARG &amp; MGEs abundances (qPCR array) of Code river in Indonesia 

Contents

This repository contains R code and data files from the study by S. Johanna Muurinen, Windi I. Muziasari, Jenni Hultman, 
Katariina Pärnänen, Vanny Narita, Christina Lyra, Lintang N. Fadlillah Ludhang P. Rizki, William Nurmi, James M. Tiedje, Iwan Dwiprahasto, 
Pramono Hadi and Marko Virta about the antibiotic resistance and river health. 

Included files


Rarefied and subsampled OTU table (Phylo.tx.1.pick.1.subsample.shared_correctedID.txt)
TSS normalized OTU table (Phylo.tx.1.pick.shared_correctedID.txt)
Metadata for OTU tables (ZnCu_metadata.txt)
Taxonomy data (Phylo.cons.taxonomy)
SmartChipTM Real-Time PCR (Takara Bio) results (7 files):
  Ind1-2017.txt
  Ind2-2017.txt
  Ind3-2017.txt
  Ind4-2017.txt
  Ind5-2017.txt
  Ind6-2017.txt
  all_Ichips (Mean and sd calculted for each sample in R using text files above)
Metadata for ARGs and MGEs: Ind_metadata.txt
ARG and MGE annotation file: Ind_assays.txt
16S rRNA gene sequencing results (OTU data table): otutable_indonesia.txt
Script for data analysis in R that uses all the files above: IndonesiaF.R
