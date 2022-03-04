
# INDONESIA

#-------------- Housekeeping ------------------# 

#R-environment is cleared

rm(list=ls())

#Working directory

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Indonesia")
#load("~/Library/Mobile Documents/com~apple~CloudDocs/Indonesia/I_Rdata.RData")


#Libraries and packages
#Before libraries can be loaded, the packages needs to be installed. 
#They can be installed like this:
#install.packages("reshape2")

library(reshape2)
library(dplyr)
library(tidyr)
library(MASS)
library(vegan)
library(gplots)
library(ggplot2)
library(ggthemes)
library(colorspace)
library(colorRamps)
library(RColorBrewer)
library(gridExtra)
library(fitdistrplus)
library(logspline)
library(car)
library(lme4) 
library(multcomp)
library(scales)
#library(ggridges)
#library(viridis)
#library(corrplot)
library(psych)

#-------------- Data & processing ------------------# 

#We define a function to read in each plate and calculate 
#the mean and SD from the three technical replicates. 
#If only 1/3 gives signal (Ct < 27), it is considered as false positive. 
#The function will be called Array_func, but can also be called 
#something else. Just be careful that you don???t overwrite some existing functions.

Array_func <- function(PATH) {
  #read in the chip
  chip <- read.table(PATH, header=TRUE, sep="\t")
  #change all Ct > 27 to NA
  chip$Ct[chip$Ct>27] <- NA
  #make a new data drame for the results
  outArray <- data.frame()
  #loop over the whole chip and calculate the mean and SD
  
  #There are 4 samples/chip and 384 different assays. 
  #The samples have 3 replcates so that rows 1-3 have replicates of Sample1
  #and the Ct values of assay called "16S new 2_2"
  #If your samples and replicates ar in different order, you need to either 
  #change the function accordingly or reorder your .txt file.
  #reordering can be done for example in excel (just save the file as .txt file)
  
  #so first we create an index over the .txt file in every 3 rows
  for (i in seq(1,nrow(chip), 3)) {
    #if 2/3 of Ct-values are NAs, write everything as NA
    #the Ct values are in the column 6 in the .txt -file
    if(is.na(chip[i,6]) && is.na(chip[i+1,6]) || is.na(chip[i+1,6]) && 
       is.na(chip[i+2,6]) || is.na(chip[i,6]) && is.na(chip[i+2,6])) {
      #sample names are in column 4 in the .txt -file and assay names in column 3
      outArray[(i+2)/3, 1] <- chip[i,4]
      outArray[(i+2)/3, 2] <- chip[i,3]
      outArray[(i+2)/3, 3] <- NA
      outArray[(i+2)/3, 4] <- NA
    }
    #If at least 2/3 of Ct values are <27, then the mean and sd are caluculated
    #Note that the loop uses the sample name of replicate 1, so in our data set 
    #the samples are named as Sample#-Rep1. We will deal with that later. 
    else {
      outArray[(i+2)/3, 1] <- chip[i,4]
      outArray[(i+2)/3, 2] <- chip[i,3]
      outArray[(i+2)/3, 3] <- mean(c(chip[i,6], chip[i+1,6], chip[i+2,6]), na.rm=TRUE)
      outArray[(i+2)/3, 4] <- sd(c(chip[i,6], chip[i+1,6], chip[i+2,6]), na.rm=TRUE)
    }
  }
  colnames(outArray) <- c("Sample", "Assay", "mean", "sd")
  return(outArray)
}

#Then we can use the new function to read in 
#all the text files. Do this for all the chips (.txt files) you have.
# Example: chip1 <- Arrray_func("PATH/TO/YOUR/RESULTS/chip1.txt")
#if your chip data is in your working directory like in our case, 
#you dont need to specify the path. 

Ichip1<- Array_func("Ind1-2017.txt")
Ichip2<- Array_func("Ind2-2017.txt")
Ichip3<- Array_func("Ind3-2017.txt")
Ichip4<- Array_func("Ind4-2017.txt")
Ichip5<- Array_func("Ind5-2017.txt")
Ichip6<- Array_func("Ind6-2017.txt")

#Then combine all chips to one data frame

all_Ichips <- rbind(Ichip1, Ichip2, Ichip3, Ichip4, Ichip5, Ichip6)

write.table(all_Ichips, "/Users/jmuurine/Library/Mobile Documents/com~apple~CloudDocs/Indonesia/all_Ichips",sep="\t")

# from now on I'll read in a file called "all_Ichips.txt" like this:

#
all_Ichips<- read.table("all_Ichips", header=TRUE,stringsAsFactors = TRUE, sep="\t")


#Do we have all?
View(all_Ichips)
AYs<-unique(all_Ichips$Assay)
length(AYs)
#[1] 384, yes we do.

Ind_assays <- read.table("Ind_assays.txt", header=TRUE, stringsAsFactors = TRUE, sep="\t")

Ind_results <- merge(all_Ichips, Ind_assays, by="Assay")

str(Ind_results)
#'data.frame':	9216 obs. of  6 variables:
#$ Assay                     : Factor w/ 384 levels "AY1","AY10","AY100",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ Sample                    : Factor w/ 72 levels "1A-rep1","1A-rep2",..: 1 4 7 55 22 25 28 31 16 19 ...
#$ mean                      : num  24.3 11.7 13.6 12 23.3 ...
#$ sd                        : num  0.596 0.14 0.576 0.513 0.662 ...
#$ Gene                      : Factor w/ 384 levels "16S rRNA_1","16S rRNA_2",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ Target.antibiotics..major.: Factor w/ 23 levels "(Flouro)Quinolone",..: 2 2 2 2 2 2 2 2 2 2 ...

#Do we have everything?
9216 / 384
#[1] 24
 24/8
#[1] 3
#yes!
 
 View(Ind_results)
 levels(Ind_results$Gene)
 
 #the tnpA genes can be linked to IS-elements, somehow this info has not been attached to their sequences. So lets change the names:
 #Currently they are:
 #"tnpA_1" "tnpA_2"  "tnpA_3" "tnpA_4"  "tnpA_5" "tnpA_6"  "tnpA_7"  
 
 #I'll change them according to previous publications:
 
 levels(Ind_results$Gene)<-c(levels(Ind_results$Gene),"tnpA-05/IS26",
                              "tnpA-04/IS6100", "tnpA-02/IS4","tnpA-01/IS21", "tnpA-07/ISEcp1B",
                              "tnpA-06/IS1216", "tnpA-03/IS6")
 
 Ind_results$Gene[Ind_results$Gene=="tnpA_5"]<-"tnpA-05/IS26"  
 Ind_results$Gene[Ind_results$Gene=="tnpA_4"]<-"tnpA-04/IS6100"  
 Ind_results$Gene[Ind_results$Gene=="tnpA_2"]<-"tnpA-02/IS4" 
 Ind_results$Gene[Ind_results$Gene=="tnpA_1"]<-"tnpA-01/IS21" 
 Ind_results$Gene[Ind_results$Gene=="tnpA_7"]<-"tnpA-07/ISEcp1B"
 Ind_results$Gene[Ind_results$Gene=="tnpA_6"]<-"tnpA-06/IS1216"
 Ind_results$Gene[Ind_results$Gene=="tnpA_3"]<-"tnpA-03/IS6" 
 
 #lets see what we have here:
 levels(Ind_results$Target.antibiotics..major.)
 
 
 #[1] "(Flouro)Quinolone"      "16S rRNA"               "Aminoglycoside"         "Amphenicol"             "Beta Lactam"            "Integron"              
 #[7] "MDR"                    "MGE"                    "MLSB"                   "Other (Antiseptic)"     "Other (Bacitracin)"     "Other (Fosfomycin)"    
 #[13] "Other (Housekeeping)"   "Other (Mercury)"        "Other (Nissin)"         "Other (Nitroimidazole)" "Other (Streptothricin)" "Other (Triclosan)"     
 #[19] "Quinolone"              "Sulfonamide"            "Tetracycline "          "Trimethoprim"           "Vancomycin"
 
 #This could be modified too but I'll do it later. 
 
 #Now I have all Ct values in one dataframe with the annotation of 
 #the different assays. Next I will make a matrix of the results 
 #that can be used to make e.g. ordination plots in vegan. 
 
 #I will combine all mean values to a matrix. 
 #Rows are samples and colums assays. 
 #If some assays have the same gene name (column Gene)
 #they will be aggregated. If this is the case, we combine the name and assay number:
 #result_matrix <- dcast(final_plates, Sample~Gene+Assay, value.var="mean")
 
 Ind_result_mat<- dcast(Ind_results, Sample~Gene, value.var="mean")
 dim(Ind_result_mat)
 #[1]  24 385
 #1st column is the sample, and we have 384 assays.
 
 #Next task is to calculate relative abundances:
 
 plot(Ind_result_mat[,2:3], main="1st 16S vs. 2nd 16S")
 plot(Ind_result_mat[,2], main="1st 16S (old)")
 plot(Ind_result_mat[,3], main="2nd 16S (new)")
 
 #Seems that the "old" 16S is a good choice for normalization,
 #it is detected in 23 samples and NA in only one "pristine sample". 
 #The Ct-values were 27.84495 and 27.72568 in two technical replicates
 #(no amplification in one)
 #I will make an exception and replace the value manually so that I can get 
 #something to sample 1B so I can calculate relative abundances 
 #(although it is highly unlikely that there are ARGs)
 (27.84495+27.72568)/2
# [1] 27.78532
 Ind_result_mat[2,2]
 Ind_result_mat[2,2]<-27.78532
 
 #Now let's calculate the relative abundance. (R = 2^-deltaCt)
 #We could also calculate the relative abundances based on the new 16S.
 #Just change the column number.

 #I have 384 assays, assays other than 16S start from column 4, the last colunm is 385. 
 #New 16S in in column 3 and the old in column 2 (We will use the old)
 Ind_results_delta <- Ind_result_mat[,4:385]- (Ind_result_mat[,2])
 Ind_results_relative <- 2^-(Ind_results_delta)
 
 #Now we have all the relative abundances calculated and we can remove 
 #all the assays that are not detected in any sample can be deleted.
 
 Ind_results_relative <- Ind_results_relative[,colSums(is.na(Ind_results_relative)) != nrow(Ind_results_relative)]
 
 #Then check how many assays were left with dim() and after that give all NAs an artificial value 
 #(NA is not the same as below detection limit, so we need a number). 
 #For example deltaCt of 20, which gives relative value of 9.536743e-07 (2^-20).
 #This can be anything, but probably should be smaller than the smallest real result.
 #That can be checked from the result_deltaCt. So the maximum deltaCt from all results.
 dim(Ind_results_relative)
# 24 and 187
 max(Ind_results_delta, na.rm=TRUE)
 #[1] 15.00762
 Ind_results_relative[is.na(Ind_results_relative)] <- 9.536743e-07
 Ind_results_relative_mat<-as.matrix(Ind_results_relative)
 
 min( Ind_results_relative_mat[ Ind_results_relative_mat!=min( Ind_results_relative_mat)])
 #[1] 3.035682e-05
 #Now I'm almost done, just need to add some meaningful 
 #row names that were lost on the way.
 row.names(Ind_results_relative_mat) <- Ind_result_mat$Sample
 
 # I need to order the samples (rows)
 
 row.names(Ind_results_relative_mat)
 
 Ind_sam_order<-c("1A-rep1", "1B-rep1","1C-rep1","2A-rep1", "2B-rep1", "2C-rep1",
             "3A-rep1", "3B-rep1", "3C-rep1", "4A-rep1","4B-rep1", "4C-rep1",
               "5A-rep1", "5B-rep1","5C-rep1", "6A-rep1", "6B-rep1","6C-rep1",
              "7A-rep1", "7B-rep1", "7C-rep1","8A-rep1","8B-rep1", "8C-rep1") 
 
 
 #-------------- Data exploration ------------------# 
 
 # I should have all the results in one matrix and I can do 
 #the first ordination plots and heatmaps. I will use packages vegan and gplots, 
 #that I loaded earlier.
 
 #A basic NMDS plot is a good way to start checking that everything 
 # is ok with the data
 
 plot(metaMDS(Ind_results_relative_mat), type="text", display="sites")
 #Warning message:
 # In metaMDS(results_relative_mat) :
 #  stress is (nearly) zero: you may have insufficient data
 
 #Let's remember that Sample 1 was taken from a pristine site. 
 Ind_results_relative_matII<-Ind_results_relative_mat[-c(1:3),]
 
 plot(metaMDS(Ind_results_relative_mat[-c(1:3),]), type="text", display="sites")

 plot(metaMDS(Ind_results_relative_matII), type="text", display="sites")
 #This looks better.
 
 
 #At this point importing some metadata might be a good idea. 
 I_metadata <- read.table("Ind_metadata.txt", header=TRUE, stringsAsFactors = TRUE, sep="\t")
 
 #The samples were sent for qPCR array analysis with names simplified names.
 #Samples 1A, 1B and 1C were taken from the "start" or the river ("pristine").
 
 str(I_metadata)
 #'data.frame':	24 obs. of  3 variables:
 #$ sample     : Factor w/ 24 levels "1A","1B","1C",..: 1 2 3 4 5 6 7 8 9 10 ...
 #$ name       : Factor w/ 8 levels "Cattle Farm",..: 8 8 8 1 1 1 2 2 2 7 ...
 #$ description: Factor w/ 1 level "Surface river": 1 1 1 1 1 1 1 1 1 1 ...
 
 #Let's make an color palette (I'll save it for later).
 #Ind_col_pal_for_sites<-c("#9fd1e9", "#595441", "#9096f0","#553bdd","#001038", "#8d9fa6",  "#00c8ff",  "#005bff")
 
 IndHMxcol<-c(rep("springgreen4", 3), rep("sienna3", 3), rep("chocolate4",3), rep("mediumorchid1",3), 
              rep("mediumpurple3",3), rep("purple4",3), rep("royalblue1",3), rep("blue4",3))
 # adding colors like this will work although you would not know the order.
 levels(I_metadata$name)
 I_metadata$name.f<-factor(I_metadata$name, levels=c("Spring Water","Cattle Farm",
 "Chicken Slaughterhouse", "Hospital","City","City Downstream", "Estuary-Freshwater",
 "Estuary-Seawater"))
 levels(I_metadata$name.f)          
 I_metadata$color<-rep(1, nrow( I_metadata))  
 I_metadata<- within(I_metadata, color[name=="Spring Water"]<-"springgreen4")
 I_metadata<- within(I_metadata, color[name=="Cattle Farm"]<-"sienna3")
 I_metadata<- within(I_metadata, color[name=="Chicken Slaughterhouse"]<-"chocolate4")
 I_metadata<- within(I_metadata, color[name=="Hospital"]<-"mediumorchid1")
 I_metadata<- within(I_metadata, color[name=="City"]<-"mediumpurple3")
 I_metadata<- within(I_metadata, color[name=="City Downstream"]<-"purple4")
 I_metadata<- within(I_metadata, color[name=="Estuary-Freshwater"]<-"royalblue1")
 I_metadata<- within(I_metadata, color[name=="Estuary-Seawater"]<-"blue4")
 
 #I need to drop samples 1A-1C from metadata because I can't use them in NMDS
 I_metadataII<-subset(I_metadata, name!="Spring Water")
 levels(I_metadataII$name.f)
 # The order of colors will not work if the unused levels are not dropped:
 I_metadataII <- droplevels(I_metadataII)
 levels(I_metadataII$name.f)
 str(I_metadataII)
 
 IND_NMDS<-metaMDS(Ind_results_relative_matII)
 
 ordiplot(IND_NMDS,type="n", xlim=c(-2.6,2.3),
          ylim=c(-1,1),cex.axis = 1.5, cex.lab = 1.5, 
          main = "ARGs and MGEs")
with(I_metadataII, ordiellipse(IND_NMDS, name.f, 
 kind = "se", conf = 0.95, col=c("sienna3", "chocolate4", "mediumorchid1",
  "mediumpurple3", "purple4","royalblue1","blue4")))
 #with(I_metadataII, ordispider(IND_NMDS, name,col=c( "sienna3", "chocolate4", "mediumorchid1",
#"mediumpurple3", "purple4","royalblue1","blue4"), label=FALSE))
 with(I_metadataII, points(IND_NMDS$points,pch=17, cex=1.4, col=I_metadataII$color))
 with(I_metadataII, legend(-2.8,-0.7, legend= levels(name.f), cex=0.9, bty= "n", 
 col=c("sienna3", "chocolate4", "mediumorchid1","mediumpurple3",
   "purple4","royalblue1","blue4"), pch=17))
 
 
 adonis(Ind_results_relative_matII ~ name.f, data=I_metadataII, permutations=9999)
 
 #Call:
  # adonis(formula = Ind_results_relative_matII ~ name.f, data = I_metadataII,      permutations = 9999) 
 
 #Permutation: free
 #Number of permutations: 9999
 
 #Terms added sequentially (first to last)
 
            #Df SumsOfSqs MeanSqs F.Model R2    Pr(>F)    
 #name.f     6    2.9337 0.48896   4.699 0.6682  1e-04 ***
# Residuals 14    1.4568 0.10406         0.3318           
 #Total     20    4.3905                 1.0000           
 #---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
 levels(I_metadataII$name.f)
 
#[1] "Cattle Farm"            "Chicken Slaughterhouse" "Hospital"              
# [4] "City"                   "City Downstream"        "Estuary-Freshwater"    
#[7] "Estuary-Seawater"  
 
# 66.82 of the variation is explained by the sampling site, that's a lot  
 
 pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
 {
   library(vegan)
   co = as.matrix(combn(unique(factors),2))
   pairs = c()
   F.Model =c()
   R2 = c()
   p.value = c()
   
   for(elem in 1:ncol(co)){
     ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                   factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
     pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
     F.Model =c(F.Model,ad$aov.tab[1,4]);
     R2 = c(R2,ad$aov.tab[1,5]);
     p.value = c(p.value,ad$aov.tab[1,6])
   }
   p.adjusted = p.adjust(p.value,method=p.adjust.m)
   pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
   return(pairw.res)
 }
 #*https://github.com/pmartinezarbizu/pairwiseAdonis, 
 #original function “published” in a ResearchGate thread in 2015 
 #(https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R)
 

 pairwise.adonis(Ind_results_relative_matII, I_metadataII$name.f, sim.method = "bray", p.adjust.m = "BH")
 #p.adjust.m = "BH" is False discovery rate
 
 # Large p-values due to small group size....
 # all not significant 
 
 #And heatmap is also a good tool for understanding data. 
 #Some transformation is probably needed, here we used square root (sqrt()).
 #You can also do it without transormation. Depends on your data.
 
 #with transformation
 heatmap.2(as.matrix(sqrt(Ind_results_relative_mat)), col = rev(heat.colors(100)), 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 #without
 heatmap.2(as.matrix(Ind_results_relative_mat), col = rev(heat.colors(100)), 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 #Heatmap from only the most abundant genes is usually more helpful. 
 #For example the ones where the maximum abundance is over 0.01. 
 Ind_result_abund001_mat <- Ind_results_relative_mat[,apply(Ind_results_relative_mat, 2, max)>0.01]
 Ind_result_abund0001_mat <- Ind_results_relative_mat[,apply(Ind_results_relative_mat, 2, max)>0.001]
 
 heatmap.2(as.matrix(sqrt(Ind_result_abund001_mat)), col = rev(heat.colors(100)), 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 heatmap.2(as.matrix(sqrt(Ind_result_abund0001_mat)), col = rev(heat.colors(100)), #Colv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 #Lets put the samples in order and forget the clustering:
 heatmap.2(as.matrix(sqrt(Ind_result_abund0001_mat)), col = rev(heat.colors(100)),
           dendrogram="none",Rowv=FALSE, Colv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 heatmap.2(as.matrix(sqrt(Ind_result_abund001_mat)), col = rev(heat.colors(100)),
           dendrogram="none",Rowv=FALSE, Colv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 #clustering by rows and columns separately:
 
 heatmap.2(as.matrix(sqrt(Ind_result_abund001_mat)), col = rev(heat.colors(100)),
           dendrogram="row",Colv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 heatmap.2(as.matrix(sqrt(Ind_result_abund001_mat)), col = rev(heat.colors(100)),
           dendrogram="column",Rowv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 # We can also make a heatmap with ggplot2. 
 # The data needs to be in a long format.
 Ind_result_abund001_mat.m<-melt(Ind_result_abund001_mat)
 
 #With ggplot2 the colors work better with NA's. 
 Ind_result_abund001_mat.m$valueII<-Ind_result_abund001_mat.m$value
 min(Ind_result_abund001_mat.m$valueII)
 Ind_result_abund001_mat.m$valueII[Ind_result_abund001_mat.m$valueII==9.536743e-07]<-NA
 
 #lets check what genes are among the most abundant ones:
 
 levels(Ind_result_abund001_mat.m$Var2)
 #[1] "ISPps"        "ISSm2"        "IncP_1alpha"  "aadA_1"       "aadA_2"       "aadA_3"       "aadA_4"       "aadA_5"       "aadA_6"       "blaACT_3"    
 #[11] "blaGES"       "blaOXA_1"     "blaOXA_2"     "cmlA_2"       "cmlA_3"       "cmlA_4"       "floR_1"       "intI1_1"      "intI1_2"      "intI1_3"     
 #[21] "intI3_1"      "intI3_2"      "mexF"         "qacEdelta1_1" "qacEdelta1_2" "qacEdelta1_3" "qacH_2"       "rpoB"         "sul1_2"       "sul2_1"      
 #[31] "sul2_2"       "tetA_2"       "tetG_1"       "tetG_2"       "tetM_1"       "tetR_2"       "tnpA-02/IS4"  "tnpA-01/IS21"
  
 #Then lets order the assays we have according to the antibiotic they give resistance
 Gene_ord<-c( "aadA_1","aadA_2", "aadA_3", "aadA_4", "aadA_5" , "aadA_6", "blaACT_3",
              "blaGES","blaOXA_1", "blaOXA_2","qacEdelta1_1", "qacEdelta1_2","qacEdelta1_3",
              "qacH_2","mexF", "cmlA_2", "cmlA_3","cmlA_4","floR_1","sul1_2", "sul2_1",      
               "sul2_2","tetA_2","tetG_1","tetG_2","tetM_1","tetR_2","IncP_1alpha" ,"intI1_1",
              "intI1_2", "intI1_3","intI3_1","intI3_2","ISPps", "ISSm2","tnpA-01/IS21","tnpA-02/IS4", "rpoB")  
 
 Ind_result_abund001_mat.m$Genes<-factor(Ind_result_abund001_mat.m$Var2, levels =Gene_ord)
 # rpoB is a housekeeping gene
 toBeRemoved<-which(Ind_result_abund001_mat.m$Var2=="rpoB")
 Ind_result_abund001_mat.m<-Ind_result_abund001_mat.m[-toBeRemoved,]
 Ind_result_abund001_mat.m<-droplevels(Ind_result_abund001_mat.m)
# I want to get rid of the "-rep1" in the sample name 
  Ind_result_abund001_mat.m$Sample<- substr(Ind_result_abund001_mat.m$Var1, 1,2)
  
 
 #and give colors to our sampling sites
  
# IndHMxcol<-c(rep("#9fd1e9", 3), rep("#595441", 3), rep("#9096f0",3), rep("#553bdd",3), 
           #   rep("#001038",3), rep("#8d9fa6",3), rep("#00c8ff",3), rep("#005bff",3))
# topo.colors(8)
 #IndHMxcol<-rev(c(rep("#4C00FFFF", 3), rep("#0019FFFF", 3), rep("#0080FFFF",3), rep( "#00E5FFFF",3), 
            #  rep("#00FF4DFF",3), rep("#E6FF00FF",3), rep("#FFFF00FF",3), rep("#FFE0B3FF",3)))
 
 IndHMxcol<-c(rep("springgreen4", 3), rep("sienna3", 3), rep("chocolate4",3), rep("mediumorchid1",3), 
              rep("mediumpurple3",3), rep("purple4",3), rep("royalblue1",3), rep("blue4",3))
 
 #and get rid of the assay number in the gene name...
 #result_abund001_mat.mO$Gene_name<-result_abund001_mat.mO$Genes
 #result_abund001_mat.mO<- separate(data = result_abund001_mat.mO, col = Gene_name, into = c("Gene_name", "AssayNO"), sep = "_")
 
 #For the color breaks:  
 #min(subset(result_abund001_mat.m$valueII, result_abund001_mat.m$valueII!="NA"))
 #[1] 1.220662e-05
 #max(subset(result_abund001_mat.m$valueII, result_abund001_mat.m$valueII!="NA"))
 #[1] 0.3530962
 #https://www.rdocumentation.org/packages/colorspace/versions/1.4-1/topics/scale_colour_continuous_sequential
 
 
 Ind_HM001 <- ggplot(Ind_result_abund001_mat.m, aes(x=Sample, y=Genes, fill=valueII))+
   geom_tile(color = "black", size=0.005)+
  scale_fill_continuous_sequential(palette =  "Oslo",
  #scale_fill_continuous_sequential(palette =  "Dark Mint",
   begin = 0.9, end = 0, na.value="white", name="Relative abundance")+
   #scale_fill_continuous_sequential(palette =  "Lajolla",
   #I don't want the most darkest color and I want to reverse the colors
   #begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
   theme(panel.background=element_rect(fill="black", colour="black")) +
   theme(axis.text.y=element_text(colour= "black",size=10))+
   theme(axis.text.x = element_text(colour= "black",angle = 90, hjust=0, vjust=0, size=10))+
   theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
   theme(panel.border=element_blank())+
   theme(axis.title.x = element_blank()) + 
   theme(axis.title.y = element_blank()) + 
   theme(legend.position="bottom")+
   theme(legend.key.size=unit(0.2, "cm"))+
   theme(legend.key.width=unit(0.5, "cm"))
 Ind_HM001
 
 
 # I'll see if I could include more genes
 
 max(Ind_result_abund001_mat.m$value)
 #[1] 0.4622248
 #With "Oslo"
 #   begin = 0.8, end = 0, na.value="white", name="Relative abundance")+
 
 
 
 Ind_result_abund0001_mat.m<-melt( Ind_result_abund0001_mat)
 
 Ind_result_abund0001_mat.m$valueII<- Ind_result_abund0001_mat.m$value
 
 Ind_result_abund0001_mat.m$valueII[Ind_result_abund0001_mat.m$valueII==9.536743e-07]<-NA
 
 #lets check what genes are among the most abundant ones:
 
 levels(Ind_result_abund0001_mat.m$Var2)
 
 Gene_ordII<-c("aac(6)-Ib_1","aac(6)-Ib_2" , "aac(6)-Ib_3", "aadA_1", "aadA_2", "aadA_3", "aadA_4",         
    "aadA_5","aadA_6","aadA5_1", "aadA5_2","aadD","aadE", "strB",
    "blaACT_2", "blaACT_3" , "blaCMY_3", "blaDHA", "blaGES","blaOXA_1" ,"blaOXA_2", "blaOXA_3",       
   "blaTLA" ,"blaVEB","cfxA",
   "qacEdelta1_1", "qacEdelta1_2", "qacEdelta1_3" ,"qacH_1" , "qacH_2",
   "ermB", "ermC", "ermF","ermX",  "lnuB_1", "lnuB_2",  
   "ereA", "matA_mel" ,"merA","mexE","mexF",      
   "catB3" ,"cmlA_2", "cmlA_3" , "cmlA_4",         
   "cmxA","floR_1",             
   "sul1_2" ,"sul2_1", "sul2_2",
   "tet(32)" ,"tetA_2","tetC_3","tetG_1" ,"tetG_2","tetH" ,          
  "tetL_2" ,"tetM_1", "tetM_2","tetM_3" ,"tetO_1", "tetO_2","tetQ",           
   "tetR_2","tetW", "tetX" ,     
   "dfrA1_1",
  "IncP_1alpha" , "intI1_1", "intI1_2","intI1_3","intI1_4", "intI3_1",        
  "intI3_2", "IS613", "ISAba3","ISPps", "ISSm2",                
  "repA",  "tnpA-02/IS4" ,"tnpA-01/IS21" ,"tnpA-07/ISEcp1B", "tnpA-06/IS1216", 
 "tnpA-03/IS6" ,
 "rpoB") 

 
 Ind_result_abund0001_mat.m$Genes<-factor(Ind_result_abund0001_mat.m$Var2, levels =Gene_ordII)
 
 toBeRemovedII<-which(Ind_result_abund0001_mat.m$Var2=="rpoB")
 
 Ind_result_abund0001_mat.m<-Ind_result_abund0001_mat.m[-toBeRemovedII,]
 Ind_result_abund0001_mat.m<-droplevels(Ind_result_abund0001_mat.m)
 # I want to get rid of the "-rep1" in the sample name 
 Ind_result_abund0001_mat.m$Sample<- substr(Ind_result_abund0001_mat.m$Var1, 1,2)
 
 
length( Gene_ordII)
 
 library(colorspace)
 
 summary(Ind_result_abund0001_mat.m$valueII)
 
 
 Ind_HM0001 <- ggplot(Ind_result_abund0001_mat.m, aes(x=Sample, y=Genes, fill=valueII))+
   geom_tile(color = "black", size=0.005)+
   scale_fill_continuous_sequential(palette =  "Oslo",
                                    #scale_fill_continuous_sequential(palette =  "Dark Mint",
                                    begin = 0.8, end = 0, na.value="white", name="Relative abundance")+
   #scale_fill_continuous_sequential(palette =  "Lajolla",
   #I don't want the most darkest color and I want to reverse the colors
   #begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
   theme(panel.background=element_rect(fill="black", colour="black")) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black",angle = 90, hjust=0, vjust=0, size=10))+
   theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
   theme(panel.border=element_blank())+
   theme(axis.title.x = element_blank()) + 
   theme(axis.title.y = element_blank()) + 
   theme(legend.position="bottom")+
   theme(legend.key.size=unit(0.2, "cm"))+
   theme(legend.key.width=unit(0.5, "cm"))
 Ind_HM0001
 
 
 # I would like to scale the colors:
 
 dat <- data.frame(values = as.numeric(Ind_result_abund0001_mat.m$valueII))
 ggplot(dat, aes(values)) + geom_density(bw = "SJ")  
 
 range(Ind_result_abund0001_mat.m$value)
 
 summary(Ind_result_abund0001_mat.m$valueII)
 
 min(subset(Ind_result_abund001_mat.m$valueII, Ind_result_abund001_mat.m$valueII!="NA"))
 #[1] 8.080552e-05 -> 0.00008080552
 #1.227083e-05 -> 0.0000122708
 # max: 4.622248e-01 -> 0.5
 
 Ind_result_abund0001_mat.m$categories <- Ind_result_abund0001_mat.m%>% .$value %>% cut(c(9.536743e-07,1.227083e-05, 8.180552e-05,
     1.227083e-04,8.180552e-04, 1.227083e-03,8.180552e-03,  1.227083e-02, 8.180552e-02, 1.227083e-01, 8.180552e-01))
 b <- Ind_result_abund0001_mat.m%>% .$value %>% cut(c(9.536743e-07,  1.227083e-05, 8.180552e-05,
  1.227083e-04,8.180552e-04, 1.227083e-03,8.180552e-03,  1.227083e-02, 8.180552e-02, 1.227083e-01, 8.180552e-01))
 levels(b) <- 1:11
 Ind_result_abund0001_mat.m$categories2 <- b
 
 Indheatcol<-heat.colors(10)
 
 IndHMcolors<- sequential_hcl(n=10, palette="Mako")
 
 
 IndpizzaHM<-ggplot(Ind_result_abund0001_mat.m,aes(x=Sample, y=Genes, fill= categories))+
    geom_tile(aes(fill=categories2))+
    scale_fill_manual(values= rev(IndHMcolors), na.value="white", name= "Relative Abundance") +
    theme(axis.text.x = element_text(colour= "black", angle = 90, hjust = 1, vjust = 0.5, size=7), axis.text.y= element_text(colour="black", size=7))+
    theme(axis.title.x = element_blank(), axis.title.y =  element_blank())+
    theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
    theme(panel.border=element_blank())+
    theme(axis.title.x = element_blank()) + 
    theme(axis.title.y = element_blank()) + 
    #theme(legend.position="bottom")+
    theme(legend.key.size=unit(0.2, "cm"))+
    theme(legend.key.width=unit(0.5, "cm"))
 
 IndpizzaHM
 
 # the color key will need to be edited in InkScape
 
 
 ####################
 
 # What genes are detected in only one treatment group and how many genes are shared?
 # venn diagrams
 
 View(Ind_results_relative_mat)
 Ind_results_relative_mat.m <- melt(Ind_results_relative_mat)
 Ind_results_relative_mat.m$Sample<- substr(Ind_results_relative_mat.m$Var1, 1,2)
 
 Ind_results_relative_mat.m$S.type<-rep(1, nrow(Ind_results_relative_mat.m))  
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="1A"| Sample=="1B"| Sample=="1C" ]<-"Spring Water")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="2A"| Sample=="2B"| Sample=="2C" ]<-"Rural")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="3A"| Sample=="3B"| Sample=="3C" ]<-"Rural")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="4A"| Sample=="4B"| Sample=="4C" ]<-"City")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="5A"| Sample=="5B"| Sample=="5C" ]<-"City")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="6A"| Sample=="6B"| Sample=="6C" ]<-"City")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="7A"| Sample=="7B"| Sample=="7C" ]<-"Estuary")
 Ind_results_relative_mat.m<- within(Ind_results_relative_mat.m, S.type[Sample=="8A"| Sample=="8B"| Sample=="8C" ]<-"Estuary")
 
 #Lets convert the data frame back into matrix:
 
 Ind_results_relative_mat_venn <- dcast(Ind_results_relative_mat.m, S.type+Sample~Var2, value.var="value")
 rownames(Ind_results_relative_mat_venn)<-interaction(Ind_results_relative_mat_venn$S.type,Ind_results_relative_mat_venn$Sample, sep="_")
 Ind_results_relative_mat_venn<-as.matrix(Ind_results_relative_mat_venn[,3:ncol(Ind_results_relative_mat_venn)])
 
 # No genes detected in Spring Water, so
 Ru_genes<-as.matrix(Ind_results_relative_mat_venn[grep("Rural",rownames(Ind_results_relative_mat_venn)),]
                     [, which(colSums(Ind_results_relative_mat_venn[grep("Rural",rownames(Ind_results_relative_mat_venn)),]) != 6*(9.536743e-07))])
 str(Ru_genes)
 Ru_genesU<-colnames(Ru_genes)
 str(Ru_genesU)
 length(Ru_genesU)
 #[1] 159
 
 City_genes<-as.matrix(Ind_results_relative_mat_venn[grep("City",rownames(Ind_results_relative_mat_venn)),]
                       [, which(colSums(Ind_results_relative_mat_venn[grep("City",rownames(Ind_results_relative_mat_venn)),]) != 9*(9.536743e-07))])
 str(City_genes)
 City_genesU<-colnames(City_genes)
 str(City_genesU)
 length(City_genesU)
 #[1] 164
 
 Est_genes<-as.matrix(Ind_results_relative_mat_venn[grep("Est",rownames(Ind_results_relative_mat_venn)),]
                      [, which(colSums(Ind_results_relative_mat_venn[grep("Est",rownames(Ind_results_relative_mat_venn)),]) != 6*(9.536743e-07))])
 str(Est_genes)
 Est_genesU<-colnames(Est_genes)
 str(Est_genesU)
 length(Est_genesU)
 #[1] 41
 
 input<-list(Rural=Ru_genesU, City=City_genesU, Estuary= Est_genesU)
 
 venn(input)
 
 universe <- unique(c(Ru_genesU,City_genesU,Est_genesU))
 Ru.l <- universe %in% Ru_genesU
 City.l <- universe %in% City_genesU
 Est.l <- universe %in% Est_genesU
 
 universe[Ru.l & !City.l &!Est.l]
 #[1] "ISEfm1"      "aacA_aphD"   "aacC4"       "acrR_1"      "blaKPC_1"    "erm(36)"    
 #[7] "ermC"        "intI2_1"     "lmrA_1"      "nisB_2"      "pbp"         "qnrB"       
 #[13] "str"         "strA"        "tetJ"        "tetPB_3"     "tolC_3"      "vanC2_vanC3"
 #[19] "vanC_1"      "vanC_4"      "vanSC_2"     "vanTC_1"     "vanTC_3"
 
 universe[!Ru.l & City.l & !Est.l]
 #[1] "IncP_1beta" "IncP_oriT"  "IncQ_oriT"  "aac(6)-II"  "aadA9_1"    "acrA_5"    
 #[7] "ampC_5"     "blaACT_2"   "blaIMP_3"   "blaSFO"     "blaSHV_2"   "blaTLA"    
 #[13] "catA1"      "dfrA12"     "ereB"       "erm(35)"    "mdtE"       "oprD"      
 #[19] "oprJ"       "qnrA"       "sul3"       "tet(36)_1"  "tetA(P)"    "tetM_3"    
 #[25] "tetT"       "vanYD_1"   
 
 universe[!Ru.l & !City.l & Est.l]
 #character(0)
 
 universe[!Ru.l & City.l & Est.l]
 #[1] "blaOXA_2" "repA"    
 
 universe[Ru.l & City.l & !Est.l]
 #[1] "IS613"         "ISAba3"        "IncN_rep"      "Tp614"         "aac(6)-Ib_1"  
 #[6] "aac(6)-Ib_2"   "aac(6)-Ib_3"   "aacC2"         "aadA5_1"       "aadA5_2"      
 #[11] "aadD"          "aadE"          "acrA_1"        "acrA_2"        "acrA_4"       
 #[16] "acrB_1"        "acrF"          "acrR_2"        "acrR_3"        "ampC_1"       
 #[21] "ampC_2"        "aphA1_7"       "aphA3_1"       "aphA3_2"       "bacA_2"       
 #[26] "blaACT_3"      "blaCMY_1"      "blaCMY_2"      "blaCMY_3"      "blaDHA"       
 #[31] "blaFOX"        "blaMOX_blaCMY" "blaOXA_3"      "blaOXA_4"      "blaPSE"       
 #[36] "blaSHV_1"      "blaTEM"        "blaVEB"        "catB3"         "catB8"        
 #[41] "cfxA"          "cmlA_3"        "cphA_1"        "cphA_2"        "dfrA1_1"      
 #[46] "dfrA1_2"       "ermB"          "ermF"          "ermX"          "gapA"         
 #[51] "lnuA_1"        "lnuB_1"        "lnuB_2"        "matA_mel"      "mdh"          
 #[56] "mdtA"          "mdtF"          "mdtG_1"        "mdtH_1"        "mdtH_2"       
 #[61] "mdtL"          "mefA"          "merA"          "mphA_1"        "mphA_2"       
 #[66] "mtrD_3"        "qacH_1"        "rarD_2"        "sat4"          "tet(32)"      
 #[71] "tetA_B_1"      "tetA_B_2"      "tetC_1"        "tetC_2"        "tetC_3"       
 #[76] "tetE"          "tetH"          "tetL_2"        "tetM_1"        "tetM_2"       
 #[81] "tetO_1"        "tetO_2"        "tetQ"          "tetR_3"        "tetS"         
 #[86] "tetW"          "tetX"          "tnpA_3"        "tnpA_4"        "tnpA_5"       
 #[91] "tnpA_6"        "tnpA_7"
 

 #######################################################################
 
######### 16S rRNA sequence analysis #################
 
  #Let's bring in OTUs!
 
 #####################################################################
 
 # this will be a large file.... 
 
 Ind_OTUtableTEMP_I <- read.table("otutable_indonesia.txt", header=TRUE, sep="\t")
 str(Ind_OTUtableTEMP_I)
 
 grep("Eukaryota",Ind_OTUtableTEMP_I$taxonomy)
 
 grep("Mitochondria",Ind_OTUtableTEMP_I$taxonomy)
 
 grep("Chloroplast",Ind_OTUtableTEMP_I$taxonomy)
 
 grep("unknown",Ind_OTUtableTEMP_I$taxonomy)
 
 grep("Archaea",Ind_OTUtableTEMP_I$taxonomy)
 
 EukarOTUrows<-grep("Eukaryota",Ind_OTUtableTEMP_I$taxonomy)
 
 MitochOTUrows<-grep("Mitochondria",Ind_OTUtableTEMP_I$taxonomy)
 
 ChloroplastOTUrows<-grep("Chloroplast",Ind_OTUtableTEMP_I$taxonomy)
 
 unknownOTUrows<-grep("unknown",Ind_OTUtableTEMP_I$taxonomy)
 
 ArchaeaOTUrows<- grep("Archaea",Ind_OTUtableTEMP_I$taxonomy)
 
 OTUstoBeRemoved<-c(EukarOTUrows, MitochOTUrows,ChloroplastOTUrows,unknownOTUrows, ArchaeaOTUrows)
 
 Ind_OTUtableTEMP_II<-Ind_OTUtableTEMP_I[-OTUstoBeRemoved,]
 
 dim(Ind_OTUtableTEMP_II)
 dim(Ind_OTUtableTEMP_I)
 
 # ... so I will reduce its size as soon & as much as possible
 
 Ind_OTUtableTEMP<-Ind_OTUtableTEMP_II
 rm(Ind_OTUtableTEMP_II)
 rm(Ind_OTUtableTEMP_I)
 
 # Above I removed everything that are not bacteria
 
 # Now I'll continue working on the file. For instance, 
 # there are samples that do not belong to this dataset
 
 row.names(Ind_OTUtableTEMP)<-Ind_OTUtableTEMP$OTU.ID
 Ind_OTUtableTEMP$OTU.ID<-NULL
 
 IndTaxonomy<- as.data.frame(Ind_OTUtableTEMP$taxonomy)
 row.names(IndTaxonomy)<- row.names(Ind_OTUtableTEMP)
 
 Ind_OTUtableTEMP$taxonomy<-NULL
 Ind_OTUtableTEMP<-Ind_OTUtableTEMP[grep("w",colnames( Ind_OTUtableTEMP))]
 dim(Ind_OTUtableTEMP) 
 
 Ind_OTUtableF<-t(Ind_OTUtableTEMP)
 str(Ind_OTUtableF)
 
 # I also want to get rid of singletons and doubletons
 
 dim(Ind_OTUtableF)
 Ind_OTUtableF<-Ind_OTUtableF[, which(colSums(Ind_OTUtableF)!= 0)]
 dim(Ind_OTUtableF)
 plot(metaMDS(Ind_OTUtableF), type="text", display="sites")
 
 Ind_OTUtableFover3<-Ind_OTUtableF[, which(colSums(Ind_OTUtableF)>3)]
 dim(Ind_OTUtableFover3)
 plot(metaMDS(Ind_OTUtableFover3), type="text", display="sites")
 
 rownames(Ind_OTUtableFover3)
 
 include_list<-c("A007.1C.w1.CGTAATAG", "A008.1C.w2.TCTAAGGA", "A009.1C.w3.CGAACTTC",
                 "A010.2A.w1.TAGCACGC", "A011.2A.w2.AGTCCGAA", "A012.2A.w3.CATACTCT",
                 "A013.2B.w1.CGATACGC","A014.2B.w2.TGTTCTCG", "A015.2B.w3.TTCCTCCG",
                 "A016.2C.w1.TCTCTATC", "A017.2C.w2.CCTGTATG", "A018.2C.w3.ATCCGACG",
                 "A019.3A.w1.GGTCCATT","A020.3A.w2.CCAGTCAC", "A021.3A.w3.AGACGAAC",
                 "A028.3B.w1.AGTTCCGC", "A029.3B.w2.CCGCCTAA", "A030.3B.w3.GTAAGCCT",
                 "A034.4B.w1.AGTCAGTT","A035.4B.w2.TACATCTG","A036.4B.w3.GGTAGAAT",
                 "A037.4C.w1.TCCACCTT","A038.4C.w2.GAGGATAT","A039.4C.w3.GAGATTGT")
 
 Ind_OTUtableFover3F2<-subset(Ind_OTUtableFover3, rownames(Ind_OTUtableFover3) %in% include_list)
 
 dim(Ind_OTUtableFover3F2)
 
 #Now I have only those samples that belong to this dataset.
 
 # I will test if this still works
 
 plot(metaMDS(Ind_OTUtableFover3F2), type="text", display="sites")
 
 # TSS-normalization:
 
 Ind_OTUtableFover3F2_rel <- apply(Ind_OTUtableFover3F2, 1, function(i) i/sum(i))
 
 Ind_OTUtableFover3F2_rel<-t(Ind_OTUtableFover3F2_rel)
 
 plot(metaMDS(Ind_OTUtableFover3F2_rel), type="text", display="sites")
 
 # The OTU data has its own metadata:
 Ind_OTU_Sample_names<- read.table("Otu_SampleNamesJM.txt", header=TRUE, sep="\t")
 
 Ind_OTU_Sample_names$Seq_name
 
 Ind_OTUtableFover3F2_rel_o<-Ind_OTUtableFover3F2_rel[match(Ind_OTU_Sample_names$Seq_name,row.names(Ind_OTUtableFover3F2_rel)),]
 
 row.names(Ind_OTUtableFover3F2_rel_o)
 
 row.names(Ind_OTUtableFover3F2_rel_o)<-Ind_OTU_Sample_names$New_nameII
 
 plot(metaMDS(Ind_OTUtableFover3F2_rel_o), type="text", display="sites")
 

 View(I_metadataII)
 View(I_metadata)
 levels(I_metadata$name.f)
 str(I_metadataII)
 str(I_metadata)
 
 
 # I want to replace zeros:
 min(Ind_OTUtableFover3F2_rel_o)
 min(Ind_OTUtableFover3F2_rel_o[Ind_OTUtableFover3F2_rel_o!=min(Ind_OTUtableFover3F2_rel_o)])
 #[1] 4.748564e-06
 
 Ind_OTUtableFover3F2_rel_o[Ind_OTUtableFover3F2_rel_o == 0] <- 9.536743e-08
 
 IND_OTU_NMDS<-(metaMDS(Ind_OTUtableFover3F2_rel_o))
 
 ordiplot(IND_OTU_NMDS,type="n", xlim=c(-2.6,2.3),
          ylim=c(-1,1),cex.axis = 1.5, cex.lab = 1.5, 
          main = "OTUs")
 with(I_metadata, ordiellipse(IND_OTU_NMDS, name.f, 
                              kind = "se", conf = 0.95, col=c("springgreen4","sienna3", "chocolate4", "mediumorchid1",
                                                              "mediumpurple3", "purple4","royalblue1","blue4")))
 with(I_metadata, points(IND_OTU_NMDS$points,pch=19, cex=1.4, col=I_metadata$color))
 with(I_metadata, legend(-2.82,1.5, legend= levels(name.f), cex=0.8, bty= "n", 
                         col=c("springgreen4","sienna3", "chocolate4", "mediumorchid1","mediumpurple3",
                               "purple4","royalblue1","blue4"), pch=19))
 
 
 # PERMANOVA
 
 adonis(Ind_OTUtableFover3F2_rel_o ~ name.f, data=I_metadata, permutations=9999)
 
 #Call:
 #  adonis(formula = Ind_OTUtableFover3F2_rel_o ~ name.f, data = I_metadata,      permutations = 9999) 
 
 #Permutation: free
 #Number of permutations: 9999
 
 #Terms added sequentially (first to last)
 
 #Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
 #name.f     7    5.7829 0.82613  20.121 0.89799 1e-04 ***
 # Residuals 16    0.6497 0.04106          0.10201            
 #Total     23    6.4398                 1.00000           
 #---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
 
 pairwise.adonis(Ind_OTUtableFover3F2_rel_o, I_metadata$name.f, sim.method = "bray", p.adjust.m = "BH")
 
 
 # I also want to test if the library size gives me trouble:
 #meta_noMCnmdsII$lib_size<-rowSums(noMC_phylOTUtable_mat)
 
 rowSums(Ind_OTUtableFover3F2)
 
 View(Ind_OTUtableFover3F2)
 
 Ind_OTUtableFover3F2_o<-Ind_OTUtableFover3F2[match(Ind_OTU_Sample_names$Seq_name,row.names(Ind_OTUtableFover3F2)),]
 
 row.names(Ind_OTUtableFover3F2_o)
 
 row.names(Ind_OTUtableFover3F2_o)<-Ind_OTU_Sample_names$New_nameII
 
 rowSums(Ind_OTUtableFover3F2_o)
 
 I_metadata$lib_size <-rowSums(Ind_OTUtableFover3F2_o)
 
 
 #General adonis for relative abundance normalized OTUs
 #that has the also the variable library size:
 adonis(Ind_OTUtableFover3F2_rel_o ~ I_metadata$name.f + I_metadata$lib_size, 
        permutations = 9999)
 #Call:
 # adonis(formula = Ind_OTUtableFover3F2_rel_o ~ I_metadata$name.f +      I_metadata$lib_size, permutations = 9999) 
 
 #Permutation: free
 #Number of permutations: 9999
 
 #Terms added sequentially (first to last)
 
 #                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 #I_metadata$name.f    7    5.7927 0.82753 19.4449 0.89915 0.0001 ***
 #I_metadata$lib_size  1    0.0113 0.01135  0.2666 0.00176 0.9756    
 #Residuals           15    0.6384 0.04256         0.09909           
 #Total               23    6.4424                 1.00000           
 #---
 #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
 # Good!! But let's see this too:
 
 adonis(Ind_OTUtableFover3F2_rel_o ~  I_metadata$lib_size, 
        permutations = 9999)
 
 #Call:
 #adonis(formula = Ind_OTUtableFover3F2_rel_o ~ I_metadata$lib_size,      permutations = 9999) 
 
 #Permutation: free
 #Number of permutations: 9999
 
 #Terms added sequentially (first to last)
 #                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 #I_metadata$lib_size  1    1.4503 1.45029  6.3913 0.22434  1e-04 ***
 #Residuals           22    4.9921 0.22692         0.77488           
 #Total               23    6.4424                 1.00000           
 
 # Even as the only term the library size accounts for only 22% of the variance so
 # We are safe :)
 
 # However, lets check the correlation with library size and Shannon's diversity
 # I will use Shannon's diversity also later on.
 # Why Shannon's diversity? 
 
 #Irrespective of its faults, the Shannon index (H′) seems a useful general diversity index that 
 #is influenced by both richness and evenness and is more sensitive to changes in abundance of the rare groups. 
 #The meaning of the number is more comprehensible when expressed as the exponential. The underestimation caused by 
 #incomplete coverage can be minimised by analysing several hundred clones per sample and by only comparing 
 #samples of equal size. The Shannon evenness index (E) is, likewise, more sensitive to rare OTUs. 
 #As with H′ its meaning is clearer when the exponential (eH′/S) is used. Low coverage causes overestimation.
 
 #https://academic.oup.com/femsec/article/43/1/1/509981
 
 I_metadata$OTU_divH<-diversity(Ind_OTUtableFover3F2_rel_o)
 
 plot(I_metadata$OTU_divH, I_metadata$lib_size)
 
 cor(I_metadata$OTU_divH, I_metadata$lib_size)
 #[1] -0.4478421
 cor(I_metadata$OTU_divH, I_metadata$lib_size, method = "spearman")
 #[1] -0.386087
 cor.test(I_metadata$OTU_divH, I_metadata$lib_size, method = "spearman")
 
 #Spearman's rank correlation rho
 
 #data:  I_metadata$OTU_divH and I_metadata$lib_size
 #S = 3188, p-value = 0.06332
 #alternative hypothesis: true rho is not equal to 0
 #sample estimates:
 #     rho 
 #-0.386087 
 
 # So no relationship here, we can proceed :))))
 
 
 ###############################
 
 # OTU data exploration
 
 #Heatmaps:
 
 
 Ind_OTU_abund001_mat <- Ind_OTUtableFover3F2_rel_o[,apply(Ind_OTUtableFover3F2_rel_o, 2, max)>0.01]
 
 dim(Ind_OTU_abund001_mat)
 
 #Lets  forget the clustering:
 heatmap.2(as.matrix(sqrt(Ind_OTU_abund001_mat)), col = rev(heat.colors(100)),
           dendrogram="none",Rowv=FALSE, Colv=FALSE, 
           trace = "none", density.info = "none",  srtCol=45, margins=c(5,8))
 
 # heatmap with ggplot2. 
 # The data needs to be in a long format.
 Ind_OTU_abund001_mat.m<-melt(Ind_OTU_abund001_mat)
 
 #With ggplot2 the colors work better with NA's. 
 Ind_OTU_abund001_mat.m$valueII<-Ind_OTU_abund001_mat.m$value
 min(Ind_OTU_abund001_mat.m$valueII)
 Ind_OTU_abund001_mat.m$valueII[Ind_OTU_abund001_mat.m$valueII==9.536743e-08]<-NA
 
 #lets check what OTUs are among the most abundant ones:
 
 levels(Ind_OTU_abund001_mat.m$Var2)
 #[1] "OTU26508" "OTU28394" "OTU37232" "OTU881"   "OTU2459"  "OTU2506"  "OTU471"   "OTU25798" "OTU40731"
 #[10] "OTU25159" "OTU27550" "OTU464"   "OTU150"   "OTU6570"  "OTU6586"  "OTU33261" "OTU1699"  "OTU6121" 
 #[19] "OTU6119"  "OTU7973"  "OTU25160" "OTU34215" "OTU6120"  "OTU28499" "OTU10218" "OTU10237" "OTU454"  
 #[28] "OTU6117"  "OTU1508"  "OTU3006"  "OTU302"   "OTU457"   "OTU462"   "OTU10912" "OTU164"   "OTU75"   
 #[37] "OTU2422"  "OTU6122"  "OTU26027" "OTU28205" "OTU37233" "OTU4990"  "OTU6123"  "OTU6578"  "OTU6134" 
 #[46] "OTU460"   "OTU484"   "OTU5000"  "OTU455"   "OTU53"    "OTU480"   "OTU6166"  "OTU486"   "OTU180"  
 #[55] "OTU28396" "OTU28403" "OTU499"   "OTU466"   "OTU28400" "OTU1717"  "OTU477"   "OTU479"   "OTU2458" 
 #[64] "OTU463"   "OTU28431" "OTU468"   "OTU25664" "OTU26794" "OTU28405" "OTU476"   "OTU28402" "OTU26604"
 #[73] "OTU475"   "OTU503"   "OTU28398" "OTU565"   "OTU10715" "OTU21109"
 
 IndTaxonomy2<-IndTaxonomy 
 colnames(IndTaxonomy2)<-"Taxonomy"
 IndTaxonomy2<- separate(data = IndTaxonomy2, col = Taxonomy, 
                         into = c("kingdom", "phylum", "class","order", "family", "genus", "species"), sep = ";")
 
 IndTaxonomy2$OTU<-row.names(IndTaxonomy2)
 
 Ind_OTU_abund001_mat.mTAX<-merge(Ind_OTU_abund001_mat.m, IndTaxonomy2, by.x="Var2",by.y="OTU")
 
 # just checking:
 unique(IndTaxonomy2$kingdom)
 
 #IndHMxcol<-c(rep("springgreen4", 3), rep("sienna3", 3), rep("chocolate4",3), rep("mediumorchid1",3), 
 #rep("mediumpurple3",3), rep("purple4",3), rep("royalblue1",3), rep("blue4",3))
 
 Ind_OTU_HM001 <- ggplot(Ind_OTU_abund001_mat.mTAX,aes(x=Var1, y=genus, fill=valueII))+
   geom_tile()+
   scale_fill_continuous_sequential(palette =  "Oslo",
                                    #scale_fill_continuous_sequential(palette =  "Dark Mint",
                                    begin = 0.8, end = 0, na.value="white", name="Relative abundance")+
   #scale_fill_continuous_sequential(palette =  "Lajolla",
   #I don't want the most darkest color and I want to reverse the colors
   #begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
   theme(panel.background=element_rect(fill="black", colour="black")) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black",angle = 90, hjust=0, vjust=0, size=10))+
   theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=0, size=10))+
   theme(panel.border=element_blank())+
   theme(axis.title.x = element_blank()) + 
   theme(axis.title.y = element_blank()) + 
   theme(legend.position="bottom")+
   theme(legend.key.size=unit(0.2, "cm"))+
   theme(legend.key.width=unit(0.5, "cm"))
 Ind_OTU_HM001
 
 Ind_HM001
 
 
 #People usually like morte abut OTU barplots.
 # I'll make one.
 
 str(Ind_OTUtableFover3F2_rel_o)
 View(Ind_OTUtableFover3F2_rel_o)
 INDotu.summary <- prop.table(as.matrix(Ind_OTUtableFover3F2_rel_o), 1) 
 INDotu_abund <- colSums(INDotu.summary)
 INDotu.summary <- rbind(INDotu_abund, INDotu.summary)
 INDotu.summary_sorted <- INDotu.summary[,order(INDotu.summary[1,], decreasing = TRUE)]
 
 INDnum_genera <- 16 # enter the number of genera you want
 
 INDmelt_otu <- melt(INDotu.summary_sorted[,c(1:INDnum_genera)])
 colnames(INDmelt_otu) <- c("Sample", "OTU", "Abundance")
 tail(INDmelt_otu)
 
 
 #Ind_OTUtableFover3F2_rel_oBAR.m<-melt(Ind_OTUtableFover3F2_rel_o)
 #Ind_OTUtableFover3F2_rel_oBAR.mTAX<-merge(Ind_OTUtableFover3F2_rel_oBAR.m, IndTaxonomy2, by.x="Var2",by.y="OTU")
 
 
 IND_meta_otu <- merge(I_metadata, INDmelt_otu, by.x = "sample", by.y = "Sample")
 IND_meta_otu_tax <- merge(IND_meta_otu, IndTaxonomy2)
 str(IND_meta_otu_tax)
 
 #Also I want to get rid of the "(100)" in the genus
 #and I want to have them ordered according to abundance
 #meta_otu_RA_tax$Genus<- substr(meta_otu_RA_tax$genus, 1,nchar(meta_otu_RA_tax$genus)-5)
 #unique(meta_otu_RA_tax$Genus)
 #meta_otu_RA_tax$Genus<-factor(meta_otu_RA_tax$Genus, 
 #  levels = rev(c("Prevotella","Lactobacillus","Megasphaera","Ruminococcaceae_unclassified",
 #               "Veillonellaceae_unclassified","Streptococcus",
 #              "Lachnospiraceae_unclassified","Dialister", "Blautia",
 #             "Clostridiales_unclassified","Roseburia","Phascolarctobacterium","Alloprevotella",
 #            "Faecalibacterium","Clostridium_sensu_stricto", "Bacteria_unclassified")))
 
 #I will again give colors to our treatments:
 #Barxcol<-c(rep("#924e7d", 5), rep("#749b9d", 5), rep("#8b58c8",6), rep("#aeb646",6))
 divergingx_hcl(n=16, palette="Cividis") 
 
 IND16OTUcols<-c("#00214E", "#002B5F", "#0D376E", "#2E436D", "#434F6F", "#545C71", "#646875", "#747579", "#85837B", "#979179", "#A89F76", "#BAAD71",
                 "#CCBC6A", "#DECB61", "#F1DA54", "#FFE93F")
 
 sequential_hcl(n=16, palette="Oslo") 
 
 IND16OTUcols<-c("#FCFCFC", "#E5EAF4", "#CED7EC", "#B6C5E3", "#9EB4DB", "#86A2D3", "#6C92CB", "#4E81C3", "#3771B3", "#2F609A", "#275182", "#1F416A",
                 "#173353", "#10243E", "#0A1728", "#040404")
 
 sequential_hcl(n=16, palette="ag_Sunset")
 
 IND16OTUcols<-c("#4B1D91", "#6A1697", "#83109A", "#990F9C", "#AD169C","#BF219A", "#D02F95", "#DE3F8D", "#EB4F81", "#EF6674", "#F27B68", "#F38E60",
                 "#F2A15E", "#F0B265", "#EDC376", "#E7D39A")
 
 Testcol<-sample(IND16OTUcols, 11)
 
 IND16OTUcols<-divergingx_hcl(n=16, palette="Zissou 1") 
 
 IND11OTUcols<-divergingx_hcl(n=11, palette="Zissou 1") 
 
 IND11OTUcols<-divergingx_hcl(n=11, palette="Earth") 
 
 #IND11OTUcols<-divergingx_hcl(n=11, palette="Cividis") 
 
 unique(IND_meta_otu_tax$genus)
 
 IND_meta_otu_tax$Genus<-factor(IND_meta_otu_tax$genus, 
                                levels = rev(unique(IND_meta_otu_tax$genus)))
 
 #IND_meta_otu_tax$Genus<-
 
 IND_meta_otu_tax$Genus<-gsub("\\s*\\([^\\)]+\\)\\s*$","",as.character(IND_meta_otu_tax$Genus))
 
 unique(IND_meta_otu_tax$Genus)
 
 IND_meta_otu_tax$Genus<-factor(IND_meta_otu_tax$Genus, 
                                levels = rev(unique(IND_meta_otu_tax$Genus)))
 
 
 
 #The plot:
 ggplot(IND_meta_otu_tax, aes(x = sample, y = Abundance, fill = Genus)) + 
   geom_bar(stat = "identity",color="black",size=0.09) +
   scale_fill_manual(values =IND11OTUcols)+
   theme_bw ()+
   # Remove x axis title
   theme(axis.title.x = element_blank()) + 
   theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
   scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
   #ylim(c(0,1)) +
   guides(fill = guide_legend(reverse = F, keywidth = .9, keyheight = .9, ncol = 1)) +
   theme(axis.text.y=element_text(colour= "black",size=10))+
   theme(legend.text=element_text(size=8)) +
   #theme(legend.position="bottom") +
   #theme(axis.text.x = element_blank()) +
   ylab(paste0("Relative Abundance (top ", INDnum_genera, " OTUs)"))
 
 #I'll check that all the OTUs sum to 1
 rowSums(INDotu.summary)
 #they do
 

 
 
 # I'll try another approach because on the first samples
 
 Ind_OTU_abund002_mat <- Ind_OTUtableFover3F2_rel_o[,apply(Ind_OTUtableFover3F2_rel_o, 2, max)>0.02]
 
 dim(Ind_OTU_abund002_mat)
 
 Ind_OTU_abund005_mat <- Ind_OTUtableFover3F2_rel_o[,apply(Ind_OTUtableFover3F2_rel_o, 2, max)>0.05]
 
 dim(Ind_OTU_abund005_mat)
 
 Ind_OTU_abund005_mat.m<-melt(Ind_OTU_abund005_mat)
 
 #With ggplot2 the colors work better with NA's. 
 Ind_OTU_abund005_mat.m$valueII<-Ind_OTU_abund005_mat.m$value
 min(Ind_OTU_abund005_mat.m$valueII)
 Ind_OTU_abund005_mat.m$valueII[Ind_OTU_abund005_mat.m$valueII==9.536743e-08]<-NA
 
 Ind_OTU_abund005_mat.mTAX<-merge(Ind_OTU_abund005_mat.m, IndTaxonomy2, by.x="Var2",by.y="OTU")
 
 ggplot(Ind_OTU_abund005_mat.mTAX, aes(x = Var1, y = valueII, fill = family)) + 
   geom_bar(stat = "identity",color="black",size=0.09) +
   scale_fill_manual(values =IND11OTUcols)+
   theme_bw ()+
   # Remove x axis title
   theme(axis.title.x = element_blank()) + 
   theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
   scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
   #ylim(c(0,1)) +
   guides(fill = guide_legend(reverse = F, keywidth = .9, keyheight = .9, ncol = 1)) +
   theme(axis.text.y=element_text(colour= "black",size=10))+
   theme(legend.text=element_text(size=8)) +
   #theme(legend.position="bottom") +
   #theme(axis.text.x = element_blank()) +
   ylab("Relative Abundance (top 12 OTUs)")
 
# not better.
 
 # Another another approach:
 
 Ind_OTUtableFover3F2_rel_o.m<- melt(Ind_OTUtableFover3F2_rel_o)
 
 TAXbarplotdf<-melt(Ind_OTUtableFover3F2_rel_o.m)
 
 TAXbarplotdf<-merge(TAXbarplotdf, IndTaxonomy2, by.x="Var2",by.y="OTU")
 
 TAXbarplotdf$valueII<-TAXbarplotdf$value
 
 TAXbarplotdf$valueII[TAXbarplotdf$valueII==9.536743e-08]<-0
 
 TAXbarplot_ord<-dcast(TAXbarplotdf, Var1 ~ order, value.var = "valueII", sum)
 
 dim(TAXbarplot_ord)
 
 str(TAXbarplot_ord)
 
 rownames(TAXbarplot_ord)<-TAXbarplot_ord$Var1
 
 TAXbarplot_ord$Var1<-NULL
 
 TAXbarplot_ord_mat<-as.matrix(TAXbarplot_ord)
 
 rowSums(TAXbarplot_ord_mat)
 
 TAXbarplot_ord005_mat <- TAXbarplot_ord_mat[,apply(TAXbarplot_ord_mat, 2, max)>0.05]
 
 
 dim(TAXbarplot_ord005_mat)
 
 TAXbarplot_ord005_mat.m<-melt(TAXbarplot_ord005_mat)
 
 
 IND17Huepalvik2<-c("#002E60", "#00487A", "#00659E", "#3D81B1", "#739EC5", "#9FBBD8", "#C9D8E9",
                    "#2686A0","#6CAA9F","#EDEAC2", "#E3D4C3", "#CEB495", "#B69467", "#9C752E", "#805600",
                    "#5E3A00", "#3E2000")
 
 TAXbarplot_ord005_mat.m$Order<-TAXbarplot_ord005_mat.m$Var2
 
 TAXbarplot_ord005_mat.m$Order<-gsub("\\s*\\([^\\)]+\\)\\s*$","",as.character(TAXbarplot_ord005_mat.m$Order))
 
 
 IND16Huepalvik2<-c("#002E60", "#00487A", "#00659E", "#3D81B1", "#739EC5", "#9FBBD8", "#C9D8E9",
                    "#2686A0","#6CAA9F","#EDEAC2", "#CEB495", "#B69467", "#9C752E", "#805600",
                    "#5E3A00", "#3E2000")
 

 
 ggplot(TAXbarplot_ord005_mat.m, aes(x = Var1, y = value, fill = Order)) + 
   geom_bar(stat = "identity",color="black",size=0.09) +
   scale_fill_manual(values =rev(IND16Huepalvik2))+
   theme_bw ()+
   # Remove x axis title
   theme(axis.title.x = element_blank()) + 
   theme(axis.text.x = element_text(colour= IndHMxcol,angle = 90, hjust=0, vjust=0, size=10))+
   scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
   #ylim(c(0,1)) +
   guides(fill = guide_legend(reverse = F, keywidth = .9, keyheight = .9, ncol = 1)) +
   theme(axis.text.y=element_text(colour= "black",size=10))+
   theme(legend.text=element_text(size=8)) +
   #theme(legend.position="bottom") +
   #theme(axis.text.x = element_blank()) +
   ylab("Relative Abundance (top 16 orders)")
 
 # This will do!
 
 #OTU venns
 
 
 View(Ind_OTUtableFover3F2_rel_o)
# 
 Ind_OTUtableFover3F2_rel_o.m<- melt(Ind_OTUtableFover3F2_rel_o)
 
 Ind_OTUtableFover3F2_rel_o.m$S.type<-rep(1, nrow(Ind_OTUtableFover3F2_rel_o.m))  
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="1A"| Var1=="1B"| Var1=="1C" ]<-"Spring Water")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="2A"| Var1=="2B"| Var1=="2C" ]<-"Rural")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="3A"| Var1=="3B"| Var1=="3C" ]<-"Rural")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="4A"| Var1=="4B"| Var1=="4C" ]<-"City")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="5A"| Var1=="5B"| Var1=="5C" ]<-"City")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="6A"| Var1=="6B"| Var1=="6C" ]<-"City")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="7A"| Var1=="7B"| Var1=="7C" ]<-"Estuary")
 Ind_OTUtableFover3F2_rel_o.m<- within(Ind_OTUtableFover3F2_rel_o.m, S.type[Var1=="8A"| Var1=="8B"| Var1=="8C" ]<-"Estuary")
 
 #Lets convert the data frame back into matrix:
 
 Ind_OTUtableFover3F2_rel_o_mat_venn <- dcast(Ind_OTUtableFover3F2_rel_o.m, S.type+Var1~Var2, value.var="value")
 rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)<-interaction(Ind_OTUtableFover3F2_rel_o_mat_venn$S.type,Ind_OTUtableFover3F2_rel_o_mat_venn$Var1, sep="_")
 Ind_OTUtableFover3F2_rel_o_mat_venn<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[,3:ncol(Ind_OTUtableFover3F2_rel_o_mat_venn)])
 
 # Rural
 Ru_otus<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Rural",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]
                    [, which(colSums(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Rural",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]) != 6*(9.536743e-08))])
 str(Ru_otus)
 Ru_otusU<-colnames(Ru_otus)
 str(Ru_otusU)
 length(Ru_otusU)
 #15352
 
 City_otus<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("City",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]
                      [, which(colSums(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("City",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]) != 9*(9.536743e-08))])
 str(City_otus)
 City_otusU<-colnames(City_otus)
 str(City_otusU)
 length(City_otusU)
 #[1] 13932
 
 Est_otus<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Est",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]
                     [, which(colSums(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Est",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]) != 6*(9.536743e-08))])
 str(Est_otus)
 Est_otusU<-colnames(Est_otus)
 str(Est_otusU)
 length(Est_otusU)
 #[1] 11000
 
 SpWat_otus<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Spr",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]
                       [, which(colSums(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Spr",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),]) != 3*(9.536743e-08))])
 str(SpWat_otus)
 SpWat_otusU<-colnames(SpWat_otus)
 str(SpWat_otusU)
 length(SpWat_otusU)
 #[1] 3719
 
 input_otus<-list(Spring= SpWat_otusU, Rural=Ru_otusU, City=City_otusU, Estuary= Est_otusU)
 
 venn(input_otus)
 
 
 
 
 ###############################################################################
 
 
 # How does the abundance of ARHs and MGEs change from one sampling point to another
 
 
 ###############################################################################
 
 
 # I will try gamma GLMs. ARGs first. 
 
 View(Ind_results_relative_mat.m)
 
 # We need to have equal amount of observations per group we are comparing
 # City downstream has outliers
 
 Ind_results_relative_mat.mII<-Ind_results_relative_mat.m
 
 CityDownstream<-which(Ind_results_relative_mat.mII$Sample=="6A" |Ind_results_relative_mat.mII$Sample=="6B" | Ind_results_relative_mat.mII$Sample=="6C" )
 Ind_results_relative_mat.mII<-Ind_results_relative_mat.mII[-CityDownstream,]
 
 dim(Ind_results_relative_mat.mII)
 dim(Ind_results_relative_mat.m)
 
 View(Ind_results_relative_mat.mII)
 
 str(Ind_results_relative_mat.mII)
 
 
 #lets check this
 table(subset(Ind_results_relative_mat.mII$Var2,Ind_results_relative_mat.mII$value>9.536743e-07))
 #yes there are several ARGs that have detected in only one sample.
 
 library(fitdistrplus)
 library(logspline)
 
 ?fitdist
 
 #All values:
 fit_gammaARGs<-fitdist(Ind_results_relative_mat.mII$value,"gamma")
 plot(fit_gammaARGs)
 summary(fit_gammaARGs)
 
 #Fitting of the distribution ' gamma ' by maximum likelihood 
 #Parameters : 
 #  estimate  Std. Error
 #shape  0.1434189 0.002433366
 #rate  61.5272299 2.794877229
 #Loglikelihood:  32923.77   AIC:  -65843.54   BIC:  -65830.99 
 #Correlation matrix:
 #  shape      rate
 #shape 1.0000000 0.3735148
 #rate  0.3735148 1.0000000
 
 
 qqp(Ind_results_relative_mat.mII$value, "gamma", 
     shape = fit_gammaARGs$estimate[[1]], rate = fit_gammaARGs$estimate[[2]])
 #poor fit
 
 #With only detected ARGs:
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$value>9.536743e-07), "gamma", 
     shape = fit_gammaARGs$estimate[[1]], rate = fit_gammaARGs$estimate[[2]])
 
 # hmmmm
 
 colnames(Ind_results_relative_mat.mII)[1]<-"Sample_name"
 colnames(Ind_results_relative_mat.mII)[2]<-"Gene"
 
 
 #Let's try with some individual genes:
 
 unique(Ind_results_relative_mat.mII$Gene)
 #Note!!!!!!!!! #check at this point if you have "Gene" or "Var2" !!!!!
 
 fit_gamma_tetws<-fitdist(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="tetW"),"gamma")
 
 plot(fit_gamma_tetws)
 
 summary(fit_gamma_tetws)
 
 #Fitting of the distribution ' gamma ' by maximum likelihood 
 #Parameters : 
 # estimate  Std. Error
 #shape   0.2515632   0.0604716
 #rate  331.4264949 164.7512418
 #Loglikelihood:  155.6829   AIC:  -307.3657   BIC:  -305.2767 
 #Correlation matrix:
 # shape      rate
 #shape 1.0000000 0.4836201
 #rate  0.4836201 1.0000000
 
 
 View(Ind_results_relative_mat.mII)
 unique(Ind_results_relative_mat.mII$Gene)
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="tetW"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 #hmmm...
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="tnpA-04/IS6100"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="qacEdelta1_2"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="aadA_4"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="tet(36)_1"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 # This gene might cause problems
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="sul1_2"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="mexF"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="ISSm2"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="IncP_1alpha"), "gamma", 
     shape = fit_gamma_tetws$estimate[[1]], rate = fit_gamma_tetws$estimate[[2]])
 
 
 # with many genes the fit might be ok but with some not. 
 # Fit is poorer with those genes that are highly abundant in certain samples.
 
 # I wonder if I could use some other distribution...
 
 qqp(log10(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$value>9.536743e-07)), "norm")
 # this could be a possibility but I dont like the transformation
 
 qqp(log10(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="sul1_2")), "norm")
 # Does'nt look so nice
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="sul1_2"), "norm")
 # nope
 
 qqp(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="sul1_2"), "lnorm")
 # nope
 
 qqp(Ind_results_relative_mat.mII$value, "lnorm")
 # nope
 
 descdist(Ind_results_relative_mat.mII$value, discrete = FALSE)
 
 descdist(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$Gene=="sul1_2"), discrete = FALSE)
 
 descdist(subset(Ind_results_relative_mat.mII$value, Ind_results_relative_mat.mII$value>9.536743e-07), discrete = FALSE)
 
 #beta distribution..? I'll consider: https://cran.r-project.org/web/packages/betareg/index.html
 #https://stats.stackexchange.com/questions/169391/beta-distribution-glm-with-categorical-independents-and-proportional-response
 # https://rcompanion.org/handbook/J_02.html
 
 
 # Another alternative might be inverse gaussian. I'll come back to these later: 
 # the second answer here: https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
 #"If you are dealing with continuous non-negative outcome, then you could consider the Gamma distribution, or Inverse Gaussian distribution."
 #https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
 #https://stats.stackexchange.com/questions/238391/what-glm-family-for-continuous-positive-data
 
 
 # after doing some research it seems that my response variable should follow 
 # 1) beta distribution 2) gamma distribution 3) Inverse Gaussian. 
 
 # Beta distribution might be somewhat challenging since this "family" is not
 # in glm package.
 
 # I'll try again with gamma
 
 # from Crawley's book:
 #If the errors are gamma distributed, then the variance is proportional 
 #to the square of the mean (recall that with Poisson errors, the variance 
 #is equal to the mean). 
 

 
 var(Ind_results_relative_mat.mII$value)
 #[1] 0.0003272925
 mean(Ind_results_relative_mat.mII$value)
 #[1] 0.002332445
 2^(mean(Ind_results_relative_mat.mII$value))
 #[1] 1.001618
 
 library(plyr)
 
 Ind_mvtab <- ddply(Ind_results_relative_mat.mII,.(S.type),
                 summarise,value.mean = mean(value),value.var = var(value))
 
 plot(Ind_mvtab$value.mean,Ind_mvtab$value.var)
 
 # so I would say close enough! 
 # (Proportional = having the same or a constant ratio)
 # Estuary has a high variance but that is due to a biological reason
 
 # So I'll start fitting gamma model
 
 fit.I_ARGgamma1<-glm(value ~ S.type, data=Ind_results_relative_mat.mII, family=Gamma (link="log"))
 #Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
 # NA/NaN/Inf in 'x'
 # In addition: Warning message:
 # step size truncated due to divergence 
 
 # :D I know this problem. 
 # Let's see if I'm lucky and it's tet(36) again.
 
 Ind_results_relative_mat.mII_NOtet36<-subset(Ind_results_relative_mat.mII, Gene!="tet(36)_1")
 
 fit.I_ARGgamma1<-glm(value ~ S.type, data=Ind_results_relative_mat.mII_NOtet36, family=Gamma (link="log"))
 #Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
 #NA/NaN/Inf in 'x'
 # In addition: Warning message:
 #step size truncated due to divergence 
 
 # Ok so it wasn't tet(36), now I need tho find where the problem is.
 
 # First I'll try with some individual genes
 
 tnpAIS4glm<-subset(Ind_results_relative_mat.mII, Gene=="tnpA-02/IS4")
 
 fit.tnpAIS4glmgamma1<-glm(value ~ S.type, data=tnpAIS4glm, family=Gamma (link="log"))
 summary(fit.tnpAIS4glmgamma1)
 # Works
 
 tetWglm<-subset(Ind_results_relative_mat.mII, Gene=="tetW")
 
 tetWglmgamma1<-glm(value ~ S.type, data=tetWglm, family=Gamma (link="log"))
 summary(tetWglmgamma1)
 # Works
 
 table(subset(Ind_results_relative_mat.mII$Gene,Ind_results_relative_mat.mII$value>9.536743e-07))
 
 blaIMP_3glm<-subset(Ind_results_relative_mat.mII, Gene=="blaIMP_3")
 
 blaIMP_3glmgamma1<-glm(value ~ S.type, data=blaIMP_3glm, family=Gamma (link="log"))
 summary(blaIMP_3glmgamma1)
 # Works
 
 
 tet36glm<-subset(Ind_results_relative_mat.mII, Gene=="tet(36)_1")
 
 tet36glmgamma1<-glm(value ~ S.type, data=tet36glm, family=Gamma (link="log"))
 summary(tet36glmgamma1)
 # Works
 
 tetTglm<-subset(Ind_results_relative_mat.mII, Gene=="tetT")
 
 tetTglmgamma1<-glm(value ~ S.type, data=tetTglm, family=Gamma (link="log"))
 summary(tetTglmgamma1)
 # Works
 
 # This has happened before.... 
 
 # I'll set up the function.
 
 str(Ind_results_relative_mat.mII)
 
 unique(Ind_results_relative_mat.mII$Sample)
 
 # need to get rid of "Spring water" !!! 
 
 # If I'm lucky, that was causing me trouble. 
 
 Ind_results_relative_mat.mII_NOSW<-subset(Ind_results_relative_mat.mII, S.type!="Spring Water")
 
 Ind_results_relative_mat.mII_NOSW$S.type.f<-as.factor(Ind_results_relative_mat.mII_NOSW$S.type)
 
 Ind_results_relative_mat.mII_NOSW$S.type.f<-factor(Ind_results_relative_mat.mII_NOSW$S.type.f, levels=c("Rural","City","Estuary"))
 
 # I'll try this trick as I would like to control the order of sites
 Ind_results_relative_mat.mII_NOSW$S.type.f<-factor(Ind_results_relative_mat.mII_NOSW$S.type.f, levels=c("Rural","City","Estuary"))
 
 Ind_results_relative_mat.mII_NOSW<-droplevels(Ind_results_relative_mat.mII_NOSW)
 levels(Ind_results_relative_mat.mII_NOSW$S.type.f)
 str(Ind_results_relative_mat.mII_NOSW)
 
 library(multcomp)
 
 JohannaIsTheBest<-function(ARGx) {
   TEMP<-subset(Ind_results_relative_mat.mII_NOSW, Gene==ARGx)
   gammafit <-glm(value ~ S.type.f, data=TEMP, family=Gamma (link="log"),maxit = 10000)
   glht.gamma <- glht(gammafit, mcp(S.type.f = "Tukey"), p.adjust="none")#, p.adjust="BH")
   gamma_glht_df<-as.data.frame(summary(glht(glht.gamma))$test[-(1:2)])
   return(gamma_glht_df) 
 }
 
 JohannaIsTheBest("tetW")
 
 gamma_tmpARG <- unique(Ind_results_relative_mat.mII_NOSW$Gene)
 gamma_tmpARG_test <- lapply(gamma_tmpARG, JohannaIsTheBest)
 #There were 50 or more warnings (use warnings() to see the first 50)
 # IT WORKS!!!!! So Spring water was my problem!!
 head(gamma_tmpARG_test)
 names(gamma_tmpARG_test)
 names(gamma_tmpARG_test) <- gamma_tmpARG
 head(gamma_tmpARG_test)
 gamma_tmpARG_test.m<- melt(gamma_tmpARG_test)
 head(gamma_tmpARG_test.m)
 dim(gamma_tmpARG_test.m)
 #[1] 2244    4
 
 levels(gamma_tmpARG_test.m$variable)
 coefDF<-subset(gamma_tmpARG_test.m, variable=="coefficients")
 stderDF<-subset(gamma_tmpARG_test.m, variable=="sigma")
 tstatDF<-subset(gamma_tmpARG_test.m, variable=="tstat")
 pvalDF<-subset(gamma_tmpARG_test.m, variable=="pvalues")
 coefDF<-coefDF[,-c(1,2)]
 colnames(coefDF)<-c("coefficients","ARG_c")
 stderDF<-stderDF[,-c(1,2)]
 colnames(stderDF)<-c("std_error","ARG_s")
 tstatDF<-tstatDF[,-c(1,2)]
 colnames(tstatDF)<-c("tstat","ARG_t")
 pvalDF<-pvalDF[,-c(1,2)]
 colnames(pvalDF)<-c("p-value","ARG_p")
 
 Gamma_ARGs_results<-cbind(coefDF,stderDF,tstatDF,pvalDF)
 
 #I need the comparisons 
 df_ARGcomp_gamma<-as.data.frame(lapply(gamma_tmpARG, JohannaIsTheBest))
 getARGcomparisons<-rownames(df_ARGcomp_gamma)
 
 getARGcomparisons
 
 #Is everyhting ok?
 dim(Gamma_ARGs_results)
 length(unique(Ind_results_relative_mat.mII_NOSW$Gene))
 561 /187
 #[1] 3, Yes I had 3 comparisons
 
 #I will add a column that has the comparison
 Gamma_ARGs_results$comparison<-rep(getARGcomparisons,187)
 #And I want to remove columns that were used for checking the data
 Gamma_ARGs_results$Gene<-Gamma_ARGs_results$ARG_c
 Gamma_ARGs_results<-Gamma_ARGs_results[,-c(2,4,6,8)]
 #Renaming columns
 colnames(Gamma_ARGs_results)<-c("Delta.Estimate", "Std.error", "z-value", "p.value","Comparison", "Gene")
 #And finally the rownames
 rownames(Gamma_ARGs_results)<-interaction(Gamma_ARGs_results$Comparison,Gamma_ARGs_results$Gene, sep="_")
 head(rownames(Gamma_ARGs_results))
 #[1] "City - Rural_aac(6)-Ib_1"    "Estuary - Rural_aac(6)-Ib_1"
 #[3] "Estuary - City_aac(6)-Ib_1"  "City - Rural_aac(6)-Ib_2"   
 #[5] "Estuary - Rural_aac(6)-Ib_2" "Estuary - City_aac(6)-Ib_2" 
 
 Gamma_ARGs_results_sigP<-subset(Gamma_ARGs_results, p.value <= 0.05)
 #sig.PC_NC$p.adjFDR<-p.adjust(sig.PC_NC$p.value,method="BH")
 Gamma_ARGs_results_sigP$p.adj<-p.adjust(Gamma_ARGs_results_sigP$p.value, method="BH")
 
 #Gamma_ARGs_results_sigP$p.adj<-p.adjust(Gamma_ARGs_results_sigP$p.value, method="bonferroni")
 
 Gamma_ARGs_results_sigP.a<-subset(Gamma_ARGs_results_sigP, p.adj <= 0.05)
 
 dim(Gamma_ARGs_results_sigP.a)
 dim(Gamma_ARGs_results_sigP)
 
 unique(Gamma_ARGs_results_sigP.a$Gene)
 
 #[1] "aac(6)-Ib_1"     "aac(6)-Ib_2"     "aac(6)-Ib_3"     "aac(6)-II"      
 #[5] "aacA_aphD"       "aacC2"           "aacC4"           "aadA_2"         
 #[9] "aadA_3"          "aadA_5"          "aadA_6"          "aadA5_1"        
 #[13] "aadA5_2"         "aadD"            "aadE"            "acrA_1"         
 #[17] "acrA_2"          "acrA_4"          "acrA_5"          "acrB_1"         
 #[21] "acrF"            "acrR_1"          "acrR_2"          "acrR_3"         
 #[25] "ampC_1"          "ampC_2"          "ampC_5"          "aphA1_7"        
 #[29] "aphA3_1"         "aphA3_2"         "bacA_2"          "blaACT_2"       
 #[33] "blaACT_3"        "blaCMY_1"        "blaCMY_2"        "blaCMY_3"       
 #[37] "blaDHA"          "blaFOX"          "blaGES"          "blaIMP_3"       
 #[41] "blaKPC_1"        "blaMOX_blaCMY"   "blaOXA_1"        "blaOXA_2"       
 #[45] "blaOXA_3"        "blaOXA_4"        "blaPSE"          "blaSFO"         
 #[49] "blaSHV_1"        "blaSHV_2"        "blaTEM"          "blaTLA"         
 #[53] "blaVEB"          "catA1"           "catB3"           "catB8"          
 #[57] "cfxA"            "cmlA_2"          "cmlA_3"          "cmlA_4"         
 #[61] "cmxA"            "cphA_1"          "cphA_2"          "dfrA1_1"        
 #[65] "dfrA1_2"         "dfrA12"          "ereA"            "ereB"           
 #[69] "erm(35)"         "erm(36)"         "ermB"            "ermC"           
 #[73] "ermF"            "ermX"            "floR_1"          "gapA"           
 #[77] "IncN_rep"        "IncP_1alpha"     "IncP_1beta"      "IncP_oriT"      
 #[81] "IncQ_oriT"       "intI1_1"         "intI1_2"         "intI1_3"        
 #[85] "intI1_4"         "intI2_1"         "intI3_1"         "intI3_2"        
 #[89] "IS613"           "ISAba3"          "ISEfm1"          "ISPps"          
 #[93] "ISSm2"           "lmrA_1"          "lnuA_1"          "lnuB_1"         
 #[97] "lnuB_2"          "matA_mel"        "mdh"             "mdtA"           
 #[101] "mdtE"            "mdtF"            "mdtG_1"          "mdtH_1"         
 #[105] "mdtH_2"          "mdtL"            "mefA"            "merA"           
 #[109] "mexF"            "mphA_1"          "mphA_2"          "mtrD_3"         
 #[113] "nisB_2"          "oprD"            "oprJ"            "pbp"            
 #[117] "qacEdelta1_1"    "qacEdelta1_2"    "qacEdelta1_3"    "qacH_1"         
 #[121] "qacH_2"          "qnrA"            "qnrB"            "rarD_2"         
 #[125] "repA"            "rpoB"            "sat4"            "str"            
 #[129] "strA"            "sul1_2"          "sul2_1"          "sul2_2"         
 #[133] "sul3"            "tet(32)"         "tetA_2"          "tetA_B_1"       
 #[137] "tetA_B_2"        "tetA(P)"         "tetC_1"          "tetC_2"         
 #[141] "tetC_3"          "tetE"            "tetG_1"          "tetG_2"         
 #[145] "tetH"            "tetJ"            "tetL_2"          "tetM_1"         
 #[149] "tetM_2"          "tetM_3"          "tetO_1"          "tetO_2"         
 #[153] "tetPB_3"         "tetQ"            "tetR_2"          "tetR_3"         
 #[157] "tetS"            "tetW"            "tetX"            "tolC_1"         
 #[161] "tolC_2"          "tolC_3"          "Tp614"           "uidA"           
 #[165] "vanC_1"          "vanC_4"          "vanC2_vanC3"     "vanSC_2"        
 #[169] "vanTC_1"         "vanTC_3"         "vanYD_1"         "vatE_1"         
 #[173] "vatE_2"          "tnpA-05/IS26"    "tnpA-04/IS6100"  "tnpA-02/IS4"    
 #[177] "tnpA-01/IS21"    "tnpA-07/ISEcp1B" "tnpA-06/IS1216"  "tnpA-03/IS6"   
 
 # I can calculate the fold changes:
 Gamma_ARGs_results_sigP.a$Fold.Change<-exp(Gamma_ARGs_results_sigP.a$Delta.Estimate)
 
 
 # About the fold changes:
 
 JohannaIsTheBest("tetW")
 
 #                   coefficients     sigma      tstat    pvalues        type
 #City - Rural      -0.9174494 0.3711964  -2.471601 0.03579508 single-step
 #Estuary - Rural   -7.5962598 0.3711964 -20.464258 0.00000000 single-step
 #Estuary - City    -6.6788104 0.3711964 -17.992657 0.00000000 single-step
 #Warning messages:
 # 1: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
 #2: In chkdots(...) : Argument(s) ‘complete’ passed to ‘...’ are ignored
 
 # The output above shows "delta estimates" under the column coefficients. 
 # By these I mean for example that the estimate of  "Estuary" minus 
 #the estimate of "Rural" is -7.5962598 
 #we can check this:
 
 tetWglm<-subset(Ind_results_relative_mat.mII_NOSW, Gene=="tetW")
 
 tetWglmgamma1<-glm(value ~ S.type, data=tetWglm, family=Gamma (link="log"))
 
 summary(tetWglmgamma1)
 
 #Coefficients:
 #Estimate Std. Error t value Pr(>|t|)    
 #(Intercept)    -7.1841     0.2625 -27.371 3.22e-14 ***
 # S.typeEstuary  -6.6788     0.3712 -17.993 1.45e-11 ***
 #S.typeRural     0.9174     0.3712   2.472   0.0259 *  
 
 -6.6788  - 0.9174 
 #[1] 7.5962
 
 # Gamma distribution is a member of the exponential family and I used the log link. 
 # in other words, to backtransform the "delta estimate",
 # I need to do this:
 exp(-7.5962)
 #[1] 0.0005023568
 #This number is the fold change. BUT WE ARE NOT DONE!
 
 # to understand what is going on, I will do some calculations with the data: 
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="tetW" & S.type.f=="Estuary")))
 #[1] 9.536743e-07
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="tetW" & S.type.f=="Rural")))
 #[1] 0.001898514
 
 # in our glht output we had Estuary - Rural =  -7.5962. This was in log scale.
 # the corresponding in normal scale would be:
 9.536743e-07/0.001898514
 #[1] 0.0005023267
 # It's the same number that we got with exp(-7.5962)!
 # so gene tetW is ~ 0.0005023267 times less abundant in  Rural than in Estuary
 
 
 # With another gene:
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="aac(6)-II" & S.type.f=="City")))
 #[1] 0.0003538067
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="aac(6)-II" & S.type.f=="Rural")))
 #[1] 9.536743e-07
 0.0003538067/9.536743e-07
 #[1] 370.9932
 
 # and third:
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="sul3" & S.type.f=="City")))
 #[1] 0.0001178401
 mean(with (Ind_results_relative_mat.mII_NOSW, subset(value, Gene=="sul3" & S.type.f=="Rural")))
 #[1] 9.536743e-07
 0.0001178401/9.536743e-07
 #[1] 123.5643
 
 # Let's see some diagnostic plots at this point
 

 # The diagnostic plots are not beautiful but the fold changes seem to match. 
 # Boxplots are revealing, I might do them later.
 
 # Also the standard error needs to be transformed back
 Gamma_ARGs_results_sigP.a$expStd.error<-exp(Gamma_ARGs_results_sigP.a$Std.error)
 
 # log 2 transformation is handy for plotting
 
 Gamma_ARGs_results_sigP.a$log2Fold.Change<-log2(Gamma_ARGs_results_sigP.a$Fold.Change)
 
 Gamma_ARGs_results_sigP.a$log2expStd.e<-log2(Gamma_ARGs_results_sigP.a$expStd.error)
 
 Gamma_ARGs_results_sigP.a[["change"]] = ifelse(Gamma_ARGs_results_sigP.a[["log2Fold.Change"]] < 0, "decrease", "increase")
 
 # I want to plot the results  
 # I will need to get the genes that had significant differences in each 
 # comparison.
 
 str(Gamma_ARGs_results_sigP.a)
 
 Gamma_ARGs_results_sigP_RuCi<-subset(Gamma_ARGs_results_sigP.a, Comparison=="City - Rural")
 Gamma_ARGs_results_sigP_RuCi<-droplevels(Gamma_ARGs_results_sigP_RuCi)
 length(unique(Gamma_ARGs_results_sigP_RuCi$Gene))
 #[1] 104 
 # a lot of genes...
 #View(FC_Ru_Ci)
 
 Gamma_ARGs_results_sigP_CiEs<-subset(Gamma_ARGs_results_sigP.a, Comparison=="Estuary - City")
 Gamma_ARGs_results_sigP_CiEs<-droplevels(Gamma_ARGs_results_sigP_CiEs)
 length(unique(Gamma_ARGs_results_sigP_CiEs$Gene))
 #[1] 135
 #View(FC_Ci_Es)
 
 Gamma_ARGs_results_sigP_RuEs<-subset(Gamma_ARGs_results_sigP.a, Comparison=="Estuary - Rural")
 Gamma_ARGs_results_sigP_RuEs<-droplevels(Gamma_ARGs_results_sigP_RuEs)
 length(unique(Gamma_ARGs_results_sigP_RuEs$Gene))
 #[1] 152
 
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change >0))
 #[1] 29 11
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change <0))
 #[1] 123  11
 
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change >2))
 #[1] 27 11
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change < -2))
 #[1] 123  11
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change < -4))
 #[1] 111  11
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change < -6))
 #[1] 74 11
 log2(100)
 
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change < -6.643856))
 #[1] 61 11
 dim(subset(Gamma_ARGs_results_sigP_RuEs, log2Fold.Change > 6.643856))
 #[1] 117  11
 
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change < 0))
 #[1] 117  11
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change > 0))
 #[1] 18 11
 
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change < -6.643856))
 #[1] 71 11
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change > 6.643856))
 #[1]  1 11
 
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change > 0))
 #[1] 66 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change < 0))
 #[1] 38 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change < -6.643856))
 #[1]  6 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change > 6.643856))
 #[1]  8 11
 log2(50)
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change < -5.643856))
 #[1]  9 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change > 5.643856))
 #[1]  8 11
 log2(10)
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change < -3.321928))
 #[1]  30 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change > 3.321928))
 #[1] 25 11
 log2(5)
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change < -2.321928))
 #[1]  35 11
 dim(subset(Gamma_ARGs_results_sigP_RuCi, log2Fold.Change > 2.321928))
 #[1] 45 11
 
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change < -2.321928))
 #[1] 117  12
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change > 2.321928))
 #[1] 14 12
 
 
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change < -3.321928))
 #[1] 109  12
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change > 3.321928))
 #[1] 7 12
 
 log2(50)
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change < -5.643856))
 #[1] 84 12
 dim(subset(Gamma_ARGs_results_sigP_CiEs, log2Fold.Change > 5.643856))
 #[1]  1 12
 
 
 FCRuCi<-ggplot(subset(Gamma_ARGs_results_sigP_RuCi,log2Fold.Change< -2.321928|log2Fold.Change>2.321928), 
                aes(x = Gene, y = log2Fold.Change, fill = change)) + 
   geom_bar(stat="identity",width=.6) +
   scale_fill_manual(values = c("increase" = "darkorchid3", "decrease" = "chocolate4")) +
   theme_bw(base_family="Helvetica")+
   geom_errorbar(aes(ymin=log2Fold.Change-log2expStd.e, ymax=log2Fold.Change+log2expStd.e), width=.05,
                 position=position_dodge(.9)) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black", hjust=0.5, vjust=0, size=10))+
   theme(aspect.ratio = 2/1)+
   # theme_classic() + theme(axis.text.x=element_text(hjust=1,angle=90,size=8),
   #panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
   #panel.border=element_blank(), panel.grid.major.y=element_blank()) +
   ggtitle("Fold change from rural to city \n> 5 fold difference (p<0.05)")+
   #scale_y_continuous(limits=c(-7,16)) + 
   #labs(y="log2fold",x="")+
   #coord_flip(ylim=c(-7,6))
   #
   scale_y_continuous(name="log2fold",limits=c(-10,15), breaks=seq(-10,15,2))+
   coord_flip()
 
 FCRuCi
 
 FCCiEs<-ggplot(subset(Gamma_ARGs_results_sigP_CiEs,log2Fold.Change< -6.643856|log2Fold.Change>6.643856), 
                aes(x = Gene, y = log2Fold.Change, fill = change)) + 
   geom_bar(stat="identity",width=.6) +
   scale_fill_manual(values = c("increase" = "royalblue1", "decrease" = "mediumpurple3")) +
   theme_bw(base_family="Helvetica")+
   geom_errorbar(aes(ymin=log2Fold.Change-log2expStd.e, ymax=log2Fold.Change+log2expStd.e), width=.05,
                 position=position_dodge(.9)) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black", hjust=0.5, vjust=0, size=10))+
   theme(aspect.ratio = 2/1)+
   # theme_classic() + theme(axis.text.x=element_text(hjust=1,angle=90,size=8),
   #panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
   #panel.border=element_blank(), panel.grid.major.y=element_blank()) +
   ggtitle("Fold change from city to estuary \n> 100 fold difference (p<0.05)")+
   #scale_y_continuous(limits=c(-7,16)) + 
   #labs(y="log2fold",x="")+
   #coord_flip(ylim=c(-7,6))
   #
   scale_y_continuous(name="log2fold",limits=c(-14,8), breaks=seq(-14,8,2))+
   coord_flip()
 FCCiEs
 
 FCRuEs<-ggplot(subset(Gamma_ARGs_results_sigP_RuEs,log2Fold.Change< -6.643856|log2Fold.Change>6.643856), 
                aes(x = Gene, y = log2Fold.Change, fill = change)) + 
   geom_bar(stat="identity",width=.6) +
   scale_fill_manual(values = c("increase" = "royalblue1", "decrease" = "chocolate4")) +
   theme_bw(base_family="Helvetica")+
   geom_errorbar(aes(ymin=log2Fold.Change-log2expStd.e, ymax=log2Fold.Change+log2expStd.e), width=.05,
                 position=position_dodge(.9)) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black", hjust=0.5, vjust=0, size=10))+
   theme(aspect.ratio = 2/1)+
   # theme_classic() + theme(axis.text.x=element_text(hjust=1,angle=90,size=8),
   #panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
   #panel.border=element_blank(), panel.grid.major.y=element_blank()) +
   ggtitle("Fold change from rural to estuary \n> 100 fold difference (p<0.05)")+
   #scale_y_continuous(limits=c(-7,16)) + 
   #labs(y="log2fold",x="")+
   #coord_flip(ylim=c(-7,6))
   #
   scale_y_continuous(name="log2fold",limits=c(-12,14), breaks=seq(-12,14,2))+
   coord_flip()
 FCRuEs
 
 grid.arrange(FCRuCi,FCCiEs,FCRuEs, nrow=1,ncol=3 )
 
 str(Gamma_ARGs_results_sigP.a)
 
 Gamma_ARGs_results_sigP.a$helpcol1<-rep(1, nrow(Gamma_ARGs_results_sigP.a))  
 
 unique(Gamma_ARGs_results_sigP.a$Comparison)
 
 #[1] "City - Rural"    "Estuary - Rural" "Estuary - City" 
 
 Gamma_ARGs_results_sigP.a<- within(Gamma_ARGs_results_sigP.a, helpcol1[Comparison=="City - Rural"]<-"from rural to city") 
 Gamma_ARGs_results_sigP.a<- within(Gamma_ARGs_results_sigP.a, helpcol1[Comparison=="Estuary - Rural"]<-"from rural to estuary") 
 Gamma_ARGs_results_sigP.a<- within(Gamma_ARGs_results_sigP.a, helpcol1[Comparison=="Estuary - City" ]<-"from city to estuary") 
 
 Gamma_ARGs_results_sigP.a$Change<-interaction(Gamma_ARGs_results_sigP.a$change, Gamma_ARGs_results_sigP.a$helpcol1, sep=" ")
 
 write.table(Gamma_ARGs_results_sigP.a, file="Gamma_ARGs_results_sigPa.txt",sep="\t")
 
 #################################
 
 
 #### NMDS plots again  ####
 
 # This is their place in the paper #
 
##################################
 
 
 ###############################################################################
 
 #             Do OTUs explain ARGs? Mantel's tests
 
 ###############################################################################
 
 View(Ind_OTUtableFover3F2_rel_o)
 
 Ind_OTUtableFover3F2_rel_o_PRO<-Ind_OTUtableFover3F2_rel_o[4:nrow(Ind_OTUtableFover3F2_rel_o),]
 
 Ind_ARGandMGE_dist<-vegdist(Ind_results_relative_matII)
 Ind_OTU_dist<-vegdist(Ind_OTUtableFover3F2_rel_o_PRO)
 
 length(Ind_ARGandMGE_dist)
 length(Ind_OTU_dist)
 
 #lets see whart the relationship looks like:
 plot(Ind_ARGandMGE_dist, Ind_OTU_dist)
 #a Mantel test with kendall (I think Kendall is the right choice here)
 mantel(Ind_ARGandMGE_dist, Ind_OTU_dist, method = "kendall")
 
 #Mantel statistic based on Kendall's rank correlation tau 
 
 #Call:
 #mantel(xdis = Ind_ARGandMGE_dist, ydis = Ind_OTU_dist, method = "kendall") 
 
 #Mantel statistic r: 0.2147 
 #Significance: 0.007 
 
 #Upper quantiles of permutations (null model):
 #  90%    95%  97.5%    99% 
 #0.0868 0.1212 0.1582 0.2012 
 #Permutation: free
 #Number of permutations: 999
 
 #spearman:
 mantel( Ind_OTU_dist, Ind_ARGandMGE_dist, method = "spearman")
 
 
 #Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = Ind_OTU_dist, ydis = Ind_ARGandMGE_dist, method = "spearman") 

#Mantel statistic r: 0.2982 
 #     Significance: 0.009 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.139 0.187 0.226 0.283 
#Permutation: free
#Number of permutations: 999

 
 # Just out of curiosity, what would happen if we would not remove spring water
 # that had no ARGs.
 
 View(Ind_results_relative_mat)
 View(Ind_OTUtableFover3F2_rel_o)
 
 TESTInd_ARGandMGE_dist<-vegdist(Ind_results_relative_mat)
 TESTInd_OTU_dist<-vegdist(Ind_OTUtableFover3F2_rel_o)
 
 mantel(TESTInd_OTU_dist, TESTInd_ARGandMGE_dist, method = "spearman")
 
 #Mantel statistic r: 0.6085 
 #Significance: 0.001 
 
 # but this is incorrect.
 

 #Pearson:
 
 mantel(Ind_ARGandMGE_dist, Ind_OTU_dist)
 
 #Mantel statistic based on Pearson's product-moment correlation 
 
 #Call:
 #mantel(xdis = Ind_ARGandMGE_dist, ydis = Ind_OTU_dist) 
 
 #Mantel statistic r: 0.3765 
 #Significance: 0.001 
 
 #Upper quantiles of permutations (null model):
 #90%   95% 97.5%   99% 
 #0.106 0.150 0.184 0.233 
 #Permutation: free
 #Number of permutations: 999
 
 #mantel tests for groups according to sample type
 
 View(Ind_results_relative_matII)
 View(Ind_OTUtableFover3F2_rel_o_PRO)
 View(Ind_results_relative_mat_venn)
 View(Ind_OTUtableFover3F2_rel_o_mat_venn)
 
 # Rural
 Ru_ARGs<-as.matrix(Ind_results_relative_mat_venn[grep("Rural",rownames(Ind_results_relative_mat_venn)),])
 Ru_OTUs<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Rural",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),])
 
 RU_ARGandMGE_dist<-vegdist(Ru_ARGs)
 RU_OTU_dist<-vegdist(Ru_OTUs)
 
 mantel(RU_ARGandMGE_dist, RU_OTU_dist, method = "spearman")
 
 #Mantel statistic based on Spearman's rank correlation rho 
 
 #Call:
 #mantel(xdis = RU_ARGandMGE_dist, ydis = RU_OTU_dist, method = "spearman") 
 
 #Mantel statistic r: 0.6536 
 #     Significance: 0.019444 
 
 #Upper quantiles of permutations (null model):
 # 90%   95% 97.5%   99% 
 #0.319 0.486 0.582 0.713 
 #Permutation: free
 #Number of permutations: 719
 
 # Urban
 City_ARGs<-as.matrix(Ind_results_relative_mat_venn[grep("City",rownames(Ind_results_relative_mat_venn)),])
 City_OTUs<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("City",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),])
 
 CI_ARGandMGE_dist<-vegdist(City_ARGs)
 CI_OTU_dist<-vegdist(City_OTUs)
 
 mantel(CI_ARGandMGE_dist, CI_OTU_dist, method = "spearman")
 
 #Mantel statistic based on Spearman's rank correlation rho 
 
 #Call:
 #mantel(xdis = CI_ARGandMGE_dist, ydis = CI_OTU_dist, method = "spearman") 
 
 #Mantel statistic r: 0.2216 
 #     Significance: 0.113 
 
 #Upper quantiles of permutations (null model):
 # 90%   95% 97.5%   99% 
 #0.238 0.328 0.374 0.442 
 #Permutation: free
 #Number of permutations: 999

 
 # Estuary
 Est_ARGs<-as.matrix(Ind_results_relative_mat_venn[grep("Estuary",rownames(Ind_results_relative_mat_venn)),])
 Est_OTUs<-as.matrix(Ind_OTUtableFover3F2_rel_o_mat_venn[grep("Estuary",rownames(Ind_OTUtableFover3F2_rel_o_mat_venn)),])
 
 EST_ARGandMGE_dist<-vegdist(Est_ARGs)
 EST_OTU_dist<-vegdist(Est_OTUs)
 
 mantel(EST_ARGandMGE_dist, EST_OTU_dist, method = "spearman")
 
 #Mantel statistic based on Spearman's rank correlation rho 
 
 #Call:
 #mantel(xdis = EST_ARGandMGE_dist, ydis = EST_OTU_dist, method = "spearman") 
 
 #Mantel statistic r: 0.175 
 #     Significance: 0.20972 
 
 #Upper quantiles of permutations (null model):
 # 90%   95% 97.5%   99% 
 #0.414 0.518 0.596 0.754 
 #Permutation: free
 #Number of permutations: 719

 ####################
 
 # Sum abundance
 
 ######################
 
 
 library(ggpubr)
 
 View(Ind_results_relative_mat)
 
 Ind_results_relative_matNAs<-Ind_results_relative_mat
 
 min(Ind_results_relative_matNAs)
 
 Ind_results_relative_matNAs[Ind_results_relative_matNAs == 9.536743e-07] <- NA
 
 View(I_metadata)
 
 I_metadata$ARGsumAbund<-rowSums(Ind_results_relative_matNAs, na.rm=TRUE)
 
 I_metadata$ARGsumAbund[ I_metadata$ARGsumAbund==0]<-NA
 
 

 I_metadata$S.type<-rep(1, nrow(I_metadata))  
 I_metadata<- within(I_metadata, S.type[sample=="1A"| sample=="1B"| sample=="1C" ]<-"Spring Water")
 I_metadata<- within(I_metadata, S.type[sample=="2A"| sample=="2B"| sample=="2C" ]<-"Rural")
 I_metadata<- within(I_metadata, S.type[sample=="3A"| sample=="3B"| sample=="3C" ]<-"Rural")
 I_metadata<- within(I_metadata, S.type[sample=="4A"| sample=="4B"| sample=="4C" ]<-"City")
 I_metadata<- within(I_metadata, S.type[sample=="5A"| sample=="5B"| sample=="5C" ]<-"City")
 I_metadata<- within(I_metadata, S.type[sample=="6A"| sample=="6B"| sample=="6C" ]<-"City")
 I_metadata<- within(I_metadata, S.type[sample=="7A"| sample=="7B"| sample=="7C" ]<-"Estuary")
 I_metadata<- within(I_metadata, S.type[sample=="8A"| sample=="8B"| sample=="8C" ]<-"Estuary")
 
 levels(I_metadata$S.type)
 
 I_metadata$S.type.f<-factor(I_metadata$S.type, levels=c("Spring Water","Rural","City","Estuary"))
 
 ARGSumAbund_comparisons <- list(c("Rural", "City"), c("City","Estuary"), c("Rural","Estuary"))
 
 
 ARGSumAbundBp<-ggplot(subset(I_metadata,S.type.f !="Spring Water"), aes(x=S.type.f, y=ARGsumAbund, colour=S.type.f)) +
    #
    #geom_point(aes(fill=factor(S.type.f)), size=3, shape=21, colour="grey20",alpha=0.7,
     #          position=position_jitter(width=0.01, height=0.01)) +
    #
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    scale_color_manual(values=c( "chocolate4","purple4","royalblue"))+
    #facet_grid(C_code ~ .) +
    theme_bw(base_size = 14)+
    theme(strip.background = element_rect(colour = "black", fill = "white"))+
    #theme(strip.background = element_blank())+
    theme(panel.border = element_rect(colour = "black"))+
   # theme(axis.text.x=element_blank())+
    #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
    #scale_y_log10()+
    #scale_y_log2(limits = c(0.0001,10))+
    #theme(strip.text.y = element_text(size=10, angle = 360))+
    #
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    #guides(color=guide_legend(title="Sample type"))+
    theme(legend.position="none")+
    geom_hline(yintercept = mean(na.omit(I_metadata$ARGsumAbund)), linetype = 2)+
    #stat_compare_means(method = "kruskal.test", label.y = 2.5, label.x=2)+ 
    stat_compare_means(comparisons = ARGSumAbund_comparisons,method = "wilcox.test",label = "p.signif")
 
 ARGSumAbundBp
 
 
#######################################################################

# Alpha diversity analysis

######################################################################

unique(I_metadata$color)

View(I_metadata)

#I_metadata$S.type<-rep(1, nrow(I_metadata))  
#I_metadata<- within(I_metadata, S.type[sample=="1A"| sample=="1B"| sample=="1C" ]<-"Spring Water")
#I_metadata<- within(I_metadata, S.type[sample=="2A"| sample=="2B"| sample=="2C" ]<-"Rural")
#I_metadata<- within(I_metadata, S.type[sample=="3A"| sample=="3B"| sample=="3C" ]<-"Rural")
#I_metadata<- within(I_metadata, S.type[sample=="4A"| sample=="4B"| sample=="4C" ]<-"City")
#I_metadata<- within(I_metadata, S.type[sample=="5A"| sample=="5B"| sample=="5C" ]<-"City")
#I_metadata<- within(I_metadata, S.type[sample=="6A"| sample=="6B"| sample=="6C" ]<-"City")
#I_metadata<- within(I_metadata, S.type[sample=="7A"| sample=="7B"| sample=="7C" ]<-"Estuary")
#I_metadata<- within(I_metadata, S.type[sample=="8A"| sample=="8B"| sample=="8C" ]<-"Estuary")


anovaShan <- aov(I_metadata$OTU_divH~ I_metadata$name.f)
anovaShan

#Call:
# aov(formula = I_metadata$OTU_divH ~ I_metadata$name.f)

#Terms:
# I_metadata$name.f Residuals
#Sum of Squares          18.605523  0.756767
#Deg. of Freedom                 7        16

#Residual standard error: 0.2174809
#Estimated effects may be unbalanced


TukeyHSD(anovaShan)

#Fit: aov(formula = I_metadata$OTU_divH ~ I_metadata$name.f)

#$`I_metadata$name.f`
#diff         lwr        upr     p adj
#Cattle Farm-Spring Water                   0.82479152  0.21000899  1.4395740 0.0050919
#Chicken Slaughterhouse-Spring Water       -1.98857835 -2.60336088 -1.3737958 0.0000001
#Hospital-Spring Water                     -1.82975794 -2.44454047 -1.2149754 0.0000004
#City-Spring Water                         -0.92266725 -1.53744978 -0.3078847 0.0017535
#City Downstream-Spring Water              -1.01655738 -1.63133991 -0.4017749 0.0006441
#Estuary-Freshwater-Spring Water           -1.22865670 -1.84343923 -0.6138742 0.0000748
#Estuary-Seawater-Spring Water             -1.41270297 -2.02748550 -0.7979204 0.0000133
#Chicken Slaughterhouse-Cattle Farm        -2.81336987 -3.42815240 -2.1985873 0.0000000
#Hospital-Cattle Farm                      -2.65454946 -3.26933198 -2.0397669 0.0000000
#City-Cattle Farm                          -1.74745877 -2.36224130 -1.1326762 0.0000008
#City Downstream-Cattle Farm               -1.84134890 -2.45613143 -1.2265664 0.0000004
#Estuary-Freshwater-Cattle Farm            -2.05344822 -2.66823075 -1.4386657 0.0000001
#Estuary-Seawater-Cattle Farm              -2.23749449 -2.85227702 -1.6227120 0.0000000
#Hospital-Chicken Slaughterhouse            0.15882041 -0.45596211  0.7736029 0.9823792
#City-Chicken Slaughterhouse                1.06591110  0.45112857  1.6806936 0.0003847
#City Downstream-Chicken Slaughterhouse     0.97202097  0.35723844  1.5868035 0.0010325
#Estuary-Freshwater-Chicken Slaughterhouse  0.75992165  0.14513912  1.3747042 0.0103820
#Estuary-Seawater-Chicken Slaughterhouse    0.57587538 -0.03890715  1.1906579 0.0749636
#City-Hospital                              0.90709069  0.29230816  1.5218732 0.0020751
#City Downstream-Hospital                   0.81320055  0.19841802  1.4279831 0.0057824
#Estuary-Freshwater-Hospital                0.60110123 -0.01368130  1.2158838 0.0577235
#Estuary-Seawater-Hospital                  0.41705496 -0.19772756  1.0318375 0.3272000
#City Downstream-City                      -0.09389013 -0.70867266  0.5208924 0.9992712
#Estuary-Freshwater-City                   -0.30598945 -0.92077198  0.3087931 0.6740387
#Estuary-Seawater-City                     -0.49003572 -1.10481825  0.1247468 0.1740548
#Estuary-Freshwater-City Downstream        -0.21209932 -0.82688185  0.4026832 0.9222283
#Estuary-Seawater-City Downstream          -0.39614559 -1.01092812  0.2186369 0.3843620
#Estuary-Seawater-Estuary-Freshwater       -0.18404627 -0.79882880  0.4307363 0.9612254

# So lots of significant differences. 


#OTUcomplist<-split(t(combn(levels(I_metadata$name.f), 2)), seq(nrow(t(combn(levels(I_metadata$name.f), 2)))))


# This is the way to do it.

#library(ggpubr)

compare_means(OTU_divH ~ name.f,  data = I_metadata,
              ref.group = ".all.", method = "t.test")
# From here: 
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

#A tibble: 8 × 8
#.y.      group1 group2                             p       p.adj p.format p.signif method
#<chr>    <chr>  <chr>                          <dbl>       <dbl> <chr>    <chr>    <chr> 
#1 OTU_divH .all.  Spring Water           0.000837      0.005       0.00084  ***      T-test
#2 OTU_divH .all.  Cattle Farm            0.00000000168 0.000000013 1.7e-09  ****     T-test
#3 OTU_divH .all.  Chicken Slaughterhouse 0.0507        0.2         0.05071  ns       T-test
#4 OTU_divH .all.  Hospital               0.000287      0.002       0.00029  ***      T-test
#5 OTU_divH .all.  City                   0.872         1           0.87222  ns       T-test
#6 OTU_divH .all.  City Downstream        0.751         1           0.75115  ns       T-test
#7 OTU_divH .all.  Estuary-Freshwater     0.161         0.48        0.16081  ns       T-test
#8 OTU_divH .all.  Estuary-Seawater       0.0290        0.15        0.02901  *        T-test


OTUdivH<-ggplot(I_metadata, aes(x=name.f, y=OTU_divH)) +
  #
  geom_point(aes(fill=factor(name.f)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  scale_fill_manual(values=c( "springgreen4","sienna3",
          "chocolate4","mediumorchid1", "mediumpurple3",
            "purple4","royalblue1","blue4"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
   theme(legend.position="none")+
  #guides(fill=guide_legend(title="Site"))+
  geom_hline(yintercept = mean(I_metadata$OTU_divH), linetype = 2)+
 # stat_compare_means(method = "anova", label.y = 9, label.x=2)+ 
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE) 

OTUdivH
  
  # stat_compare_means(method = "t.test",comparisons =OTUcomplist)+
  #stat_compare_means(method = "anova",label.y = 9)
#geom_signif(comparisons = split(t(combn(levels(I_metadata$name.f), 2)), seq(nrow(t(combn(levels(I_metadata$name.f), 2))))), 
            # test = "t.test",map_signif_level = FALSE)
#theme(legend.title=element_blank())


# for the fun of it, I'll do this for ARGs too

I_metadata$ARG_divH<-diversity(Ind_results_relative_mat)

I_metadataARDdiv<-subset(I_metadata, name.f!="Spring Water")

I_metadataARDdiv<-droplevels(I_metadataARDdiv)

compare_means(ARG_divH ~ name.f,  data = I_metadataARDdiv,
               ref.group = ".all.", method = "t.test")


# A tibble: 7 × 8
#.y.      group1 group2                        p  p.adj p.format p.signif method
#<chr>    <chr>  <chr>                     <dbl>  <dbl> <chr>    <chr>    <chr> 
#1 ARG_divH .all.  Cattle Farm            0.00283  0.017  0.00283  **       T-test
#2 ARG_divH .all.  Chicken Slaughterhouse 0.132    0.26   0.13218  ns       T-test
#3 ARG_divH .all.  Hospital               0.0544   0.19   0.05437  ns       T-test
#4 ARG_divH .all.  City                   0.0163   0.081  0.01626  *        T-test
#5 ARG_divH .all.  City Downstream        0.195    0.26   0.19509  ns       T-test
#6 ARG_divH .all.  Estuary-Freshwater     0.000242 0.0017 0.00024  ***      T-test
#7 ARG_divH .all.  Estuary-Seawater       0.0485   0.19   0.04851  *        T-test

I_metadataARGdivplot<-I_metadata

I_metadataARGdivplot<- within(I_metadataARGdivplot, ARG_divH[sample=="1A"| sample=="1B"| sample=="1C" ]<-NA)

I_metadata<-I_metadataARGdivplot

rm(I_metadataARGdivplot)

ARGdivH<-ggplot(I_metadata, aes(x=name.f, y=ARG_divH)) +
  #
  geom_point(aes(fill=factor(name.f)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  #
  scale_fill_manual(values=c("springgreen4","sienna3",
                              "chocolate4","mediumorchid1", "mediumpurple3",
                              "purple4","royalblue1","blue4"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(fill=guide_legend(title="Site"))+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(I_metadata$ARG_divH)), linetype = 2)+
 # stat_compare_means(method = "anova", label.y = 5, label.x=2)+ 
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE) 

ARGdivH

plot(I_metadata$OTU_divH, I_metadata$ARG_divH)

library(ggrepel)

DivScatter<-ggplot(I_metadataARDdiv, aes(x= OTU_divH, y= ARG_divH))+
  geom_point(aes(fill=factor(name.f)), size=3, shape=21, colour="grey20",alpha=0.7,
              position=position_jitter(width=0.01, height=0.01)) +
  geom_text_repel(min.segment.length = 0,aes(label=name.f))+
theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  scale_fill_manual(values=c("sienna3",
                             "chocolate4","mediumorchid1", "mediumpurple3",
                             "purple4","royalblue1","blue4"))+
  #theme(axis.title.x=element_blank(),
       # axis.title.y=element_blank())+
  theme(legend.position="none")+
  labs(x ="OTU diversity", y = "ARG diversity")

DivScatter
  #guides(fill=guide_legend(title="Site"))

ggarrange(ARGSumAbundBp, ARGdivH, OTUdivH,ncol = 3,nrow = 1)
          
#DivScatter, ncol = 2, nrow = 2,heights = c(3, 3,), widths = c(4,1))

# these plots are again finalized with Inkscape


# Let's see ARGs and MGEs separately (Marko asked)

View(Ind_results_relative_mat.m)
str(Ind_results_relative_mat.m)


ARGassayAnnot<-Ind_results[,c(1,5,6)]


ARGassayAnnot<-unique(ARGassayAnnot)


names(Ind_results_relative_mat.m)[2]<-"Gene"

Ind_results_relative_mat.m_annot<-merge(Ind_results_relative_mat.m,ARGassayAnnot,by="Gene")

names(Ind_results_relative_mat.m_annot)[7]<-"Classification"

which(Ind_results_relative_mat.m_annot$Classification=="Integron")
Ind_results_relative_mat.m_annot$Classification[2017:2184]<-"MGE"

levels(Ind_results_relative_mat.m_annot$Classification)

Ind_results_relative_mat.m_annot<-droplevels(Ind_results_relative_mat.m_annot)

View(I_metadata)

Ind_results_relative_mat.m_annotOnlyARGs<- subset(Ind_results_relative_mat.m_annot, Classification != "MGE" )

Ind_results_relative_mat.m_annotOnlyARGs<-droplevels(Ind_results_relative_mat.m_annotOnlyARGs)

Ind_result_matOnlyARGs <- dcast(Ind_results_relative_mat.m_annotOnlyARGs, Sample~Gene, value.var="value")

View(Ind_result_matOnlyARGs)

rownames(Ind_result_matOnlyARGs)<-Ind_result_matOnlyARGs$Sample
Ind_result_matOnlyARGs$Sample<-NULL

Ind_result_matOnlyARGs<-as.matrix(Ind_result_matOnlyARGs)

I_metadata$ONLY_ARG_divH<-diversity(Ind_result_matOnlyARGs)


I_metadata<- within(I_metadata, ONLY_ARG_divH[sample=="1A"| sample=="1B"| sample=="1C" ]<-NA)



ONLY_ARG_divHBp<-ggplot(I_metadata, aes(x=name.f, y=ONLY_ARG_divH)) +
  #
  geom_point(aes(fill=factor(name.f)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  #
  scale_fill_manual(values=c("springgreen4","sienna3",
                             "chocolate4","mediumorchid1", "mediumpurple3",
                             "purple4","royalblue1","blue4"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(fill=guide_legend(title="Site"))+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(I_metadata$ONLY_ARG_divH)), linetype = 2)+
  # stat_compare_means(method = "anova", label.y = 5, label.x=2)+ 
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE) 

ONLY_ARG_divHBp


Ind_results_relative_mat.m_annotOnlyMGEs<- subset(Ind_results_relative_mat.m_annot, Classification == "MGE" )

Ind_results_relative_mat.m_annotOnlyMGEs<-droplevels(Ind_results_relative_mat.m_annotOnlyMGEs)

Ind_result_matOnlyMGEs <- dcast(Ind_results_relative_mat.m_annotOnlyMGEs, Sample~Gene, value.var="value")

View(Ind_result_matOnlyMGEs)

rownames(Ind_result_matOnlyMGEs)<-Ind_result_matOnlyMGEs$Sample
Ind_result_matOnlyMGEs$Sample<-NULL

Ind_result_matOnlyMGEs<-as.matrix(Ind_result_matOnlyMGEs)

I_metadata$ONLY_MGE_divH<-diversity(Ind_result_matOnlyMGEs)


I_metadata<- within(I_metadata, ONLY_MGE_divH[sample=="1A"| sample=="1B"| sample=="1C" ]<-NA)



ONLY_MGE_divHBp<-ggplot(I_metadata, aes(x=name.f, y=ONLY_MGE_divH)) +
  #
  geom_point(aes(fill=factor(name.f)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  #
  scale_fill_manual(values=c("springgreen4","sienna3",
                             "chocolate4","mediumorchid1", "mediumpurple3",
                             "purple4","royalblue1","blue4"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(fill=guide_legend(title="Site"))+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(I_metadata$ONLY_MGE_divH)), linetype = 2)+
  # stat_compare_means(method = "anova", label.y = 5, label.x=2)+ 
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE) 



####  Sum abundance for ARGs and MGEs separately


Ind_result_matOnlyARGsNA<-Ind_result_matOnlyARGs

min(Ind_result_matOnlyARGsNA)

Ind_result_matOnlyARGsNA[Ind_result_matOnlyARGsNA == 9.536743e-07] <- NA

View(I_metadata)

I_metadata$ONLY_ARGsumAbund<-rowSums(Ind_result_matOnlyARGsNA, na.rm=TRUE)

I_metadata$ONLY_ARGsumAbund[I_metadata$ONLY_ARGsumAbund==0]<-NA




ONLY_ARGSumAbundBp<-ggplot(subset(I_metadata,S.type.f !="Spring Water"), aes(x=S.type.f, y=ONLY_ARGsumAbund, colour=S.type.f)) +
  #
  #geom_point(aes(fill=factor(S.type.f)), size=3, shape=21, colour="grey20",alpha=0.7,
  #          position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_color_manual(values=c( "chocolate4","purple4","royalblue"))+
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  # theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  #
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(color=guide_legend(title="Sample type"))+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(I_metadata$ONLY_ARGsumAbund)), linetype = 2)+
  #stat_compare_means(method = "kruskal.test", label.y = 2.5, label.x=2)+ 
  stat_compare_means(comparisons = ARGSumAbund_comparisons,method = "wilcox.test",label = "p.signif")



Ind_result_matOnlyMGEsNAs<-Ind_result_matOnlyMGEs

min(Ind_result_matOnlyMGEsNAs)

Ind_result_matOnlyMGEsNAs[Ind_result_matOnlyMGEsNAs == 9.536743e-07] <- NA

View(I_metadata)

I_metadata$ONLY_MGEsumAbund<-rowSums(Ind_result_matOnlyMGEsNAs, na.rm=TRUE)

I_metadata$ONLY_MGEsumAbund[I_metadata$ONLY_MGEsumAbund==0]<-NA




ONLY_MGESumAbundBp<-ggplot(subset(I_metadata,S.type.f !="Spring Water"), aes(x=S.type.f, y=ONLY_MGEsumAbund, colour=S.type.f)) +
  #
  #geom_point(aes(fill=factor(S.type.f)), size=3, shape=21, colour="grey20",alpha=0.7,
  #          position=position_jitter(width=0.01, height=0.01)) +
  #
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_color_manual(values=c( "chocolate4","purple4","royalblue"))+
  #facet_grid(C_code ~ .) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  #theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  # theme(axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  #scale_y_log10()+
  #scale_y_log2(limits = c(0.0001,10))+
  #theme(strip.text.y = element_text(size=10, angle = 360))+
  #
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(color=guide_legend(title="Sample type"))+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(I_metadata$ONLY_MGEsumAbund)), linetype = 2)+
  #stat_compare_means(method = "kruskal.test", label.y = 2.5, label.x=2)+ 
  stat_compare_means(comparisons = ARGSumAbund_comparisons,method = "wilcox.test",label = "p.signif")


ggarrange(ARGSumAbundBp, ARGdivH, OTUdivH,ncol = 3,nrow = 1)


ggarrange(ONLY_ARGSumAbundBp, ONLY_MGESumAbundBp,ONLY_ARG_divHBp, ONLY_MGE_divHBp, ncol = 2,nrow = 2)

ggarrange(ONLY_ARGSumAbundBp, ONLY_ARG_divHBp, ONLY_MGESumAbundBp,ONLY_MGE_divHBp, ARGSumAbundBp, ARGdivH, 
          labels = c("Only ARGs SA", "Only ARGs H", "Only MGEs SA", "Only MGEs H", "ARGs&MGEs SA", "ARGs&MGEs H"), ncol = 2,nrow = 3)



##################################################################################

#  Can we say something about the river health according to indicator bacteria?

#################################################################################

View(TAXbarplotdf)

#Enterobacterales
grep("Enterobacterales",TAXbarplotdf$order)

grep("Enterobacteri",TAXbarplotdf$order)

TAXbarplotdf[6049,] 

grep("Enterobacteriales",TAXbarplotdf$order)

#old order name in Silva....

grep("Enterobacteriaceae",TAXbarplotdf$family)

grep("Streptococcaceae",TAXbarplotdf$family)

grep("Staphylococcaceae",TAXbarplotdf$family)

grep("Campylobacteraceae",TAXbarplotdf$family)

grep("Moraxellaceae",TAXbarplotdf$family)

grep("Enterococcaceae",TAXbarplotdf$family)

grep("Pseudomonadales",TAXbarplotdf$order)
#Moraxellaceae belong to these

grep("Clostridiaceae",TAXbarplotdf$family)

#Ok, there are a lot of these...

TAXbarplotdf$Family<-TAXbarplotdf$family

TAXbarplotdf$Family<-gsub("\\s*\\([^\\)]+\\)\\s*$","",as.character(TAXbarplotdf$Family))

head(TAXbarplotdf$Family)

Family_include_list<-c(" Enterobacteriaceae"," Streptococcaceae"," Staphylococcaceae"," Campylobacteraceae"," Moraxellaceae",
" Enterococcaceae"," Pseudomonadaceae"," Clostridiaceae", " Aeromonadaceae")

EskapeEtAl<-subset(TAXbarplotdf, TAXbarplotdf$Family %in% Family_include_list)

dim(TAXbarplotdf)

dim(EskapeEtAl)

unique(EskapeEtAl$genus)

# I wonder if I could get rid of unclassified and uncultured stuff

unc_rows<-grep("unc",EskapeEtAl$genus)

EskapeEtAl_II<-EskapeEtAl[-unc_rows,]

EskapeEtAl_II<-droplevels(EskapeEtAl_II)

unique(EskapeEtAl_II$genus)

dim(EskapeEtAl_II)

EskapeEtAl_II$valueII[EskapeEtAl_II$valueII==0]<-NA

EskapeEtAl_II$Genus<-gsub("\\s*\\([^\\)]+\\)\\s*$","",as.character(EskapeEtAl_II$genus))



 ggplot(EskapeEtAl_II,aes(x=Var1, y=Genus, fill=valueII))+
   geom_tile()+
   scale_fill_continuous_sequential(palette =  "Oslo",
     #scale_fill_continuous_sequential(palette =  "Dark Mint",
      begin = 0.8, end = 0, na.value="white", name="Relative abundance")+
   #scale_fill_continuous_sequential(palette =  "Lajolla",
   #I don't want the most darkest color and I want to reverse the colors
   #begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
   theme(panel.background=element_rect(fill="black", colour="black")) +
   theme(axis.text.y=element_text(colour= "black",size=8))+
   theme(axis.text.x = element_text(colour= "black",angle = 90, hjust=0, vjust=0, size=10))+
   theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=0, size=10))+
   theme(panel.border=element_blank())+
   theme(axis.title.x = element_blank()) + 
   theme(axis.title.y = element_blank()) + 
   theme(legend.position="bottom")+
   theme(legend.key.size=unit(0.2, "cm"))+
   theme(legend.key.width=unit(0.5, "cm"))
 
 unique(EskapeEtAl_II$Genus)
 
 Pathogens_and_traffickers_list<-c(" Tolumonas"," Thiopseudomonas", " Staphylococcus"," Campylobacter",
                                  " Kosakonia", " Serratia", " Enterococcus"," Klebsiella", " Acinetobacter",
                                  " Plesiomonas", " Providencia", " Samsonia"," Arcobacter"," Moraxella",
                                  " Pseudomonas", " Dickeya"," Aeromonas",  " Streptococcus", " Pectobacterium",
                                  " Escherichia-Shigella")
 
 EskapeEtAl.F<-subset(EskapeEtAl_II, EskapeEtAl_II$Genus %in%  Pathogens_and_traffickers_list)
 
 
 ggplot(EskapeEtAl.F,aes(x=Var1, y=Genus, fill=valueII))+
    geom_tile()+
    scale_fill_continuous_sequential(palette =  "Oslo",
                                     #scale_fill_continuous_sequential(palette =  "Dark Mint",
                                     begin = 0.8, end = 0, na.value="white", name="Relative abundance")+
    #scale_fill_continuous_sequential(palette =  "Lajolla",
    #I don't want the most darkest color and I want to reverse the colors
    #begin = 1, end = 0.1, na.value="white", name="Relative abundance")+
    theme(panel.background=element_rect(fill="black", colour="black")) +
    theme(axis.text.y=element_text(colour= "black",size=8))+
    theme(axis.text.x = element_text(colour= "black",angle = 90, hjust=0, vjust=0, size=10))+
    theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=0, size=10))+
    theme(panel.border=element_blank())+
    theme(axis.title.x = element_blank()) + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.position="bottom")+
    theme(legend.key.size=unit(0.2, "cm"))+
    theme(legend.key.width=unit(0.5, "cm"))
 
 
 EskapeEtAl.F[is.na(EskapeEtAl.F)]<-0
 
 
 EskapeEtAl.F_mat<-dcast(EskapeEtAl.F, Var1 ~ Genus, value.var = "valueII", sum)
 
rownames( EskapeEtAl.F_mat)<- EskapeEtAl.F_mat$Var1

EskapeEtAl.F_mat.m<-melt(EskapeEtAl.F_mat) 


EskapeEtAl.F_mat.m$value[EskapeEtAl.F_mat.m$value==0]<-NA

EskapeEtAl.F_mat.m<-merge(EskapeEtAl.F_mat.m, I_metadata, by.x="Var1", by.y="sample")


ggplot(EskapeEtAl.F_mat.m, aes(x=name.f, y=value, colour=name.f)) + 
   geom_boxplot()+
   geom_jitter() +
   scale_color_manual(values=c("springgreen4","sienna3",
                               "chocolate4","mediumorchid1", "mediumpurple3",
                               "purple4","royalblue1","blue4"))+
   theme_bw(base_size = 14)+
   theme(strip.background = element_rect(colour = "black", fill = "white"))+
   #theme(strip.background = element_blank())+
   theme(panel.border = element_rect(colour = "black"))+
   theme(axis.title.x=element_blank(),
         axis.title.y=element_blank())+
  geom_hline(yintercept = mean(EskapeEtAl.F_mat.m$value, na.rm = TRUE), linetype = 2)+
   stat_compare_means(label = "p.signif", method = "wilcox.test",
                      ref.group = ".all.", hide.ns = TRUE, na.rm = TRUE)+
   facet_wrap(variable ~ ., scales='free')


unique(EskapeEtAl.F_mat.m$variable)

EskapeEtAl.F_mat.m.s<-subset(EskapeEtAl.F_mat.m, variable != " Campylobacter" & 
                                variable != " Moraxella" & variable != " Providencia" & variable !=" Thiopseudomonas")


EskapeEtAl.F_mat.m.s<-droplevels(EskapeEtAl.F_mat.m.s)

unique(EskapeEtAl.F_mat.m.s$variable)

compare_means(value ~ S.type.f,  data = EskapeEtAl.F_mat.m, ref.group = "Spring Water",
              method = "wilcox.test", na.rm = TRUE)

ggplot(EskapeEtAl.F_mat.m.s, aes(x=S.type.f, y=value, colour=S.type.f)) + 
   geom_boxplot()+
   geom_jitter() +
   scale_color_manual(values=c("springgreen4", "chocolate4",
                               "purple4","royalblue"))+
   theme_bw(base_size = 14)+
   theme(strip.background = element_rect(colour = "black", fill = "white"))+
   #theme(strip.background = element_blank())+
   theme(panel.border = element_rect(colour = "black"))+
   theme(axis.title.x=element_blank(),
         axis.title.y=element_blank())+
   theme(axis.text.x=element_blank())+
   facet_wrap(variable ~ ., scales='free')+
   #geom_hline(yintercept = mean(EskapeEtAl.F_mat.m.s$value, na.rm = TRUE), linetype = 2)+
   #takes this for all the data.
   stat_compare_means(label = "p.signif", method = "wilcox.test",
                      ref.group = "Spring Water", hide.ns = TRUE, na.rm = TRUE,label.y = 0.00065)

fancy_scientific <- function(l) {
   # turn in to character string in scientific notation
   l <- format(l, scientific = TRUE)
   # quote the part before the exponent to keep all the digits
   l <- gsub("^(.*)e", "'\\1'e", l)
   # turn the 'e+' into plotmath format
   l <- gsub("e", "%*%10^", l)
   # return this as an expression
   parse(text=l)
}
  


ggplot(EskapeEtAl.F_mat.m.s, aes(x=S.type.f, y=value, colour=S.type.f)) + 
   geom_boxplot()+
   geom_jitter() +
   scale_color_manual(values=c("springgreen4", "chocolate4",
                               "purple4","royalblue"))+
   theme_bw(base_size = 14)+
   theme(strip.background = element_rect(colour = "black", fill = "white"))+
   #theme(strip.background = element_blank())+
   theme(panel.border = element_rect(colour = "black"))+
   theme(axis.title.x=element_blank(),
         axis.title.y=element_blank())+
   theme(axis.text.x=element_blank())+
   facet_wrap(variable ~ ., scales="free_y")+
  scale_y_continuous(expand = c(0.1,0),labels=fancy_scientific)+
   #facet_wrap(variable ~ ., scales='free')+
   #geom_hline(yintercept = mean(EskapeEtAl.F_mat.m.s$value, na.rm = TRUE), linetype = 2)+
   #takes this for all the data.
   stat_compare_means(label = "p.signif", method = "wilcox.test",
                      ref.group = "Spring Water", hide.ns = TRUE, na.rm = TRUE)
  

################################################################

# Supporting information: Mobility potential, Correlation between ARGs and MGEs

################################################################

# I need ARGs and MGEs in separate matrices by sampling area

str(Ru_ARGs) # I have these matrices, but it might be easier to take these again from a data frame.
str(Ind_assays)

View(Ind_assays)
View(Ind_results) # this dataframe has my corrections to tnpA's 

str(Ind_results_relative_mat.mII_NOSW) # as well as this

ARGassayAnnot<-Ind_results[,c(1,5,6)]

ARGassayAnnot<-unique(ARGassayAnnot)

Ind_results_relative_mat.mII_NOSW_annot<-merge(Ind_results_relative_mat.mII_NOSW,ARGassayAnnot,by="Gene")

names(Ind_results_relative_mat.mII_NOSW_annot)[8]<-"Classification"

which(Ind_results_relative_mat.mII_NOSW_annot$Classification=="Integron")
Ind_results_relative_mat.mII_NOSW_annot$Classification[1513:1638]<-"MGE"

Ru_ARGsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="Rural" & Classification != "MGE" )
Ru_ARGsCorr<-droplevels(Ru_ARGsCorr)

str(Ru_ARGsCorr)

levels(Ru_ARGsCorr$Classification)

Ru_ARGsCorr$Sam_Stype<-interaction(Ru_ARGsCorr$Sample, Ru_ARGsCorr$S.type, sep= "_") 

Ru_ARGsCorr_mat<-dcast(Ru_ARGsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(Ru_ARGsCorr_mat)<-Ru_ARGsCorr_mat$Sam_Stype

Ru_ARGsCorr_mat[1]<-NULL

Ru_ARGsCorr_mat<-as.matrix(Ru_ARGsCorr_mat)

dim(Ru_ARGsCorr_mat)

Ru_ARGsCorr_mat<-Ru_ARGsCorr_mat[,which(colSums(Ru_ARGsCorr_mat)>6*9.536743e-07)]

min(Ru_ARGsCorr_mat[Ru_ARGsCorr_mat!=min(Ru_ARGsCorr_mat)])
# 5.399767e-05
Ru_ARGsCorr_mat<-Ru_ARGsCorr_mat[,which(colSums(Ru_ARGsCorr_mat)>2*5.399767e-05)]

Ru_MGEsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="Rural" & Classification == "MGE")

Ru_MGEsCorr<-droplevels(Ru_MGEsCorr)

str(Ru_MGEsCorr)

levels(Ru_MGEsCorr$Gene)

Ru_MGEsCorr$Sam_Stype<-interaction(Ru_MGEsCorr$Sample, Ru_MGEsCorr$S.type, sep= "_") 

Ru_MGEsCorr_mat<-dcast(Ru_MGEsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(Ru_MGEsCorr_mat)<-Ru_MGEsCorr_mat$Sam_Stype

Ru_MGEsCorr_mat[1]<-NULL

Ru_MGEsCorr_mat<-as.matrix(Ru_MGEsCorr_mat)

dim(Ru_MGEsCorr_mat)

Ru_MGEsCorr_mat<-Ru_MGEsCorr_mat[,which(colSums(Ru_MGEsCorr_mat)>6*9.536743e-07)]

min(Ru_MGEsCorr_mat[Ru_MGEsCorr_mat!=min(Ru_MGEsCorr_mat)])
# 5.99961e-05
Ru_MGEsCorr_mat<-Ru_MGEsCorr_mat[,which(colSums(Ru_MGEsCorr_mat)>2*5.99961e-05)]

Ru_ARGsCorr_mat.df<-as.data.frame(Ru_ARGsCorr_mat)

Ru_MGEsCorr_mat.df<-as.data.frame(Ru_MGEsCorr_mat)

library(tidyverse)
library(corrr)
library(corrplot)

# I'll do as I have done and find a way to visualize.

##For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctRu_ARG_MGEfdr<-corr.test(Ru_ARGsCorr_mat, Ru_MGEsCorr_mat, use="pairwise", method="spearman",
                           adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 
str(ctRu_ARG_MGEfdr)
#I will put the results in a data frame
ctRu_ARG_MGEfdr_ci<-ctRu_ARG_MGEfdr$ci

#corr.test combines the names and shortens them as it returns the results 
# in a "melted format". I will get the Assay names through the function "cor":

corRu_ARG_MGE<-cor(Ru_ARGsCorr_mat, Ru_MGEsCorr_mat,use="pairwise.complete.obs", method="spearman")
dim(corRu_ARG_MGE)
corRu_ARG_MGE.m<-melt(corRu_ARG_MGE)
View(corRu_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corRu_ARG_MGE.m)
dim(ctRu_ARG_MGEfdr_ci)
#And they have the same length. 
#I will add r and p columns to the df that has assay names:
corRu_ARG_MGE.m$r<-ctRu_ARG_MGEfdr_ci$r
corRu_ARG_MGE.m$p.adj<-ctRu_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corRu_ARG_MGE.m)
#They do so I can remove the column named value
corRu_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corRu_ARG_MGE.m_sig<-filter(corRu_ARG_MGE.m, p.adj<= 0.05)
corRu_ARG_MGE.m_sig<-filter(corRu_ARG_MGE.m_sig, r>= 0.8)

corRu_ARG_MGE.m_sig_mat<-dcast(corRu_ARG_MGE.m_sig, Var2 ~ Var1, value.var = "r")

rownames(corRu_ARG_MGE.m_sig_mat)<-corRu_ARG_MGE.m_sig_mat$Var2

corRu_ARG_MGE.m_sig_mat$Var2<-NULL

corRu_ARG_MGE.m_sig_mat<-as.matrix(corRu_ARG_MGE.m_sig_mat)

heatcol3<- hcl.colors(10, "YlOrRd", rev = TRUE)

makocolcorr<-hcl.colors(10, "Mako", rev=TRUE)


corrplot(corRu_ARG_MGE.m_sig_mat, na.label = "square",na.label.col = "white", is.corr = FALSE,col.lim = c(0.8, 1), 
         col = heatcol3, cl.pos = 'b')

corrplot(corRu_ARG_MGE.m_sig_mat, method = 'color',na.label = "square",na.label.col = "white", 
         is.corr = FALSE,col.lim = c(0.8, 1),addgrid.col="grey", 
         col = makocolcorr, cl.pos = 'b',tl.col = "black", tl.cex = 0.8)

# Same for City

City_ARGsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="City" & Classification != "MGE")
City_ARGsCorr<-droplevels(City_ARGsCorr)

City_ARGsCorr$Sam_Stype<-interaction(City_ARGsCorr$Sample, City_ARGsCorr$S.type, sep= "_") 

City_ARGsCorr_mat<-dcast(City_ARGsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(City_ARGsCorr_mat)<-City_ARGsCorr_mat$Sam_Stype

City_ARGsCorr_mat[1]<-NULL

City_ARGsCorr_mat<-as.matrix(City_ARGsCorr_mat)

dim(City_ARGsCorr_mat)

City_ARGsCorr_mat<-City_ARGsCorr_mat[,which(colSums(City_ARGsCorr_mat)>6*9.536743e-07)]

min(City_ARGsCorr_mat[City_ARGsCorr_mat!=min(City_ARGsCorr_mat)])
# 3.035682e-05
City_ARGsCorr_mat<-City_ARGsCorr_mat[,which(colSums(City_ARGsCorr_mat)>2*3.035682e-05)]

City_MGEsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="City" & Classification == "MGE")

City_MGEsCorr<-droplevels(City_MGEsCorr)

str(City_MGEsCorr)

levels(City_MGEsCorr$Gene)

City_MGEsCorr$Sam_Stype<-interaction(City_MGEsCorr$Sample, City_MGEsCorr$S.type, sep= "_") 

City_MGEsCorr_mat<-dcast(City_MGEsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(City_MGEsCorr_mat)<-City_MGEsCorr_mat$Sam_Stype

City_MGEsCorr_mat[1]<-NULL

City_MGEsCorr_mat<-as.matrix(City_MGEsCorr_mat)

dim(City_MGEsCorr_mat)

City_MGEsCorr_mat<-City_MGEsCorr_mat[,which(colSums(City_MGEsCorr_mat)>6*9.536743e-07)]

min(City_MGEsCorr_mat[City_MGEsCorr_mat!=min(City_MGEsCorr_mat)])
# [1] 3.210817e-05
City_MGEsCorr_mat<-City_MGEsCorr_mat[,which(colSums(City_MGEsCorr_mat)>2*3.210817e-05)]

##For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctCity_ARG_MGEfdr<-corr.test(City_ARGsCorr_mat, City_MGEsCorr_mat, use="pairwise", method="spearman",
                             adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 
str(ctCity_ARG_MGEfdr)
#I will put the results in a data frame
ctCity_ARG_MGEfdr_ci<-ctCity_ARG_MGEfdr$ci

#corr.test combines the names and shortens them as it returns the results 
# in a "melted format". I will get the Assay names through the function "cor":

corCity_ARG_MGE<-cor(City_ARGsCorr_mat, City_MGEsCorr_mat,use="pairwise.complete.obs", method="spearman")
dim(corCity_ARG_MGE)
corCity_ARG_MGE.m<-melt(corCity_ARG_MGE)
View(corCity_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corCity_ARG_MGE.m)
dim(ctCity_ARG_MGEfdr_ci)
#And they have the same length. 
#I will add r and p columns to the df that has assay names:
corCity_ARG_MGE.m$r<-ctCity_ARG_MGEfdr_ci$r
corCity_ARG_MGE.m$p.adj<-ctCity_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corCity_ARG_MGE.m)
#They do so I can remove the column named value
corCity_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corCity_ARG_MGE.m_sig<-filter(corCity_ARG_MGE.m, p.adj<= 0.05)
corCity_ARG_MGE.m_sig<-filter(corCity_ARG_MGE.m_sig, r>= 0.8)

corCity_ARG_MGE.m_sig_mat<-dcast(corCity_ARG_MGE.m_sig, Var2 ~ Var1, value.var = "r")

rownames(corCity_ARG_MGE.m_sig_mat)<-corCity_ARG_MGE.m_sig_mat$Var2

corCity_ARG_MGE.m_sig_mat$Var2<-NULL

corCity_ARG_MGE.m_sig_mat<-as.matrix(corCity_ARG_MGE.m_sig_mat)


corrplot(corCity_ARG_MGE.m_sig_mat, na.label = "square",na.label.col = "white", is.corr = FALSE,col.lim = c(0.8, 1), 
         col = heatcol3, cl.pos = 'b', tl.cex = 0.5)

corrplot(corCity_ARG_MGE.m_sig_mat, method = 'color',na.label = "square",na.label.col = "white", 
         is.corr = FALSE,col.lim = c(0.8, 1),addgrid.col="grey", 
         col = makocolcorr, cl.pos = 'b',tl.col = "black", tl.cex = 0.5)


# same for Estuary


Est_ARGsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="Estuary" & Classification != "MGE")
Est_ARGsCorr<-droplevels(Est_ARGsCorr)

Est_ARGsCorr$Sam_Stype<-interaction(Est_ARGsCorr$Sample, Est_ARGsCorr$S.type, sep= "_") 

Est_ARGsCorr_mat<-dcast(Est_ARGsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(Est_ARGsCorr_mat)<-Est_ARGsCorr_mat$Sam_Stype

Est_ARGsCorr_mat[1]<-NULL

Est_ARGsCorr_mat<-as.matrix(Est_ARGsCorr_mat)

dim(Est_ARGsCorr_mat)

Est_ARGsCorr_mat<-Est_ARGsCorr_mat[,which(colSums(Est_ARGsCorr_mat)>6*9.536743e-07)]

min(Est_ARGsCorr_mat[Est_ARGsCorr_mat!=min(Est_ARGsCorr_mat)])
# 0.0001091332
Est_ARGsCorr_mat<-Est_ARGsCorr_mat[,which(colSums(Est_ARGsCorr_mat)>2*0.0001091332)]

Est_MGEsCorr<- subset(Ind_results_relative_mat.mII_NOSW_annot, S.type.f=="Estuary" & Classification == "MGE")

Est_MGEsCorr<-droplevels(Est_MGEsCorr)


Est_MGEsCorr$Sam_Stype<-interaction(Est_MGEsCorr$Sample, Est_MGEsCorr$S.type, sep= "_") 

Est_MGEsCorr_mat<-dcast(Est_MGEsCorr,Sam_Stype ~ Gene, value.var = "value")

rownames(Est_MGEsCorr_mat)<-Est_MGEsCorr_mat$Sam_Stype

Est_MGEsCorr_mat[1]<-NULL

Est_MGEsCorr_mat<-as.matrix(Est_MGEsCorr_mat)

dim(Est_MGEsCorr_mat)

Est_MGEsCorr_mat<-Est_MGEsCorr_mat[,which(colSums(Est_MGEsCorr_mat)>6*9.536743e-07)]

min(Est_MGEsCorr_mat[Est_MGEsCorr_mat!=min(Est_MGEsCorr_mat)])
#  0.0001192026
Est_MGEsCorr_mat<-Est_MGEsCorr_mat[,which(colSums(Est_MGEsCorr_mat)>2* 0.0001192026)]



##For getting correlation coefficients and false discovery rate adjusted p-values I will use
#function corr.test from package psych:
ctEst_ARG_MGEfdr<-corr.test(Est_ARGsCorr_mat, Est_MGEsCorr_mat, use="pairwise", method="spearman",
                            adjust="fdr",alpha=.05,ci=TRUE)
#Spearmans correlation is used because we might have 
#nonlinear dependencies. 
str(ctEst_ARG_MGEfdr)
#I will put the results in a data frame
ctEst_ARG_MGEfdr_ci<-ctEst_ARG_MGEfdr$ci

#corr.test combines the names and shortens them as it returns the results 
# in a "melted format". I will get the Assay names through the function "cor":

corEst_ARG_MGE<-cor(Est_ARGsCorr_mat, Est_MGEsCorr_mat,use="pairwise.complete.obs", method="spearman")
dim(corEst_ARG_MGE)
corEst_ARG_MGE.m<-melt(corEst_ARG_MGE)
View(corEst_ARG_MGE.m)
#Now I have ARG, MGE and their correlation in the data frame. But this is missing the 
#p-values, which I can get from corr.test. 
dim(corEst_ARG_MGE.m)
dim(ctEst_ARG_MGEfdr_ci)
#And they have the same length. 
#I will add r and p columns to the df that has assay names:
corEst_ARG_MGE.m$r<-ctEst_ARG_MGEfdr_ci$r
corEst_ARG_MGE.m$p.adj<-ctEst_ARG_MGEfdr_ci$p
#I will compare the r's and values to make sure they match
View(corEst_ARG_MGE.m)
#They do so I can remove the column named value
corEst_ARG_MGE.m$value<-NULL
#Now I want to have only positive correletions >0.8 and p.adj <0.05:
corEst_ARG_MGE.m_sig<-filter(corEst_ARG_MGE.m, p.adj<= 0.05)
corEst_ARG_MGE.m_sig<-filter(corEst_ARG_MGE.m_sig, r>= 0.8)

corEst_ARG_MGE.m_sig_mat<-dcast(corEst_ARG_MGE.m_sig, Var2 ~ Var1, value.var = "r")

rownames(corEst_ARG_MGE.m_sig_mat)<-corEst_ARG_MGE.m_sig_mat$Var2

corEst_ARG_MGE.m_sig_mat$Var2<-NULL

corEst_ARG_MGE.m_sig_mat<-as.matrix(corEst_ARG_MGE.m_sig_mat)


corrplot(corEst_ARG_MGE.m_sig_mat, method = 'color',na.label = "square",na.label.col = "white", 
         is.corr = FALSE,col.lim = c(0.8, 1), addgrid.col="grey", 
         col = makocolcorr, cl.pos = 'b',tl.col = "black", tl.cex = 1)



# These plots are then merged and made nice using Inkscape





### contact me at johanna.muurinen@onehealth.fi if you have questions ####



###########         The end           #######################

