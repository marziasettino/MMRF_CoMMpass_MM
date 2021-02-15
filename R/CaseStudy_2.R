library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(stringr)
library(ggplot2)
library(MMRFBiolinks)
library(EnsDb.Hsapiens.v86)



#HMGB3 LMNB1 KIF22 TSPEAR 



# IMPORTANT 1! If you want to reproduce exaclty the case study 2 by using complete MMRF Researcher Gateway (RG) clinical data, 
# you should gain access to MMRF Researcher Gateway (RG) for downloading the following files:
# i)MMRF_CoMMpass_IA14a_All_Canonical_Variants.csv
#ii)MMRF_CoMMpass_IA14a_All_Canonical_NS_Variants_ENSG_Mutation_Counts.csv
#iii)MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv



# Once this is done, you should import them into your R environmnet for feeding the following Workflow.
#Note1: pay attention to import MMRF-RG dataset into R environment including the heading of columns.
variant.ann<- MMRF_CoMMpass_IA14a_All_Canonical_Variants
trt<- MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP




names(variant.ann)[52]<-"Gene_name"
trt<-trt[trt$bestresp!="",]

#--------------HMGB3----------------------------

variant.ann1<-variant.ann[(variant.ann$Gene_name=="HMGB3"), ]


summary.var1<-MMRFRG_VariantCountPlot(variant.ann1,trt,topN=50,filenm=NULL,height=20, width=40)



#--------------LMNB1----------------------------No Variant----------

variant.ann2<-variant.ann[variant.ann$Gene_name=="LMNB1", ]


summary.var2<-MMRFRG_VariantCountPlot(variant.ann2,trt,topN=50,filenm=NULL,height=20, width=40)



#--------------KIF22----------------------------

variant.ann3<-variant.ann[(variant.ann$Gene_name=="KIF22") , ]


summary.var3<-MMRFRG_VariantCountPlot(variant.ann3,trt,topN=50,filenm=NULL,height=20, width=40)




#--------------TSPEAR----------------------------

variant.ann4<-variant.ann[variant.ann$Gene_name=="TSPEAR", ]

summary.var4<-MMRFRG_VariantCountPlot(variant.ann4,trt,topN=50,filenm=NULL,height=20, width=40)


#-------------ALL------------------

variant.ann<-variant.ann[(variant.ann$Gene_name=="HMGB3") |
                           (variant.ann$Gene_name=="LMNB1") |
                           (variant.ann$Gene_name=="KIF22") |
                           (variant.ann$Gene_name=="TSPEAR"), ]






summary.var<-MMRFRG_VariantCountPlot(variant.ann,trt,topN=50,filenm=NULL,height=20, width=40)

#complete.response<- summary.var[summary.var$bestresp=="Complete Response",]





