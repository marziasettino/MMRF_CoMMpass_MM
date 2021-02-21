library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(stringr)
library(ggplot2)
library(MMRFBiolinks)
library(EnsDb.Hsapiens.v86)



#KCNK2 ACTL8 EPHA5 SAGE1



# IMPORTANT 1! If you want to reproduce exaclty the case study 2 by using complete MMRF Researcher Gateway (RG) clinical data, 
# you should gain access to MMRF Researcher Gateway (RG) for downloading the following files:
# i)MMRF_CoMMpass_IA14a_All_Canonical_Variants.csv
#ii)MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv



# Once this is done, you should import them into your R environmnet for feeding the following Workflow.
#Note1: pay attention to import MMRF-RG dataset into R environment including the heading of columns.
variant.ann<- MMRF_CoMMpass_IA14a_All_Canonical_Variants
trt<- MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP




names(variant.ann)[52]<-"Gene_name"
trt<-trt[trt$bestresp!="",]

#--------------KCNK2----------------------------

variant.ann1<-variant.ann[(variant.ann$Gene_name=="KCNK2"), ]


summary.var1<-MMRFRG_VariantCountPlot2(variant.ann1,trt,topN=50,filenm=NULL,height=20, width=1500)



#--------------ACTL8-------------------------------------NO Treatment-Resppnse id found

variant.ann2<-variant.ann[variant.ann$Gene_name=="ACTL8", ]


summary.var2<-MMRFRG_VariantCountPlot2(variant.ann2,trt,topN=8,filenm=NULL,height=20, width=1500)



#--------------EPHA5----------------------------

variant.ann3<-variant.ann[(variant.ann$Gene_name=="EPHA5") , ]


summary.var3<-MMRFRG_VariantCountPlot2(variant.ann3,trt,topN=50,filenm=NULL,height=20, width=1500)




#--------------SAGE1----------------------------

variant.ann4<-variant.ann[variant.ann$Gene_name=="SAGE1", ]

summary.var4<-MMRFRG_VariantCountPlot2(variant.ann4,trt,topN=50,filenm=NULL,height=20, width=1500)


#-------------ALL------------------

variant.ann<-variant.ann[(variant.ann$Gene_name=="KCNK2") |
                           (variant.ann$Gene_name=="EPHA5") |
                           (variant.ann$Gene_name=="SAGE1"), ]






summary.var<-MMRFRG_VariantCountPlot2(variant.ann,trt,topN=50,filenm=NULL,height=20, width=1500)

#complete.response<- summary.var[summary.var$bestresp=="Complete Response",]





