
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(MMRFBiolinks)
library(EnsDb.Hsapiens.v86)


#CD138+ antigen is a plasma cell marker that is highly expressed on plasma cell membranes and it allows 
#an excellent assessment of numbers
#and distribution of plasma cell in bone marrow biopsies \cite{CD138,CD138_2}.

# IMPORTANT 1! If you want to reproduce exaclty the case study you can use embedded data:load("../data/CaseStudy_1.rda")) 

clin.mm<-MMRFGDC_QueryClinic(type = "clinical")

#------------do not delete--------------
query.exp.count <- GDCquery(project = "MMRF-COMMPASS", 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            experimental.strategy = "RNA-Seq",
                            workflow.type="HTSeq - Counts")








#downloading samples
GDCdownload(query.exp.count)


ListBarcode1<-MMRFGDC_QuerySamples(query=query.exp.count,typesample="TBM") 
ListBarcode2<-MMRFGDC_QuerySamples(query=query.exp.count,typesample="TRBM")



ListBarcode1<-setdiff(ListBarcode1,ListBarcode2)

ListSamples<-c(ListBarcode1,ListBarcode2)

query.exp.count.sub <- GDCquery(project = "MMRF-COMMPASS", 
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification",
                                experimental.strategy = "RNA-Seq",
                                workflow.type="HTSeq - Counts",
                                barcode = ListSamples)









GDCdownload(query.exp.count.sub)



MMRFdata.prep.sub <- MMRFGDC_prepare(query.exp.count.sub,
                                     save = TRUE ,
                                     save.filename = "data/RNASeqSE.rda" ,
                                     directory = "GDCdata",
                                     summarizedExperiment = TRUE)






MMRFdataPrepro.sub <- TCGAanalyze_Preprocessing(object = MMRFdata.prep.sub,
                                                cor.cut = 0.6,  # cor.cut = 0,
                                                datatype = "HTSeq - Counts",
                                                filename ="img/MMRF_Preprocessing_OK.png")







df<-assay(MMRFdata.prep.sub)
temp<-df[, unique(colnames(df))]


dataNorm<-TCGAanalyze_Normalization(temp,
                                    geneInfo = geneInfoHT,
                                    method = "gcContent")


dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   




#conversion-----------------

G_list<-rownames(dataFilt)
symbol.gene <- ensembldb::select(EnsDb.Hsapiens.v86, keys= G_list, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(dataFilt)<-symbol.gene$SYMBOL


#-------------------------DEA----------------

dataFilt.TBM<-dataFilt[,colnames(dataFilt) %in% ListBarcode1]
dataFilt.TRBM <- dataFilt[,colnames(dataFilt) %in% ListBarcode2]



dataDEGs <- TCGAanalyze_DEA(dataFilt.TBM,
                            dataFilt.TRBM,
                            metadata =FALSE,
                            Cond1type = "TBM ",
                            Cond2type = "TRBM")






dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,
                                          "TBM",
                                          "TRB",
                                          dataFilt.TBM,
                                          dataFilt.TRBM)




dataDEGsFiltLevel.rank<-dataDEGsFiltLevel[abs(dataDEGsFiltLevel$logFC)>3,]
dataDEGsFiltLevel.rank<-dataDEGsFiltLevel.rank[dataDEGsFiltLevel.rank$FDR<0.05,]
dataDEGsFiltLevel.rank<-dataDEGsFiltLevel.rank %>% arrange(FDR, -abs(logFC)) 
rownames(dataDEGsFiltLevel.rank)<-dataDEGsFiltLevel.rank$mRNA


TCGAVisualize_volcano(dataDEGsFiltLevel.rank$logFC, dataDEGsFiltLevel.rank$FDR,
                      filename = "img/DEGs_volcano_OK2.png",
                      x.cut =1,
                      y.cut = 0.05,
                      names = rownames(dataDEGsFiltLevel.rank),
                      color = c("black","red","dodgerblue3"),
                      names.size = 2,
                      show.names = "highlighted",
                      highlight = c("gene1","gene2"),
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (TBM sample vs TRBM sample)",
                      width = 10)





#-------------------------PCA----------------

group1<-ListBarcode1[ListBarcode1 %in% colnames(dataFilt)] 
group2<-ListBarcode2[ListBarcode2 %in% colnames(dataFilt)] 





pca <- TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel.rank, ntopgenes =10,group1, group2)#ggbiplot with varname.size=6

#group1=blue TBM
#group2=red TRBM













#-------------Enrichnment Analysis----------------------------------------

Genelist <- rownames(dataDEGsFiltLevel.rank)

system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA TBM sample vs TRBM sample",Genelist))

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist, 
                        nBar = 5,
                        text.size = 2.4)










#---------------------------------------------------PREPARE SURVIVAL-------------





dataMMcomplete <- log2(dataFilt)






#----------------------------------------------------------





TabSurvKM <- TCGAanalyze_SurvivalKM(clin.mm,
                                dataMMcomplete,
                                Genelist = rownames(dataDEGsFiltLevel.rank),
                                Survresult = FALSE,
                                ThreshTop = 0.67,
                                ThreshDown = 0.33,
                                p.cut = 0.05)







TabSurvKM <- TabSurvKM[order(TabSurvKM$pvalue, decreasing=F),]



surv.genes<-rownames(TabSurvKM)



#--------------------

TabSurvKM.gene1 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = surv.genes[1],
                                          Survresult = TRUE,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          p.cut = 0.05)











TabSurvKM.gene2 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = surv.genes[2],
                                          Survresult = TRUE,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          p.cut = 0.05)


TabSurvKM.gene3 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = surv.genes[3],
                                          Survresult = TRUE,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          p.cut = 0.05)


TabSurvKM.gene4 <- TCGAanalyze_SurvivalKM(clin.mm,
                                          dataMMcomplete,
                                          Genelist = surv.genes[4],
                                          Survresult = TRUE,
                                          ThreshTop = 0.67,
                                          ThreshDown = 0.33,
                                          p.cut = 0.05)
