
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





TCGAVisualize_volcano(dataDEGs$logFC, dataDEGs$FDR,
                      filename = "img/DEGs_volcano_OK2.png",
                      x.cut =1,
                      y.cut = 0.05,
                      names = rownames(dataDEGs),
                      color = c("black","red","dodgerblue3"),
                      names.size = 2,
                      show.names = "highlighted",
                      highlight = c("gene1","gene2"),
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (TBM sample vs TRBM sample)",
                      width = 10)










#Enrichnment Analysis----------------------------------------

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
                        nBar = 10,
                        text.size = 1.5)










#---------------SURVIVAL--------------------------------------


dataMMcomplete <- log2(dataFilt)
dataMMcomplete<-dataFilt













TabSurvKM <- MMRFGDC_SurvivalKM(clin.mm,
                                dataMMcomplete,
                                ListGenes = rownames(dataDEGs),
                                Results = FALSE,
                                High = 0.67,
                                Low = 0.33,
                                p.cut = 0.05,
                                ListBarcode1,
                                ListBarcode2)


TabSurvKM <- TabSurvKM[order(TabSurvKM$pvalue, decreasing=F),]


TabSurvKM.10<-TabSurvKM[1:10,]
surv.genes<-rownames(TabSurvKM.10)





TabSurvKM.gene1 <- TCGAanalyze_SurvivalKM2(clin.mm,
                                           dataMMcomplete,
                                           ListGenes = surv.genes[1],
                                           Results = TRUE,
                                           High = 0.67,
                                           Low = 0.33,
                                           p.cut = 0.05,
                                           ListBarcode1,
                                           ListBarcode2)



TabSurvKM.gene2 <- TCGAanalyze_SurvivalKM2(clin.mm,
                                           dataMMcomplete,
                                           ListGenes = surv.genes[2],
                                           Results = TRUE,
                                           High = 0.67,
                                           Low = 0.33,
                                           p.cut = 0.05,
                                           ListBarcode1,
                                           ListBarcode2)


TabSurvKM.gene3 <- TCGAanalyze_SurvivalKM2(clin.mm,
                                           dataMMcomplete,
                                           ListGenes = surv.genes[3],
                                           Results = TRUE,
                                           High = 0.67,
                                           Low = 0.33,
                                           p.cut = 0.05,
                                           ListBarcode1,
                                           ListBarcode2)



TabSurvKM.gene4 <- TCGAanalyze_SurvivalKM2(clin.mm,
                                           dataMMcomplete,
                                           ListGenes = surv.genes[4],
                                           Results = TRUE,
                                           High = 0.67,
                                           Low = 0.33,
                                           p.cut = 0.05,
                                           ListBarcode1,
                                           ListBarcode2)




#-------------------------PCA----------------

group1<-ListBarcode1[ListBarcode1 %in% colnames(dataFilt)] 
group2<-ListBarcode2[ListBarcode2 %in% colnames(dataFilt)] 


pca <- TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel.rank, ntopgenes = 10, ListBarcode1, ListBarcode2)



















