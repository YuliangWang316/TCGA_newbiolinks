library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(EDASeq)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)
library(tidyverse)
library(GSVA)
library(IDConverter)
request_cancer=c("ACC","CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS")#,"CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS"
i<-"LIHC"
cancer_type=paste("TCGA",i,sep="-")
print(cancer_type)
clinical<-GDCquery_clinic(project = cancer_type, type = "clinical")
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
count_matrix<-as.data.frame(assay(expdat))
count_gc<-TCGAanalyze_Normalization(count_matrix, geneInfoHT,method =  'gcContent')
count_gl<-TCGAanalyze_Normalization(count_matrix, geneInfoHT,method =  'geneLength')
remove(count_matrix,expdat,query)
genename<-rownames(count_gl)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
count_gl<-count_gl[e$ENSEMBL,]
rownames(count_gl)<-e$SYMBOL
remove(e,genename)
genename<-rownames(count_gc)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
count_gc<-count_gc[e$ENSEMBL,]
rownames(count_gc)<-e$SYMBOL
remove(e,genename)
samplesTP <- TCGAquery_SampleTypes(colnames(count_gc), typesample = c("TP"))
count_gc_new <- as.data.frame(count_gc["FOXP3",samplesTP])
remove(samplesTP,count_gc)
samplesTP <- TCGAquery_SampleTypes(colnames(count_gl), typesample = c("TP"))
count_gl_new <- as.data.frame(count_gl["FOXP3",samplesTP])
remove(samplesTP,count_gl)
colnames(count_gc_new)<-"FOXP3"
colnames(count_gl_new)<-"FOXP3"
count_gc_new$sample <- sapply(strsplit(rownames(count_gc_new),'-'),function(x) paste0(x[1:3],collapse="-"))
count_gl_new$sample <- sapply(strsplit(rownames(count_gl_new),'-'),function(x) paste0(x[1:3],collapse="-"))
newname_gl<-filter_tcga_barcodes(rownames(count_gl_new),analyte_target = "RNA")
newname_gc<-filter_tcga_barcodes(rownames(count_gc_new),analyte_target = "RNA")
count_gc_new_new<-count_gc_new[newname_gc,]
count_gl_new_new<-count_gl_new[newname_gl,]
remove(count_gc_new,count_gl_new,newname_gc,newname_gl)
clinical$"FOXP3_gl" <- count_gl_new_new[match(clinical$submitter_id,count_gl_new_new$sample),][,"FOXP3"]
clinical$"FOXP3_gc" <- count_gc_new_new[match(clinical$submitter_id,count_gc_new_new$sample),][,"FOXP3"]
remove(count_gc_new_new,count_gl_new_new)
df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,FOXP3_gl,FOXP3_gc,age_at_index,gender,ajcc_pathologic_stage))
df <- df[!is.na(df$FOXP3_gc),]
df<-df[which(df$vital_status!="NA"),]
for (j in 1:length(rownames(df))) {
  if(is.na(df$days_to_death[j])){
    df$Time[j] <- df$days_to_last_follow_up[j]
  }else if(is.na(df$days_to_last_follow_up[j]) ){
    df$Time[j] <- df$days_to_death[j]
  }
  else if(df$days_to_death[j] >=df$days_to_last_follow_up[j]){
    df$Time[j] <-df$days_to_death[j]
  }
}
df<-df[which(df$Time != 0),]
for (j in 1:length(rownames(df))) {
  if(df$vital_status[j] == "Alive"){
    df$events[j]<-0
  }else if(df$vital_status[j] == "Dead"){
    df$events[j]<-1
  }
}
res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = "FOXP3_gc" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(Time,events)~ FOXP3_gc,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(Time,events) ~ FOXP3_gc,data=df)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)
