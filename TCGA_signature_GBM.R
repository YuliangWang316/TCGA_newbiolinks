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
genset<-read.table("c:/Users/xjmik/Downloads/Treg_1C_signature.txt",sep = "\t",header = TRUE)
geneset_list<-as.list(genset)
request_cancer=c("ACC","CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS")#,"CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS"
for (i in request_cancer) {
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
  gsva_gl<-gsva(expr = count_gl,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
  genename<-rownames(count_gc)
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  e<-e[!duplicated(e$SYMBOL),]
  count_gc<-count_gc[e$ENSEMBL,]
  rownames(count_gc)<-e$SYMBOL
  remove(e,genename)
  gsva_gc<-gsva(expr = count_gc,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
  gsva_gl_2<-gsva(expr = count_gl,gset.idx.list = geneset_list,kcdf="Gaussian",parallel.sz=20)
  gsva_gc_2<-gsva(expr = count_gc,gset.idx.list = geneset_list,kcdf="Gaussian",parallel.sz=20)
  gsva_gl_3<-gsva(expr = count_gl,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="ssgsea")
  gsva_gc_3<-gsva(expr = count_gc,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="ssgsea")
  # gsva_gl_4<-gsva(expr = count_gl,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="zscore")
  # gsva_gc_4<-gsva(expr = count_gc,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="zscore")
  # gsva_gl_5<-gsva(expr = count_gl,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="plage")
  # gsva_gc_5<-gsva(expr = count_gc,gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20,method="plage")
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gc), typesample = c("TP"))
  gsva_gc_new <- as.data.frame(gsva_gc["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gc)
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gl), typesample = c("TP"))
  gsva_gl_new <- as.data.frame(gsva_gl["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gl)
  remove(count_gc,count_gl,geneset_list,genset)
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gc_2), typesample = c("TP"))
  gsva_gc_2_new <- as.data.frame(gsva_gc_2["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gc_2)
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gc_3), typesample = c("TP"))
  gsva_gc_3_new <- as.data.frame(gsva_gc_3["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gc_3)
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gl_2), typesample = c("TP"))
  gsva_gl_2_new <- as.data.frame(gsva_gl_2["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gl_2)
  samplesTP <- TCGAquery_SampleTypes(colnames(gsva_gl_3), typesample = c("TP"))
  gsva_gl_3_new <- as.data.frame(gsva_gl_3["Treg_1C",samplesTP])
  remove(samplesTP,gsva_gl_3)
  colnames(gsva_gc_new)<-"Treg_1C"
  colnames(gsva_gl_new)<-"Treg_1C"
  colnames(gsva_gc_2_new)<-"Treg_1C"
  colnames(gsva_gl_2_new)<-"Treg_1C"
  colnames(gsva_gc_3_new)<-"Treg_1C"
  colnames(gsva_gl_3_new)<-"Treg_1C"
  gsva_gc_new$sample <- sapply(strsplit(rownames(gsva_gc_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  gsva_gl_new$sample <- sapply(strsplit(rownames(gsva_gl_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  gsva_gc_2_new$sample <- sapply(strsplit(rownames(gsva_gc_2_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  gsva_gl_2_new$sample <- sapply(strsplit(rownames(gsva_gl_2_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  gsva_gc_3_new$sample <- sapply(strsplit(rownames(gsva_gc_3_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  gsva_gl_3_new$sample <- sapply(strsplit(rownames(gsva_gl_3_new),'-'),function(x) paste0(x[1:3],collapse="-"))
  newname_gl<-filter_tcga_barcodes(rownames(gsva_gl_new),analyte_target = "RNA")
  newname_gc<-filter_tcga_barcodes(rownames(gsva_gc_new),analyte_target = "RNA")
  newname_gl_2<-filter_tcga_barcodes(rownames(gsva_gl_2_new),analyte_target = "RNA")
  newname_gc_2<-filter_tcga_barcodes(rownames(gsva_gc_2_new),analyte_target = "RNA")
  newname_gl_3<-filter_tcga_barcodes(rownames(gsva_gl_3_new),analyte_target = "RNA")
  newname_gc_3<-filter_tcga_barcodes(rownames(gsva_gc_3_new),analyte_target = "RNA")
  gsva_gc_new_new<-gsva_gc_new[newname_gc,]
  gsva_gl_new_new<-gsva_gl_new[newname_gl,]
  gsva_gc_2_new_new<-gsva_gc_2_new[newname_gc_2,]
  gsva_gl_2_new_new<-gsva_gl_2_new[newname_gl_2,]
  gsva_gc_3_new_new<-gsva_gc_3_new[newname_gc_3,]
  gsva_gl_3_new_new<-gsva_gl_3_new[newname_gl_3,]
  remove(gsva_gc_new,gsva_gl_new,newname_gc,newname_gl)
  remove(gsva_gc_2_new,gsva_gl_2_new,newname_gc_2,newname_gl_2)
  remove(gsva_gc_3_new,gsva_gl_3_new,newname_gc_3,newname_gl_3)
  clinical$"Treg_1C_gl" <- gsva_gl_new_new[match(clinical$submitter_id,gsva_gl_new_new$sample),][,"Treg_1C"]
  clinical$"Treg_1C_gc" <- gsva_gc_new_new[match(clinical$submitter_id,gsva_gc_new_new$sample),][,"Treg_1C"]
  clinical$"Treg_1C_gl_2" <- gsva_gl_2_new_new[match(clinical$submitter_id,gsva_gl_2_new_new$sample),][,"Treg_1C"]
  clinical$"Treg_1C_gc_2" <- gsva_gc_2_new_new[match(clinical$submitter_id,gsva_gc_2_new_new$sample),][,"Treg_1C"]
  clinical$"Treg_1C_gl_3" <- gsva_gl_3_new_new[match(clinical$submitter_id,gsva_gl_3_new_new$sample),][,"Treg_1C"]
  clinical$"Treg_1C_gc_3" <- gsva_gc_3_new_new[match(clinical$submitter_id,gsva_gc_3_new_new$sample),][,"Treg_1C"]
  remove(gsva_gc_new_new,gsva_gl_new_new,gsva_gl_2_new_new,gsva_gc_2_new_new,gsva_gl_3_new_new,gsva_gc_3_new_new)
  df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,Treg_1C_gl,Treg_1C_gc,age_at_index,gender,ajcc_pathologic_stage,Treg_1C_gl_2,Treg_1C_gc_2,Treg_1C_gl_3,Treg_1C_gc_3))
  df <- df[!is.na(df$Treg_1C_gc),]
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
  res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = "Treg_1C_gc" )
  summary(res.cut)
  res.cat<-surv_categorize(res.cut)
  fit<-survfit(Surv(Time,events)~ Treg_1C_gc,data = res.cat)
  ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
  res.cox<-coxph(Surv(Time,events) ~ Treg_1C_gc,data=res.cat)
  summary(res.cox)
  test.ph<-cox.zph(res.cox)
  ggcoxzph(test.ph)
  newrescat<-cbind(res.cat,df$age_at_index,df$gender,df$ajcc_pathologic_stage)
  # res.cox_new<-coxph(Surv(Time,events) ~ Treg_1C_gc + df$age_at_index + df$gender + df$ajcc_pathologic_stage,data=newrescat)
  colnames(newrescat)[4:6]<-c("Age","gender","stage")
  res.cox_new<-coxph(Surv(Time,events) ~ Treg_1C_gc + Age ,data=newrescat)
  summary(res.cox_new)
  test.ph_new<-cox.zph(res.cox_new)
  ggcoxzph(test.ph_new)
  
  # for (j in 1:length(rownames(newrescat))) {
  #   if(newrescat$gender[j] == "male"){
  #     newrescat$gender_new[j]<-0
  #   }else if(newrescat$gender[j] == "female"){
  #     newrescat$gender_new[j]<-1
  #   }
  # }
  # newrescat<-newrescat[which(!is.na(newrescat$stage)),]
  # for (j in 1:length(rownames(newrescat))) {
  #    if(newrescat$stage[j] == "Stage IA"){
  #     newrescat$stage_new[j]<-0
  #   }else if(newrescat$stage[j] == "Stage IB"){
  #     newrescat$stage_new[j]<-1
  #   }else if(newrescat$stage[j] == "Stage II"){
  #     newrescat$stage_new[j]<-2.5
  #   }else if(newrescat$stage[j] == "Stage IIA"){
  #     newrescat$stage_new[j]<-2
  #   }else if(newrescat$stage[j] == "Stage IIB"){
  #     newrescat$stage_new[j]<-3
  #   }else if(newrescat$stage[j] == "Stage IIIA"){
  #     newrescat$stage_new[j]<-4
  #   }else if(newrescat$stage[j] == "Stage IIIB"){
  #     newrescat$stage_new[j]<-5
  #   }else if(newrescat$stage[j] == "Stage IV"){
  #     newrescat$stage_new[j]<-6
  #   }
  # }
  # res.cox_new_new<-coxph(Surv(Time,events) ~ Treg_1C_gc + Age + gender_new + stage_new+tt(stage_new)+tt(gender_new),data=newrescat,
  #                tt=function(x,t,...) {x*log(t)}
  # )
  # summary(res.cox_new_new)
  # res.cox_new_new_new<-coxph(Surv(Time,events) ~ Treg_1C_gc + Age + gender + stage,data=newrescat)
  ggadjustedcurves(res.cox_new,data = newrescat,variable = "Treg_1C_gc",risk.table = TRUE,pval = TRUE)
}








