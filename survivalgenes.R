library(RegParallel)
library(e1071)
library(caret)
library(tidyverse)
library(RTCGA.clinical) # survival times
library(RTCGA.rnaseq) # genes' expression
library(survminer)
library(survival)
library(ggfortify)
set.seed(28898)

clinical_list = list(ACC.clinical, BLCA.clinical, BRCA.clinical, CESC.clinical, 
                     CHOL.clinical, COAD.clinical, DLBC.clinical, ESCA.clinical, 
                     GBM.clinical, GBMLGG.clinical, HNSC.clinical, KICH.clinical,
                     KIPAN.clinical, KIRC.clinical, KIRP.clinical, LAML.clinical,
                     LGG.clinical, LIHC.clinical, LUAD.clinical, LUSC.clinical,
                     OV.clinical, PAAD.clinical, PCPG.clinical, PRAD.clinical,
                     READ.clinical, SARC.clinical, SKCM.clinical, STAD.clinical,
                     STES.clinical, TGCT.clinical, THCA.clinical, THYM.clinical,
                     UCEC.clinical, UCS.clinical, UVM.clinical)
rnaseq_list = list(ACC.rnaseq, BLCA.rnaseq, BRCA.rnaseq, CESC.rnaseq, 
                     CHOL.rnaseq, COAD.rnaseq, DLBC.rnaseq, ESCA.rnaseq, 
                     GBM.clinical, GBMLGG.rnaseq, HNSC.rnaseq, KICH.rnaseq,
                     KIPAN.rnaseq, KIRC.rnaseq, KIRP.rnaseq, LAML.rnaseq,
                     LGG.rnaseq, LIHC.rnaseq, LUAD.rnaseq, LUSC.rnaseq,
                     OV.rnaseq, PAAD.rnaseq, PCPG.rnaseq, PRAD.rnaseq,
                     READ.rnaseq, SARC.rnaseq, SKCM.rnaseq, STAD.rnaseq,
                     STES.rnaseq, TGCT.rnaseq, THCA.rnaseq, THYM.rnaseq,
                     UCEC.rnaseq, UCS.rnaseq, UVM.rnaseq)
studies = read.table("studies.txt")

#positional arg: study name in all caps (like BRCA); can automate this
args = commandArgs(trailingOnly = TRUE)
study = args[1]
clinicaldf = clinical_list[[which(studies==study)]] 
survivalTCGA(clinicaldf, barcode.name = "patient.bcr_patient_barcode",
             event.name = "patient.vital_status",
             days.to.followup.name = "patient.days_to_last_followup",
             days.to.death.name = "patient.days_to_death") -> surv
rnaseq = rnaseq_list[[which(studies==study)]] 
rm(clinical_list)
rm(rnaseq_list)

#a function to format tcga survival/clinical data
tcga_wrangle = function(surv, rnaseq){
  #formatting rnaseq df
  #remove numeric ids from columns
  for (i in 1:ncol(rnaseq)){
    colnames(rnaseq)[i] = gsub(r"(\|.*)", "", colnames(rnaseq)[i])
  }
  colnames(rnaseq)[1] = "full_patient_barcode"
  
  #format tcga id for merging
  bcr_patient_barcode=rnaseq[,1]
  bcr_patient_barcode=as.character(bcr_patient_barcode)
  for (i in 1:length(bcr_patient_barcode)){
    bcr_patient_barcode[i] = substr(bcr_patient_barcode[i], 1, 12)
  }
  rnaseq=cbind(bcr_patient_barcode, rnaseq)
  
  #merge dfs and fix hyphens in column names
  surv_rnaseq = merge(surv, rnaseq, by="bcr_patient_barcode")
  surv_rnaseq = surv_rnaseq %>% rename(LST3TM12=`LST-3TM12`)
  surv_rnaseq = surv_rnaseq %>% rename(RP1177G6.2=`RP1-177G6.2`)
  return(surv_rnaseq)
}

#a function to generate all possible cox models
#can ignore warnings from coxph.fit step, suppressing doesn't work
cox_models = function(surv, rnaseq){
  surv_rnaseq = tcga_wrangle(surv, rnaseq)
  
  #finding all univariate cox ph models
  genes=colnames(surv_rnaseq[which(colnames(surv_rnaseq)=="A1BG"):ncol(surv_rnaseq)])
  survobj = Surv(surv_rnaseq$times, surv_rnaseq$patient.vital_status)
  coxfun = function(x){coxph(survobj ~ x)}
  univ_models = vector(mode='list', length=length(genes))
  #have to add constant to go between columns and gene list
  constant=(which(colnames(surv_rnaseq)=="A1BG"))-1
  #can ignore warnings
  for (i in 1:length(genes)){
    model = coxfun(surv_rnaseq[,i+constant])
    univ_models[[i]] = model
    univ_models[[i]]["gene"]=genes[i]
  }
  options(scipen=999)
  #getting rid of NA coefficient models
  nonas_models = vector(mode='list', length=length(univ_models))
  for (i in length(univ_models):1){ #backwards to not mess with indices
    if (is.na(univ_models[[i]][["coefficients"]][["x"]])){
      nonas_models[i] = NULL
    } else {
      nonas_models[i] = univ_models[i]
    }
  }
  rm(univ_models)
  rm(coxfun)
  return(nonas_models)
}
  
#a function to evaluate assumptions of cox-ph models based on this website
#http://www.sthda.com/english/wiki/cox-model-assumptions
cox_assumptions = function(surv, rnaseq){
  surv_rnaseq = tcga_wrangle(surv, rnaseq)
  nonas_models = cox_models(surv, rnaseq)
  #proportional hazards assumption
  prophazards = vector(mode='list', length=length(nonas_models))
  survobj = Surv(surv_rnaseq$times, surv_rnaseq$patient.vital_status)
  for (i in length(nonas_models):1){
    tbl = cox.zph(nonas_models[[i]])
    if (is.na(tbl[["table"]][1,3])) {
      prophazards[i] = NULL
    } else if (tbl[["table"]][1,3]>0.05) {
      prophazards[i] = NULL
    } else {
      prophazards[i] = nonas_models[i]
    }
  }
  rm(tbl)
  
  #checking linearity assumption
  linear = vector(mode='list', length=length(prophazards))
  for (i in length(prophazards):1){
    gene = prophazards[[i]][["gene"]]
    genecol = which(colnames(surv_rnaseq)==gene)
    plot = ggcoxfunctional(Surv(times, patient.vital_status) ~ surv_rnaseq[,genecol], data = surv_rnaseq)
    data = plot[["surv_rnaseq[, genecol]"]]
    plotdf=data.frame(plot[["surv_rnaseq[, genecol]"]][["data"]][["explanatory"]], plot[["surv_rnaseq[, genecol]"]][["data"]][["martingale_resid"]])
    colnames(plotdf) = c("Gene_values", "Martingale_residuals")
    lm = lm(data=plotdf, Martingale_residuals~Gene_values)
    p = anova(lm)[1,5]
    if (p<0.05){
      linear[i] = prophazards[i]
    } else {
      linear[i] = NULL
    }
  }
  return(linear)
  rm(gene)
  rm(prophazards)
  rm(genecol)
  rm(plot)
  rm(data)
  rm(plotdf)
  rm(lm)
}

#can ignore warnings about deprecated feature
cox_results = function(surv, rnaseq){
  models = cox_assumptions(surv, rnaseq) #try to shorten runtime
  #checking various values
  univ_results = lapply(models,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coefficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })
  genelist = vector(length=length(models))
  for (i in 1:length(models)){
    g = as.character(models[[i]]["gene"])
    genelist[i] = g
    rm(g)
  }
  #extracting and organizing significant genes
  res_surv = t(as.data.frame(univ_results, check.names = FALSE))
  row.names(res_surv) = genelist
  res_surv = as.data.frame(res_surv)
  res_survSig = subset(res_surv, p.value < 0.01)
  siggenes = res_survSig[order(res_survSig$p.value),]
  return(siggenes)
}

#running fxns and outputting results
cancergenes = cox_results(surv, rnaseq)
filename = paste0(study,"genes.txt")
write.table(cancergenes, filename)

#cancer info: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations 


