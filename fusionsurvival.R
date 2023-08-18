library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(dplyr)
library(survival)
library(ggsurvfit)
library(survminer)
`%nin%` = Negate(`%in%`)

args = commandArgs(trailingOnly = TRUE)
#first arg: study name in all caps (like BRCA); can automate this
study = args[1]
#second arg: directory with chimerseq and studies.txt file that the output can also go to
setwd(args[2])

#organizing chimerseq data
chimerseq = read.table("chimerseq.txt", header=TRUE, sep="\t", row.names=NULL)
chimer = dplyr::select(chimerseq, row.names, Breakpoint.type, Cancer.type, Barcode.ID, Highly_Reliable_Seq)
colnames(chimer)[1] = "Fusion_pair"
rm(chimerseq)

#organizing tcga data
clinical_list = list(ACC.clinical, BLCA.clinical, BRCA.clinical, CESC.clinical, 
                     CHOL.clinical, COAD.clinical, ESCA.clinical, 
                     GBM.clinical, GBMLGG.clinical, HNSC.clinical, KICH.clinical,
                     KIPAN.clinical, KIRC.clinical, KIRP.clinical, LAML.clinical,
                     LGG.clinical, LIHC.clinical, LUAD.clinical, LUSC.clinical,
                     OV.clinical, PAAD.clinical, PCPG.clinical, PRAD.clinical,
                     READ.clinical, SARC.clinical, SKCM.clinical, STAD.clinical,
                     STES.clinical, TGCT.clinical, THCA.clinical, THYM.clinical,
                     UCEC.clinical, UCS.clinical, UVM.clinical)
rnaseq_list = list(ACC.rnaseq, BLCA.rnaseq, BRCA.rnaseq, CESC.rnaseq, 
                   CHOL.rnaseq, COAD.rnaseq, ESCA.rnaseq, 
                   GBM.rnaseq, GBMLGG.rnaseq, HNSC.rnaseq, KICH.rnaseq,
                   KIPAN.rnaseq, KIRC.rnaseq, KIRP.rnaseq, LAML.rnaseq,
                   LGG.rnaseq, LIHC.rnaseq, LUAD.rnaseq, LUSC.rnaseq,
                   OV.rnaseq, PAAD.rnaseq, PCPG.rnaseq, PRAD.rnaseq,
                   READ.rnaseq, SARC.rnaseq, SKCM.rnaseq, STAD.rnaseq,
                   STES.rnaseq, TGCT.rnaseq, THCA.rnaseq, THYM.rnaseq,
                   UCEC.rnaseq, UCS.rnaseq, UVM.rnaseq)
studies = read.table(paste0(args[2], "studies.txt"))

#function to get fusion data for a cancer
fusion_wrangle = function(study){
  cancerchimer = filter(chimer, Cancer.type==study)
  fusiongenes = data.frame(do.call('rbind', strsplit(cancerchimer$Fusion_pair, '-', fixed=TRUE)))
  cancerchimer = cbind(fusiongenes, cancerchimer)
  colnames(cancerchimer)[1] = "Gene1"
  colnames(cancerchimer)[2] = "Gene2"
  bcr_patient_barcode=as.character(cancerchimer[,6])
  for (i in 1:length(bcr_patient_barcode)){
    bcr_patient_barcode[i] = substr(bcr_patient_barcode[i], 1, 12)
  }
  cancerchimer = cbind(bcr_patient_barcode, cancerchimer)
  cancerchimer = cancerchimer[,-7]
  cancerchimer = cancerchimer[!duplicated(cancerchimer), ]
  return(cancerchimer)
}

#function to get tcga survival data for a cancer
tcga_wrangle = function(study){
  #extract data from tcga
  rnaseq = rnaseq_list[[which(studies==study)]] 
  clinicaldf = clinical_list[[which(studies==study)]] 
  survivalTCGA(clinicaldf, barcode.name = "patient.bcr_patient_barcode",
               event.name = "patient.vital_status",
               days.to.followup.name = "patient.days_to_last_followup",
               days.to.death.name = "patient.days_to_death") -> surv
  
  #formatting rnaseq df
  for (i in 1:ncol(rnaseq)){
    colnames(rnaseq)[i] = gsub(r"(\|.*)", "", colnames(rnaseq)[i])
  }
  colnames(rnaseq)[1] = "full_patient_barcode"
  
  #format tcga id for merging
  bcr_patient_barcode=as.character(rnaseq[,1])
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

#combine survival and fusion info for analysis
sigfusion_info = function(study) {
  #wrangling using siggenes
  cancerchimer = fusion_wrangle(study)
  surv_rnaseq = tcga_wrangle(study)
  
  #getting tcga patients with fusions of significant genes
  sigfusion = cancerchimer[cancerchimer$Gene1 %in% genesvec | cancerchimer$Gene2 %in% genesvec,]
  sigfusion = cbind(sigfusion[,1:2], rep(NA,nrow(sigfusion)), rep(NA,nrow(sigfusion)), sigfusion[,3], rep(NA,nrow(sigfusion)), rep(NA,nrow(sigfusion)), sigfusion[,6])
  colnames(sigfusion)[3:8] = c("Gene1.survival", "Gene1.expression", "Gene2", "Gene2.survival", "Gene2.expression", "Cancer")
  
  #merging relevant patients with rnaseq expression data
  #issue: same pt, different samples? which expression to use? solution: average :)
  #i might need to go to jail for this function lowkey
  getexp = function(sigfusion, num){
    if(num == 1){
      value = tryCatch({sigfusion[i,4] = surv_rnaseq[which(surv_rnaseq$bcr_patient_barcode==sigfusion[i,1]), which(colnames(surv_rnaseq)==sigfusion[i,2])]}, 
                       error = function(e) {
                         expdata = surv_rnaseq[which(surv_rnaseq$bcr_patient_barcode==sigfusion[i,1]),
                                               which(colnames(surv_rnaseq)==sigfusion[i,2])]
                         value = mean(expdata)
                       })
    } else if(num == 2){
      value = tryCatch({sigfusion[i,7] = surv_rnaseq[which(surv_rnaseq$bcr_patient_barcode==sigfusion[i,1]),which(colnames(surv_rnaseq)==sigfusion[i,5])]}, 
                       error = function(e) {
                         expdata = surv_rnaseq[which(surv_rnaseq$bcr_patient_barcode==sigfusion[i,1]),
                                               which(colnames(surv_rnaseq)==sigfusion[i,5])]
                         value = mean(expdata)
                       })
    }
    return(value)
  } 
  
  #adding data for better/worse survival and expression
  for(i in 1:nrow(sigfusion)){
    if(sigfusion[i,2]%in%siggenes$gene){
      sigfusion[i,3] = siggenes[which(siggenes$gene==sigfusion[i,2]),3]
      sigfusion[i,4] = getexp(sigfusion, 1)
    }
    if(sigfusion[i,5]%in%siggenes$gene){
      sigfusion[i,6] = siggenes[which(siggenes$gene==sigfusion[i,5]),3]
      sigfusion[i,7] = getexp(sigfusion, 2)
    }
  }

  sigfusion = cbind(sigfusion, rep(NA, nrow(sigfusion)))
  colnames(sigfusion)[9] = "Fusion"
  for (i in 1:nrow(sigfusion)){
    gene1 = sigfusion[i,2]
    gene2 = sigfusion[i,5]
    genesort = sort(c(gene1, gene2))
    sigfusion[i,9] = paste0(genesort[1], "-", genesort[2])
  }
  sigfusionobj = list(sigfusion, surv_rnaseq)
  return(sigfusionobj)
}

surv_analysis = function(){
  sigfusionobj = sigfusion_info(study)
  sigfusion = sigfusionobj[[1]] 
  surv_rnaseq = sigfusionobj[[2]] 
  
  ptsurv = surv_rnaseq[c(1,2,3)]
  ptsurv = ptsurv[!duplicated(ptsurv),]
  rm(sigfusionobj)
  #add steps for desired analysis - maybe more fxns?
}  

for (i in 1:nrow(studies)){
  study=studies[i,]
  print(paste0("study is: ",study))
  geneinfo= read.table(paste0("/Users/lainemarrah/Desktop/LiLab/TCGAgenes/",study,"genes.txt"), row.names = NULL)
  colnames(geneinfo)[1] = "gene"
  siggenes = geneinfo%>%dplyr::select(gene, p.value, survival.corr, cancer)
  genesvec = siggenes$gene
  rm(geneinfo)
  #surv_analysis(study)
  tryCatch({surv_analysis(study)}, error = function(e) {print("No survival-associated fusions available, or another error occurred")})
}



