## source("~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/SyapseExport_RetrospectiveAnalysis/ClinicalTrial_Matching_ARGS.R")
rm(list=ls())

library(rio)
STAMPv2_Annotation <- suppressMessages(import_list(
  "~/Documents/ClinicalDataScience_Fellowship/STAMP/2016-08-23_STAMPv2_Annotation.xlsx",setclass = "tbl"))

# Generate individual DF per worksheet in excel file
invisible(capture.output(lapply(names(STAMPv2_Annotation), 
                                function(x) assign(x, STAMPv2_Annotation[[x]], envir = .GlobalEnv))))

# Re-format protein length INFO 
colnames(Genes) <- unlist(Genes[2,])
colnames_keep <- c("Name","Transcript","Tiles","Size","Coverage","MUT","CNV","Fusions","AA Length")   
Genes <- data.frame(Genes[c(3:nrow(Genes)),colnames_keep])

remove(STAMPv2_Annotation,Domains,Tiles,Changes,colnames_keep)

#----------------------------------------------
## Iteration of difference scenarios 
#----------------------------------------------
data.root="~/Documents/ClinicalDataScience_Fellowship/"
script.root=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/",sep="")
outdir.root="~/Desktop/Trials_Iterate/"

STAMP.file=paste(data.root,"STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv",sep="")

OnCore.file=paste(data.root,"Biomarker_Report/Biomarker_Report_2019-03.csv",sep="")
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

NCI.file=paste(data.root,"PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-22.xlsx",sep="")
NCI.ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

stamp_reference_transcripts=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx.key=paste(script.root,"HistologicalDx_CTEP.tsv",sep="")
AA_key=paste(data.root,"STAMP/AminoAcid_Conversion.csv",sep="")

adult_FILTER = "TRUE"

## Replication of patients identified in NCI-MATCH Designated Lab Application: SNV/Indels
NCI.file_LabRep="~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-01-01 EDIT.xlsx"
NCI.ArmRemove_LabRep="NULL"
OnCore.file_LabRep="FALSE"

setwd(script.root)
source("SyapseExport_RetrospectiveAnalysis/ClinicalTrial_Matching_PIPELINE.R")

closeAllConnections()
