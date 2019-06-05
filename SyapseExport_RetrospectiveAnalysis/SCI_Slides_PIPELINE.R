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

remove(STAMPv2_Annotation,Domains,Tiles,Changes)

#----------------------------------------------
## Replication of patients identified in NCI-MATCH Designated Lab Application: SNV/Indels
#----------------------------------------------
## Need to input "Patient ID information" and MOIs" tabs for parse ARMS
## Switch REF and ALT columns**
## Rename "gene" > "Gene Name" for "Inclusion Non-Hotspot Rules" in ARM-H, ARM-U
## Modifications to HistologicalDx_CTEP.csv: append "adenocarcinoma,lung,,,Non-small cell lung cancer,"

data.root="~/Desktop/"
script.root="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/SyapseExport_RetrospectiveAnalysis/"
outdir.root="~/Desktop/trials_LabApp_Replication"

STAMP.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv"

OnCore.file="FALSE"
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

NCI.file="~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-01-01 EDIT.xlsx"
NCI.ArmRemove="NULL"

stamp_reference_transcripts="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt"
exons_ensembl="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt"
histoDx.key=paste(script.root,"HistologicalDx_CTEP.csv",sep="")

adult_FILTER = "FALSE"
pathogenic_FILTER = "TRUE"
disease_FILTER = "TRUE"

setwd(script.root)
source("ClinicalTrial_Matching_PIPELINE_LabAppRep.R")

## Extract fields for presentation
#----------------------------------------------
VarMatch <- read.csv(file=paste(tempdir,"/NCI_SNVIndel_Variant_Matched_2018-01-01_histologicaldxON_pathogenicON_AdultGroupOFF.tsv",sep=""),sep="\t")
NHSMatch <- read.csv(file = paste(tempdir,"/NCI_SNVIndel_NonHotspot_Matched2018-01-01_histologicaldxON_pathogenicON_AdultGroupOFF.tsv",sep=""),sep="\t")

colname_extract <- c("Arm_Name","PatientID","HistologicalDx","PrimaryTumorSite","VariantLabel","VariantPathogenicityStatus","var.type","VariantHGVSProtein")
colname_extract_new <- c("Arm_Name","Patient_ID","Histological_Dx","Primary_Tumor_Site","Variant_Label","Pathogenicity_Status","Variant_Type","HGVS_Protein")

Final <- rbind(VarMatch[,colname_extract], NHSMatch[,colname_extract])
colnames(Final) <- colname_extract_new
Final <- Final[order(Final$Arm_Name),]

write.table(Final, file = paste(outdir,"Output_Patient_List.tsv",sep=""), 
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)
View(Final)

remove(VarMatch,NHSMatch,Final,colname_extract,colname_extract_new)

# Delete temporary directory
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)}

closeAllConnections()

#----------------------------------------------
## Iteration of difference scenarios 
#----------------------------------------------
data.root="~/Documents/ClinicalDataScience_Fellowship/"
script.root=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/SyapseExport_RetrospectiveAnalysis/",sep="")
outdir.root="~/Desktop/trials_Iterate"

STAMP.file=paste(data.root,"STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv",sep="")

OnCore.file=paste(data.root,"Biomarker_Report/Biomarker_Report_2019-03.csv",sep="")
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

NCI.file=paste(data.root,"PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-22.xlsx",sep="")
NCI.ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

stamp_reference_transcripts=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx.key=paste(script.root,"HistologicalDx_CTEP.csv",sep="")

adult_FILTER = "TRUE"

setwd(script.root)
source("ClinicalTrial_Matching_PIPELINE_Iterate.R")

closeAllConnections()
