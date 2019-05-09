##  Test parameters
rm(list=ls())
data.root="~/Documents/ClinicalDataScience_Fellowship/"
script.root=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/SyapseExport_RetrospectiveAnalysis_20190507/",sep="")
outdir.root="~/Desktop/trials"

STAMP.file=paste(data.root,"STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv",sep="")

OnCore.file=paste(data.root,"Biomarker_Report/Biomarker_Report_2019-03.csv",sep="")
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

NCI.file=paste(data.root,"PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-22.xlsx",sep="")
NCI.ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

stamp_reference_transcripts=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx.key=paste(script.root,"HistologicalDx_CTEP.csv",sep="")


#----------------------------------------------
## Replication of patients identified in NCI-MATCH Designated Lab Application: SNV/Indels
#----------------------------------------------
## Need to input "Patient ID information" and MOIs" tabs for parse ARMS
## Switch REF and ALT columns**
## Rename "gene" > "Gene Name" for "Inclusion Non-Hotspot Rules" in ARM-H, ARM-U
## Modifications to HistologicalDx_CTEP.csv: append "adenocarcinoma,lung,,,Non-small cell lung cancer,"
## adult.group_FILTER == FALSE

rm(list=ls())
data.root="~/Desktop/"
script.root="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/SyapseExport_RetrospectiveAnalysis_20190507/"
outdir.root="~/Desktop/trials"

STAMP.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv"

OnCore.file="FALSE"
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

NCI.file="~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-01-01 EDIT.xlsx"
NCI.ArmRemove="NULL"

stamp_reference_transcripts="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt"
exons_ensembl="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt"
histoDx.key=paste(script.root,"HistologicalDx_CTEP.csv",sep="")

## Extract fields for presentation
#----------------------------------------------
VarMatch <- read.csv(file="~/Desktop/temp/NCI_SNVIndel_Variant_Matched_2018-01-01_histologicaldxON_pathogenicON_AdultGroupOFF.tsv",sep="\t")
NHSMatch <- read.csv(file = "~/Desktop/temp/NCI_SNVIndel_NonHotspot_Matched2018-01-01_histologicaldxON_pathogenicON_AdultGroupOFF.tsv",sep="\t")

colname_extract <- c("Arm_Name","PatientID","HistologicalDx","PrimaryTumorSite","VariantLabel","VariantPathogenicityStatus","var.type","VariantHGVSProtein")
colname_extract_new <- c("Arm_Name","Patient_ID","Histological_Dx","Primary_Tumor_Site","Variant_Label","Pathogenicity_Status","Variant_Type","HGVS_Protein")

Final <- rbind(VarMatch[,colname_extract], NHSMatch[,colname_extract])
colnames(Final) <- colname_extract_new
Final <- Final[order(Final$Arm_Name),]

View(Final)
write.table(Final, file = "~/Desktop/Output_Patient_List.tsv", 
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)

remove(VarMatch,NHSMatch,Final,colname_extract,colname_extract_new)
