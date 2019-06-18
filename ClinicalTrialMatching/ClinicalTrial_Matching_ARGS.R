rm(list=ls())

##  Test parameters
#----------------------------------------------
# Directory location of "REPORTS" folder.
data.root = "/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/STAMP_v2.4_reports"
# Directory location of pipeline scripts.
script.root = "/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts"
# File location of OUTPUT directory. 
outdir.root = paste(data.root,"/trials",sep="")

# File location of OnCore Report (Stanford OnCore Clinical Trials). To turn off matching, set args == FALSE.
OnCore.file = "~/Documents/ClinicalDataScience_Fellowship/Biomarker_Report/Biomarker_Report_2019-03.csv"
# Names of OnCore Arms to remove
OnCore.ArmRemove = "ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"

# File location of Patient Variant Report (NCI-MATCH Clinical Trials). To turn off matching, set args == FALSE.
NCI.file = "~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-12-11.xlsx"
# Names of NCI-MATCH Arms to remove
NCI.ArmRemove = "ARM-Z1C"

# File location of stamp_reference_transcripts file.
stamp_reference.file = paste(script.root,"/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
# File location of exons_ensembl file.
exons_ensembl.file = paste(script.root,"/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
# File location of disease exclusion key file.
histoDx.key = paste(script.root,"/HistologicalDx_CTEP.tsv",sep="")
# File location of amino acid 3-letter / 1-letter conversion key
AA.key=paste(data.root,"ClinicalTrialMatching/PIPELINE_scripts/AminoAcid_Conversion.csv",sep="")

# Patient ID - string format. 
patient.id = "PatientName_UniqueID"
# File location of SNV/Indel entries.
STAMP.file = paste(data.root,"/reports/",patient.id,".variant_report.txt",sep="")
# File location of CNV entries.
CNV.file = paste(data.root,"/reports/",patient.id,".cnvs",sep="")
# File location of Fusion entries.
Fusion.file = paste(data.root,"/reports/",patient.id,".fusions.filtered.txt",sep="")
