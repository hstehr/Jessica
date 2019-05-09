##  Test parameters
rm(list=ls())
data.root="~/Documents/ClinicalDataScience_Fellowship/STAMP_v2.4_reports"
script.root="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts"
outdir.root=paste(data.root,"/trials",sep="")

OnCore.file="~/Documents/ClinicalDataScience_Fellowship/Biomarker_Report/Biomarker_Report_2019-03.csv"
OnCore.ArmRemove="ECOG-ACRIN-EAY131-L,ECOG-ACRIN-EAY131-M"
NCI.file="~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-12-11.xlsx"
## Specify arms that have been officially suspended
NCI.ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

stamp_reference.file=paste(script.root,"/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl.file=paste(script.root,"/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx.key=paste(script.root,"/HistologicalDx_CTEP.csv",sep="")

patient.id="PatientName_UniqueID"
STAMP.file=paste(data.root,"/reports/",patient.id,".variant_report.txt",sep="")
CNV.file=paste(data.root,"/reports/",patient.id,".cnvs",sep="")
Fusion.file=paste(data.root,"/reports/",patient.id,".fusions.filtered.txt",sep="")
