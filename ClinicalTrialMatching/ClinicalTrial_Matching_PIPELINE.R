#!/usr/bin/env Rscript

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(Biobase))
suppressMessages(library(eeptools))
suppressMessages(library(splitstackshape))
suppressMessages(library(reshape))
suppressMessages(library(rio))
suppressMessages(library(stringr))
suppressMessages(library(openxlsx))
suppressMessages(library(tidyr))

## Script parameters
#----------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1. Directory location of "REPORTS" folder.
data.root = args[1]
# 2. Patient ID - string format. 
patient.id = args[2]
# 3. File location of SNV/Indel entries.
STAMP.file = args[3]
# 4. File location of CNV entries.
CNV.file = args[4]
# 5. File location of Fusion entries.
Fusion.file = args[5]
# 6. File location of OnCore Report (Stanford OnCore Clinical Trials). To turn off matching, set args == FALSE.
OnCore.file = args[6]
# 7. Names of OnCore Arms to remove
OnCore.ArmRemove = args[7]
# 8. File location of Patient Variant Report (NCI-MATCH Clinical Trials). To turn off matching, set args == FALSE.
NCI.file = args[8]
# 9. Names of NCI-MATCH Arms to remove
NCI.ArmRemove = args[9]
# 10. Directory location of pipeline scripts.
script.root = args[10]
# 11. File location of OUTPUT directory. 
outdir.root = args[11]
# 12. File location of stamp_reference_transcripts file.
stamp_reference.file = args[12]
# 13. File location of exons_ensembl file.
exons_ensembl.file = args[13]
# 14. File location of disease exclusion key file.
histoDx.key = args[14]

## Customize trial output
if (isTRUE(OnCore.file == "FALSE")) {Internal_match <- as.logical("FALSE")
} else {Internal_match <- as.logical("TRUE")}

if (isTRUE(NCI.file == "FALSE")) {NCI_match <- as.logical("FALSE")
} else {NCI_match <- as.logical("TRUE")}

## Directories
#----------------------------------------------
outdir = paste(outdir.root,"/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(outdir.root,"/../temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.output = paste(outdir,patient.id,".out",sep="")
err.output = paste(outdir,patient.id,".err",sep="")

## Load files and specify timestamps
#----------------------------------------------
STAMP_DF <- read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = "\t")

STAMP_CNV <- read.delim(file = CNV.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), 
                        stringsAsFactors = FALSE, sep = "\t", comment.char = '#')

STAMP_Fusion <- read.csv(file = Fusion.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = "\t")

if (isTRUE(Internal_match)) {
  OnCore_Biomarker_Report <- read.csv(file = OnCore.file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, sep = ",")
  
  OnCore_Biomarker_Report_timestamp <- format(as.Date(paste(gsub("([[:digit:]]{4}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*_", "", OnCore.file)), "-01",sep="")), format= "%Y-%m")
}

if (isTRUE(NCI_match)) {
  PATIENT_VARIANT_REPORT <- suppressMessages(import_list(NCI.file, setclass = "tbl"))
  
  Patient_Variant_Report_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1", sub(".*_", "", NCI.file))), format= "%Y-%m-%d")
}

stamp_reference_transcripts <- read.csv(file = stamp_reference.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <- read.csv(file = exons_ensembl.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl, by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)

## Specify script location
setwd(script.root)

## Disease Ontology 
source("DiseaseGroupCategories.R")
source("HistologicalDx_CTEP_Match.R")

## Filter Parameters
#----------------------------------------------
## Customizeable filters
adult.group_FILTER = FALSE # Adult = age >= 18+yo)
pathogenic_FILTER = FALSE
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

## OnCore-specific customizeable filters
disease.group_FILTER = FALSE
disease.site_FILTER = FALSE     # Dependent on disease.group_FILTER == TRUE

## NCI-MATCH-specific customizeable filters
disease.code_FILTER = FALSE

## Filters APPLIED
if (isTRUE(adult.group_FILTER)) {ageName <- "AdultGroupON"} else {ageName <- "AdultGroupOFF"}
if (isTRUE(pathogenic_FILTER)) {pathoName <- "pathogenicON"} else {pathoName <- "pathogenicOFF"}
if (isTRUE(disease.group_FILTER)) {groupName <- "diseaseFILTER_groupON"} else {groupName <- "diseaseFILTER_groupOFF"}
if (isTRUE(disease.group_FILTER & disease.site_FILTER)) {siteName <- "siteON"} else {siteName <- "siteOFF"}
if (isTRUE(disease.code_FILTER)) {dxName <- "histologicaldxON"} else {dxName <- "histologicaldxOFF"}
filterName <- paste(groupName,siteName, "_", dxName, "_", pathoName, "_",  ageName, sep="")

## Generate error file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

## Write output to file
#----------------------------------------------
sink(file = out.output, append = FALSE, split = FALSE)
options(max.print=999999)

## Print parameters to output
#----------------------------------------------
if (isTRUE(Internal_match & NCI_match)) {
  cat("Patient ID:", patient.id, "\n", "\n",
      "Stanford OnCore Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
      "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
      "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
      "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else if (isTRUE(Internal_match)) {
  cat("Patient ID:", patient.id, "\n", "\n",
      "Stanford OnCore Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
      "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else if (isTRUE(NCI_match)) {
  cat("Patient ID:", patient.id, "\n", "\n",
      "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
      "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("NO TRIAL has been specified to be matched. Check parameter specification.", "\n","\n")
  
  sink()
  
  sink(file = out.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
}

## PIPELINE
#----------------------------------------------
# Extract relevant patient data
source("Patient_Export_SNVIndel.R")
source("Patient_Export_CNV.R")
source("Patient_Export_Fusion.R")

# Extract patient_id
patient.list <- unique(unlist(STAMP_DF$PatientID, STAMP_CNV$PatientID, STAMP_Fusion$PatientID))
remove(STAMP_DF,STAMP_Fusion,STAMP_CNV,stamp_reference_full)

if (length(patient.list) > 0) {
  
  if (isTRUE(Internal_match)) {
    source("Biomarker_Report_QC.R")
    
    source("Biomarker_Report_Match_SNVIndel.R")
    source("Biomarker_Report_Match_CNV.R")
    source("Biomarker_Report_Match_Fusion.R") 
    
    remove(OnCore_Biomarker_Report,OnCore_Biomarker_QC)
  }
  
  if (isTRUE(NCI_match)) {
    source("Patient_Variant_Report_QC.R")
    
    source("Patient_Variant_Report_InclusionMatch_SNVIndel.R")
    source("Patient_Variant_Report_InclusionMatch_CNV.R")
    source("Patient_Variant_Report_InclusionMatch_Fusion.R")
    
    source("Patient_Variant_Report_NonHotspotMatch_SNVIndel.R")
    source("Patient_Variant_Report_NonHotspotMatch_CNV.R")
    
    remove(PATIENT_VARIANT_REPORT,Inclusion_NonHotspot_Rules,Inclusion_Variants)
  }
  
  ## Generate OUTPUT files 
  source("ClinicalTrial_Output_Details.R")
  source("ClinicalTrial_Output_tsv.R")
  source("ClinicalTrial_Output_Candidates.R")
  
  remove(Comments,Disease_Exclusion_Codes,Exclusion_NonHotspot_Rules,Exclusion_Variants,IHC_Results)
}

sink()

# Delete temporary directory
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)} 
