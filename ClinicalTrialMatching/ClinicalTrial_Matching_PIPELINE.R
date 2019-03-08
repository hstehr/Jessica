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

args = commandArgs(trailingOnly=TRUE)
# 1. Directory of "reports" folder.
data.root = args[1]
# 2. Location of STAMP entries.
STAMP.file = args[2]
# 3. Patient ID.
patient.id = args[3]
# 4. Location of OnCore Report (Stanford Internal Clinical Trials). To turn off matching, set args == FALSE.
OnCore.file = args[4]
# 5. Location of Patient Variant Report (NCI-MATCH Clinical Trials). To turn off matching, set args == FALSE.
NCI.file = args[5]
# 6. Directory of pipeline scripts.
script.root = args[6]
# 7. Directory to save output files to. 
outdir.root = args[7]
# 8. Location of stamp_reference_transcripts.
stamp_reference.file = args[8]
# 9. Location of exons_ensembl.
exons_ensembl.file = args[9]
# 10. Location of disease exclusion key.
histoDx.key = args[10]

## Customize trial output
if (isTRUE(OnCore.file == "FALSE")) {
  Internal_match = FALSE  
} else {
  Internal_match = TRUE
}

if (isTRUE(NCI.file == "FALSE")) {
  NCI_match = FALSE  
} else {
  NCI_match = TRUE
}

## Directories
outdir = outdir.root
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(outdir.root,"/../temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.ouput = paste(outdir,patient.id,".out",sep="")
err.output = paste(outdir,patient.id,".err",sep="")

## Load files and specify timestamps
#----------------------------------------------
STAMP_DF <- 
  read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = "\t")

if (isTRUE(Internal_match)) {
  OnCore_Biomarker_Report <- 
    read.csv(file = OnCore.file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, sep = ",")

    OnCore_Biomarker_Report_timestamp <- 
    format(as.Date(paste(gsub("([[:digit:]]{4}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*_", "", OnCore.file)), 
                         "-01",sep="")), format= "%Y-%m")
}

if (isTRUE(NCI_match)) {
  PATIENT_VARIANT_REPORT <- suppressMessages(import_list(NCI.file, setclass = "tbl"))
  
  Patient_Variant_Report_timestamp <- 
    format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",
                        sub(".*_", "", NCI.file))), format= "%Y-%m-%d")
}

stamp_reference_transcripts <- 
  read.csv(file = stamp_reference.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <-
  read.csv(file = exons_ensembl.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl,
                                  by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)

## Specify script location
setwd(script.root)

## Disease Ontology 
source("DiseaseGroupCategories.R")
source("HistologicalDx_CTEP_Match.R")

## Default filters 
adult.group_FILTER = FALSE # Adult = age >= 18+yo)
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

## Customizeable filters
pathogenic_FILTER = FALSE
disease.group_FILTER = FALSE
disease.site_FILTER = FALSE # Dependent on disease.group_FILTER == TRUE
disease.code_FILTER = FALSE 

AgePlot.FILTER = FALSE
VariantPlot.FILTER = FALSE

## Filters APPLIED
if (isTRUE(adult.group_FILTER)) { ageName <- "AdultGroupON" } else { ageName <- "AdultGroupOFF" }
if (isTRUE(pathogenic_FILTER)) { pathoName <- "pathogenicON" } else { pathoName <- "pathogenicOFF" }
if (isTRUE(disease.group_FILTER)) { groupName <- "diseaseFILTER_groupON" } else { groupName <- "diseaseFILTER_groupOFF" }
if (isTRUE(disease.group_FILTER & disease.site_FILTER)) { siteName <- "siteON" } else { siteName <- "siteOFF" }
if (isTRUE(disease.code_FILTER)) { dxName <- "histologicaldxON" } else { dxName <- "histologicaldxOFF" }
filterName <- paste(groupName,siteName, "_", dxName, "_", pathoName, "_",  ageName, sep="")

## Generate error file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

## Write output to file
#----------------------------------------------
sink(file = out.ouput, append = FALSE, split = FALSE)
options(max.print=999999)

## Print parameters to output
#----------------------------------------------
if (isTRUE(Internal_match & NCI_match)) {
  cat("Patient ID:", patient.id, "\n", "\n",
      "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
      "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
      "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
      "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else if (isTRUE(Internal_match)) {
  cat("Patient ID:", patient.id, "\n", "\n",
      "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
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
  
  sink(file = out.ouput, append = TRUE, split = FALSE)
  options(max.print=999999)
  
}

## PIPELINE
#----------------------------------------------
# Clean up patient data
source("Patient_Export_QC.R")

# Extract patient_id
patient.list <- sort(unique(STAMP_DF$PatientID))

if (length(patient.list) > 0) {
  
  if (isTRUE(Internal_match)) {
    source("Biomarker_Report_QC.R")
    source("Biomarker_Report_Match.R") 
  }
  
  if (isTRUE(NCI_match)) {
    source("Patient_Variant_Report_QC.R")
    source("Patient_Variant_Report_InclusionMatch.R")
    source("Patient_Variant_Report_NonHotspotMatch.R")
  }
  
  ## Generate OUTPUT file 
  source("ClinicalTrial_MatchOutput.R")
}

sink()

# Delete temporary directory
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)} 
