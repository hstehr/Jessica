#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(Biobase)
library(eeptools)
library(splitstackshape)
library(reshape)
library(rio)
library(stringr)
library(openxlsx)

args = commandArgs(trailingOnly=TRUE)
# 1. Directory to save output to.
outdir = args[1]
# 2. Annotation of Output folder.
outdir_anno = args[2]
# 3. Location of STAMP entries (syapse export).
STAMP.file = args[3]
# 4. Location of OnCore Report (Stanford Internal Clinical Trials).
OnCore.file = args[4]
# 5. Location of Patient Variant Report (NCI-MATCH Clinical Trials).
NCI.file = args[5]

## Load files 
#----------------------------------------------
STAMP_DF <- 
  read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")

OnCore_Biomarker_Report <- 
  read.csv(file = OnCore.file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, sep = ",")

PATIENT_VARIANT_REPORT <- import_list(NCI.file, setclass = "tbl")

## Specify Parameters 
#----------------------------------------------
setwd(outdir)

## Disease Ontology 
source("DiseaseGroupCategories.R")

## Timestamp
Syapse_Export_timestamp <- 
  format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",
                      sub(".*/", "", STAMP.file))), format= "%Y-%m-%d")

OnCore_Biomarker_Report_timestamp <- 
  format(as.Date(paste(gsub("([[:digit:]]{4}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*_", "", OnCore.file)), 
                       "-01",sep="")), format= "%Y-%m")

Patient_Variant_Report_timestamp <- 
  format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",
                      sub(".*_", "", NCI.file))), format= "%Y-%m-%d")

## Default filters 
adult.group_FILTER = TRUE # Adult = age >= 18+yo)
pathogenic_FILTER = TRUE
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

## Customizeable filters
disease.group_FILTER = TRUE
disease.site_FILTER = TRUE # Dependent on disease.group_FILTER == TRUE

## Filters APPLIED
if (isTRUE(disease.group_FILTER)) { groupName <- "diseaseFILTER_groupON" } else { groupName <- "diseaseFILTER_groupOFF" }
if (isTRUE(disease.group_FILTER & disease.site_FILTER)) { siteName <- "siteON" } else { siteName <- "siteOFF" }
if (isTRUE(pathogenic_FILTER)) { pathoName <- "pathogenicON" } else { pathoName <- "pathogenicOFF" }
if (isTRUE(adult.group_FILTER)) { ageName <- "AdultGroupON" } else { ageName <- "AdultGroupOFF" }
filterName <- paste(groupName,siteName, "_", pathoName, "_",  ageName, sep="")

## Directories
outdir_int = paste(getwd(), "/Temporary/", sep="")

outdir_patient = paste(getwd(), "/Results/", Sys.Date(), "_", outdir_anno, sep="")
if (!dir.exists(outdir_patient)){dir.create(outdir_patient)} 

## Write output to file
#----------------------------------------------
sink(file = paste(getwd(), "/Results/", Sys.Date(), "_",outdir_anno,"_Output.txt", sep=""), 
     append = FALSE, split = FALSE)
options(max.print=999999)

## Print parameters to output
#----------------------------------------------
cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
    "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
    "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
    "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
    "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
    "Outdirectory: ", getwd(), "\n",
    "Temporary Files within outdir: ", gsub(paste(getwd(),"/",sep=""), "", outdir_int), "\n",
    "Patient Match Results within outdir: ", gsub(paste(getwd(),"/",sep=""), "", outdir_patient), "\n","\n")

## PIPELINE
#----------------------------------------------
# Clean up patient data from Syapse
source("Syapse_Export_QC.R")
source("Syapse_VariantAnnotate.R")

# Extract patient_id
patient.list <- sort(unique(STAMP_DF$PatientID))

# Clean up Clinical Trial data
source("Biomarker_Report_QC.R")
source("Patient_Variant_Report_QC.R")

# Clinical trial Match
source("Biomarker_Report_Match.R")
source("Patient_Variant_Report_InclusionMatch.R")
source("Patient_Variant_Report_NonHotspotMatch.R")

## Generate OUTPUT file
source("ClinicalTrial_MatchOutput.R")

sink()
