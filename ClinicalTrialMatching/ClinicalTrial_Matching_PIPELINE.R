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
data.root = args[1]
# 2. Annotation of Output folder.
outdir_anno = args[2]
# 3. Location of STAMP entries (syapse export).
STAMP.file = args[3]
# 4. Location of OnCore Report (Stanford Internal Clinical Trials).
OnCore.file = args[4]
# 5. Location of Patient Variant Report (NCI-MATCH Clinical Trials).
NCI.file = args[5]
# 6. Directory of pipeline scripts.
script.root = args[6]
# 7. Location of stamp_reference_transcripts.
stamp_reference_transcripts = args[7]
# 8. Location of most recent exons_ensembl.
exons_ensembl = args[8]
# 9. Location of amino acid conversion key.
aminoAcid_conversion = args[9]
  
## Directories
outdir = paste(data.root,"trials/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(data.root,"temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.ouput = paste(outdir,Sys.Date(),"_",outdir_anno,".out",sep="")
err.output = paste(outdir,Sys.Date(),"_",outdir_anno,".err",sep="")

## Load files 
#----------------------------------------------
STAMP_DF <- 
  read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")

OnCore_Biomarker_Report <- 
  read.csv(file = OnCore.file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, sep = ",")

PATIENT_VARIANT_REPORT <- import_list(NCI.file, setclass = "tbl")

HistologicalDxCategory <- 
  read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-0223_HistologicalDx_CTEP.csv", 
           header = TRUE, na.strings = c(""," ","NA","."), stringsAsFactors = FALSE, sep = ",")

AminoAcid_Conversion <- 
  read.csv(file = aminoAcid_conversion, header = TRUE, na.strings = c(""," "), stringsAsFactors = FALSE, sep = ",")

stamp_reference_transcripts <- 
  read.csv(file = stamp_reference_transcripts, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <-
  read.csv(file = exons_ensembl, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl,
                                by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)
                                
## Specify script location
setwd(script.root)

## Disease Ontology 
source("DiseaseGroupCategories.R")

## Histological Dx Ontology
source("HistologicalDx_CTEP_Match.R")

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

## Customize trial output
Internal_match = FALSE
NCI_match = TRUE

## Default filters 
adult.group_FILTER = FALSE # Adult = age >= 18+yo)
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

## Customizeable filters
pathogenic_FILTER = FALSE
disease.group_FILTER = TRUE
disease.site_FILTER = TRUE # Dependent on disease.group_FILTER == TRUE
disease.code_FILTER = TRUE 

AgePlot.FILTER = TRUE
VariantPlot.FILTER = TRUE

## Filters APPLIED
if (isTRUE(adult.group_FILTER)) { ageName <- "AdultGroupON" } else { ageName <- "AdultGroupOFF" }
if (isTRUE(pathogenic_FILTER)) { pathoName <- "pathogenicON" } else { pathoName <- "pathogenicOFF" }
if (isTRUE(disease.group_FILTER)) { groupName <- "diseaseFILTER_groupON" } else { groupName <- "diseaseFILTER_groupOFF" }
if (isTRUE(disease.group_FILTER & disease.site_FILTER)) { siteName <- "siteON" } else { siteName <- "siteOFF" }
if (isTRUE(disease.code_FILTER)) { dxName <- "histologicaldxON" } else { dxName <- "histologicaldxOFF" }
filterName <- paste(groupName,siteName, "_", dxName, "_", pathoName, "_",  ageName, sep="")

## Write output to file
#----------------------------------------------
sink(file = out.ouput, append = FALSE, split = FALSE)
options(max.print=999999)

## Print parameters to output
#----------------------------------------------
cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
    "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
    "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
    "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
    "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
    "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
    "Outdirectory: ", data.root, "\n",
    "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
    "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")

## PIPELINE
#----------------------------------------------
# Clean up patient data from Syapse
source("Syapse_Export_QC.R")
source("Syapse_VariantAnnotate.R")

# Replication of Appendix A MATCH Designated Lab Application
# Filter patients from 2017-07-01 to 2017-12-31, inclusive
STAMP_DF$AssayReportDateReviewed <- as.Date(STAMP_DF$AssayReportDateReviewed, format = "%m/%d/%y")
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayReportDateReviewed >= "2017-07-01" &
                       STAMP_DF$AssayReportDateReviewed <= "2017-12-31"),]
write.table(STAMP_DF, file = paste(tempdir, Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

DF = STAMP_DF
cohort = "all"
source("Syapse_Visualizations.R")

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
source("ClinicalTrial_Graphics.R")
source("ClinicalTrial_MatchOutput.R")

sink()

# Delete temporary directory
# if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)} 

DF_Output_Patient_Variant_export <- 
  DF_Output_Patient_Variant[,c("Arm_Name","PatientID","PatientDOB","AssayReportDateReviewed",
                               "HistologicalDx","PrimaryTumorSite","VariantLabel","VariantPathogenicityStatus",
                               "Variant_Type","Protein")]

DF_Output_Patient_Variant_export <- 
  DF_Output_Patient_Variant_export[order(DF_Output_Patient_Variant_export$PatientDOB, decreasing = FALSE),]

DF_Output_Patient_Variant_export <- 
  DF_Output_Patient_Variant_export[order(DF_Output_Patient_Variant_export$Arm_Name, decreasing = FALSE),]

write.table(DF_Output_Patient_Variant_export, file = paste(outdir, "Output_Patient_Variant.csv", sep=""),
            append = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

remove(DF_Output_Patient_Variant_export)
