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

## Customize trial output
if (isTRUE(OnCore.file == "FALSE")) {Internal_match = FALSE  
} else {Internal_match = TRUE}

if (isTRUE(NCI.file == "FALSE")) {NCI_match = FALSE  
} else {NCI_match = TRUE}

## Directories
#----------------------------------------------
outdir = paste(outdir.root,"/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(outdir.root,"/temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.output = paste(outdir,Sys.Date(),".out",sep="")
err.output = paste(outdir,Sys.Date(),".err",sep="")

## Load files and specify timestamps
#----------------------------------------------
STAMP_DF <- read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")
Syapse_Export_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*/", "", STAMP.file))), format= "%Y-%m-%d")

if (isTRUE(Internal_match)) {
  OnCore_Biomarker_Report <- read.csv(file = OnCore.file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, sep = ",")
  OnCore_Biomarker_Report_timestamp <- format(as.Date(paste(gsub("([[:digit:]]{4}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*_", "", OnCore.file)), "-01",sep="")), format= "%Y-%m")
}

if (isTRUE(NCI_match)) {
  PATIENT_VARIANT_REPORT <- suppressMessages(import_list(NCI.file, setclass = "tbl"))
  Patient_Variant_Report_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1", sub(".*_", "", NCI.file))), format= "%Y-%m-%d")
}

stamp_reference_transcripts <- read.csv(file = stamp_reference_transcripts, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <- read.csv(file = exons_ensembl, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl, by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)

## Specify script location
setwd(script.root)

## Disease Ontology 
source("DiseaseGroupCategories.R")
source("HistologicalDx_CTEP_Match.R")
source("CustomPalettes.R")

## Default filters 
static.plots_FILTER = TRUE # Generate static visual graphics
adult.group_FILTER = as.logical(adult_FILTER) # Adult = age >= 18+yo)
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")
benign_removed <- c("Likely Benign")

## Customizeable filters
pathogenic_FILTER = as.logical(pathogenic_FILTER)
disease.group_FILTER = as.logical(disease_FILTER)
disease.site_FILTER = as.logical(disease_FILTER) # Dependent on disease.group_FILTER == TRUE
disease.code_FILTER = as.logical(disease_FILTER)

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

# ## Write output to file
# #----------------------------------------------
# sink(file = out.output, append = FALSE, split = FALSE)
# options(max.print=999999)

## Print parameters to output
#----------------------------------------------
if (isTRUE(Internal_match & NCI_match)) {
  cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
      "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
      "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
      "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
      "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else if (isTRUE(Internal_match)) {
  cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
      "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
      "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
      "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
      "Outdirectory: ", data.root, "\n",
      "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
      "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
  
} else if (isTRUE(NCI_match)) {
  cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
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
# Merge exports from timestamp = c("2018-10-18") and timestamp = c("2019-04-30") - Filtered for "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"
STAMP_DF_old <- read.csv(file ="~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-10-18_syapse_export_all_variants_patientNameAndMrnRemoved.csv", 
                         header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")
STAMP_DF_old <- 
  STAMP_DF_old[which(STAMP_DF_old$smpl.CancerSomaticMutationReport...smpl.assayName != "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"),]

# Incorporate additional columns
DF_TumorSite <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-01-31_syapse_export_all.csv",
                         header = TRUE, na.strings = c("NA","None"), stringsAsFactors = FALSE,sep = ",")
# Remove extraneous columns
DF_TumorSite <- DF_TumorSite[,c("UNIQUE_ID","PRIMARY_TUMOR_SITE","HISTOLOGICAL_DIAGNOSIS")]
colnames(DF_TumorSite) <- c("smpl.TestRequest...sys.uniqueId","smpl.Patient...smpl.primaryTumorSite","smpl.Patient...smpl.histologicalDiagnosis")

# Remove duplicate rows   
DF_TumorSite <- DF_TumorSite %>% dplyr::distinct(smpl.TestRequest...sys.uniqueId, .keep_all = TRUE)

# Merge with STAMP entries 
STAMP_DF_old <- left_join(STAMP_DF_old, DF_TumorSite, by = c("smpl.TestRequest...sys.uniqueId"))

# Format dob entries
STAMP_DF$smpl.CancerSomaticMutationReport...base.dob <- 
  format(as.Date(STAMP_DF$smpl.CancerSomaticMutationReport...base.dob, format = "%Y-%m-%d"), "%m/%d/%y")

# Merge datasets 
column_common <- intersect(colnames(STAMP_DF),colnames(STAMP_DF_old))

STAMP_DF <- rbind(STAMP_DF[,column_common],STAMP_DF_old[,column_common])
remove(STAMP_DF_old,column_common,DF_TumorSite)

# Clean up patient data from Syapse
#----------------------------------------------
source("Syapse_Export_QC.R") # Includes manual edits for Syapse export
source("Syapse_VariantAnnotate.R")

cat("STAMP - Solid Tumor Actionable Mutation Panels","\n")
cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

#----------------------------------------------
## # Replication of Appendix A MATCH Designated Lab Application: SNV/Indels
#----------------------------------------------
# Filter patients from 2017-07-01 to 2017-12-31, inclusive
STAMP_DF$AssayReportDateReviewed <- as.Date(STAMP_DF$AssayReportDateReviewed, format = "%Y-%m-%d")
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayReportDateReviewed >= "2017-07-01" &
                             STAMP_DF$AssayReportDateReviewed <= "2017-12-31"),]

# Filter for STAMP v2
#----------------------------------------------
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_DF$AssayName))

# Remove benign mutations 
#----------------------------------------------
STAMP_DF <- STAMP_DF[!(STAMP_DF$VariantPathogenicityStatus %in% benign_removed), ]
# sort(unique(STAMP_DF$VariantPathogenicityStatus))

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

DF = STAMP_DF
cohort = "all"
source("Syapse_Visualizations.R")

#----------------------------------------------
## Pipeline resumes
#----------------------------------------------
cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

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
  
  if (isTRUE(static.plots_FILTER)) {source("ClinicalTrial_Graphics.R")}
  # source("ClinicalTrial_Output_Details.R")
  # source("ClinicalTrial_Output_tsv.R")
  # source("ClinicalTrial_Output_Candidates.R")
  
}

closeAllConnections() 
