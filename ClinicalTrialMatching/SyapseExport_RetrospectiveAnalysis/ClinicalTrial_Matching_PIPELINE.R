suppressMessages(library("easypackages"))
suppressMessages(libraries("plyr","dplyr","Biobase","eeptools","splitstackshape",
                           "reshape","rio","stringr","openxlsx","tidyr"))

## Customize trial output
if (isTRUE(OnCore.file == "FALSE")) {Internal_match = FALSE  
} else {Internal_match = TRUE}

if (isTRUE(NCI.file == "FALSE")) {NCI_match = FALSE  
} else {NCI_match = TRUE}

## Directories
#----------------------------------------------
outdir = outdir.root
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(outdir.root,"temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.output = paste(outdir,Sys.Date(),".out.txt",sep="")
err.output = paste(outdir,Sys.Date(),".err.txt",sep="")

## Load files and specify timestamps
#----------------------------------------------
AA_key_table <- read.csv(file = AA_key, sep = ",")

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

# Functions
source("SyapseExport_RetrospectiveAnalysis/ClinicalTrial_Matching_Functions.R")
source(paste(data.root,"STAMPEDE_Visualizations/STAMPEDE_Functions.R",sep=""))

## Default filters 
saveStaticPlots = TRUE
saveDynamicPlots = FALSE
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

## Customizeable filters
adult.group_FILTER = as.logical(adult_FILTER) # Adult = age >= 18+yo)

## Generate error file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

## Write output to file
#----------------------------------------------
sink(file = out.output, append = FALSE, split = FALSE)
options(max.print=999999)

#----------------------------------------------
## PIPELINE: SNV Indel
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
source("SyapseExport_RetrospectiveAnalysis/Syapse_Export_QC.R") # Includes manual edits for Syapse export
source("SyapseExport_RetrospectiveAnalysis/Syapse_VariantAnnotate.R")

#----------------------------------------------
## PIPELINE: Fusion = STAMP_Fusion
#----------------------------------------------
source("SyapseExport_RetrospectiveAnalysis/Fusion_Export_QC.R")

#----------------------------------------------
## PIPELINE: CNV
#----------------------------------------------
source("SyapseExport_RetrospectiveAnalysis/STAMP_TestOrder_QC.R")
source("SyapseExport_RetrospectiveAnalysis/CNV_Export_QC.R")

#################################
# Starting QC-filtered (general) data sources: STAMP v2
#################################
assay_select = "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"
STAMP_DF_QC <- STAMP_DF[which(STAMP_DF$AssayName == assay_select),]
STAMP_Fusion_QC <- STAMP_Fusion[which(STAMP_Fusion$AssayName == assay_select),]
STAMP_CNV_QC <- STAMP_CNV[which(STAMP_CNV$AssayName == assay_select),]
remove(assay_select)

# Apply abbreviated primary tumor sites
shorthand_visual_fxn (DF = STAMP_DF_QC)
STAMP_DF_plot <- DF
shorthand_visual_fxn (DF = STAMP_Fusion_QC)
STAMP_Fusion_plot <- unique(DF[,!names(DF) %in% "Gene"])
shorthand_visual_fxn (DF = STAMP_CNV_QC)
STAMP_CNV_plot <- DF

# Data visualizations
cohort_id = "STAMPv2"
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Visualizations.R",sep=""))
remove(DF,STAMP_DF_plot,STAMP_Fusion_plot,STAMP_CNV_plot,cohort_id)

#################################
## Exclude missing primary tumor site information
#################################
# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# sort(unique(STAMP_DF$HistologicalDx))

# Filter for entries with primary tumor site
assay_select = c("unknown","none","other primary site")

STAMP_DF_complete <- STAMP_DF_QC[!is.na(STAMP_DF_QC$PrimaryTumorSite),]
STAMP_DF_complete <- STAMP_DF_complete[!(STAMP_DF_complete$PrimaryTumorSite %in% assay_select), ]

STAMP_Fusion_complete <- STAMP_Fusion_QC[!is.na(STAMP_Fusion_QC$PrimaryTumorSite),]
STAMP_Fusion_complete <- STAMP_Fusion_complete[!(STAMP_Fusion_complete$PrimaryTumorSite %in% assay_select), ]

STAMP_CNV_complete <- STAMP_CNV_QC[!is.na(STAMP_CNV_QC$PrimaryTumorSite),]
STAMP_CNV_complete <- STAMP_CNV_complete[!(STAMP_CNV_complete$PrimaryTumorSite %in% assay_select), ]

remove(assay_select)

# Apply abbreviated primary tumor sites
shorthand_visual_fxn (DF = STAMP_DF_complete)
STAMP_DF_plot <- DF
shorthand_visual_fxn (DF = STAMP_Fusion_complete)
STAMP_Fusion_plot <- unique(DF[,!names(DF) %in% "Gene"])
shorthand_visual_fxn (DF = STAMP_CNV_complete)
STAMP_CNV_plot <- DF

# Data visualizations
cohort_id = "STAMPv2complete"
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Visualizations.R",sep=""))
remove(DF,STAMP_DF_plot,STAMP_Fusion_plot,STAMP_CNV_plot,cohort_id)

#################################
## AgeFILTER = Filter for adults
#################################
STAMP_DF_age <- STAMP_DF_complete[which(STAMP_DF_complete$PatientAge >= 18),]
STAMP_Fusion_age <- STAMP_Fusion_complete[which(STAMP_Fusion_complete$PatientAge >= 18),]
STAMP_CNV_age <- STAMP_CNV_complete[which(STAMP_CNV_complete$PatientAge >= 18),]

# Apply abbreviated primary tumor sites
shorthand_visual_fxn (DF = STAMP_DF_age)
STAMP_DF_plot <- DF
shorthand_visual_fxn (DF = STAMP_Fusion_age)
STAMP_Fusion_plot <- unique(DF[,!names(DF) %in% "Gene"])
shorthand_visual_fxn (DF = STAMP_CNV_age)
STAMP_CNV_plot <- DF

# Data visualizations
cohort_id = "STAMPv2ageFILTER"
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Visualizations.R",sep=""))
remove(DF,STAMP_DF_plot,STAMP_Fusion_plot,STAMP_CNV_plot,cohort_id)

closeAllConnections()

#----------------------------------------------
# ITERATE through each case
#----------------------------------------------
outdir.list <- c("CompleteDx_NoFilter/","CompleteDx_PathoON/","CompleteDx_DiseaseON/","CompleteDx_AllFilter/")
pathogenic_FILTER.list <- c("FALSE","TRUE","FALSE","TRUE")
disease_FILTER.list <- c("FALSE","FALSE","TRUE","TRUE")

# Specify dataframes to iterate
STAMP_DF_iterate = STAMP_DF_age
STAMP_Fusion_iterate = STAMP_Fusion_age
STAMP_CNV_iterate = STAMP_CNV_age

for (case_No in 1:length(outdir.list)) {
  outdir = paste(outdir.root, outdir.list[case_No], sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  setwd(outdir)
  
  out.output.sub = paste(outdir,Sys.Date(),".out.txt",sep="")
  sink(file = out.output.sub, append = FALSE, split = FALSE)
  options(max.print=999999)
  
  cat(paste("Outdirectory: ", outdir, sep=""),"\n")
  pathogenic_FILTER = as.logical(pathogenic_FILTER.list[case_No])
  disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = as.logical(disease_FILTER.list[case_No])
  
  parameter_filter_fxn ()
  
  outdir = paste(outdir,"Match_Results/",sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  
  Iterate_Fxn(STAMP_DF_iterate = STAMP_DF_age,
              STAMP_Fusion_iterate = STAMP_Fusion_age,
              STAMP_CNV_iterate = STAMP_CNV_age)
  
  outdir = paste(outdir.root, outdir.list[case_No], sep="")
  ## Generate match rate visualizations
  patient.list=patient.list.full
  source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/ClinicalTrial_MatchRate_Graphics.R",sep=""))
  
  closeAllConnections() 
  
  ## Extract fields of Case No. 4 for presentation
  if (isTRUE(case_No == "4")) {OnCore_Tally_fxn ()}
  
  remove(case_No)
}
remove(STAMP_DF_iterate,STAMP_CNV_iterate)

#################################
## Exclude lung cases
#################################
sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

STAMP_DF_noLung <- STAMP_DF_age[which(tolower(STAMP_DF_age$PrimaryTumorSite) != "lung"),]
STAMP_Fusion_noLung <- STAMP_Fusion_age[which(tolower(STAMP_Fusion_age$PrimaryTumorSite) != "lung"),]
STAMP_CNV_noLung <- STAMP_CNV_age[which(tolower(STAMP_CNV_age$PrimaryTumorSite) != "lung"),]

# Apply abbreviated primary tumor sites
shorthand_visual_fxn (DF = STAMP_DF_noLung)
STAMP_DF_plot <- DF
shorthand_visual_fxn (DF = STAMP_Fusion_noLung)
STAMP_Fusion_plot <- unique(DF[,!names(DF) %in% "Gene"])
shorthand_visual_fxn (DF = STAMP_CNV_noLung)
STAMP_CNV_plot <- DF

# Data visualizations
cohort_id = "STAMPv2noLung"
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Visualizations.R",sep=""))
remove(DF,STAMP_DF_plot,STAMP_Fusion_plot,STAMP_CNV_plot,cohort_id)

closeAllConnections()

# ITERATE through each case
#----------------------------------------------
# Specify dataframes to iterate
STAMP_DF_iterate = STAMP_DF_noLung
STAMP_Fusion_iterate = STAMP_Fusion_noLung
STAMP_CNV_iterate = STAMP_CNV_noLung

for (case_No in 1:length(outdir.list)) {
  outdir = paste(outdir.root, gsub("/$","",outdir.list[case_No]),"_noLung/", sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  setwd(outdir)
  
  out.output.sub = paste(outdir,Sys.Date(),".out",sep="")
  sink(file = out.output.sub, append = FALSE, split = FALSE)
  options(max.print=999999)
  
  cat(paste("Outdirectory: ", outdir, sep=""),"\n")
  pathogenic_FILTER = as.logical(pathogenic_FILTER.list[case_No])
  disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = as.logical(disease_FILTER.list[case_No])
  
  parameter_filter_fxn ()
  
  outdir = paste(outdir,"Match_Results/",sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  
  Iterate_Fxn(STAMP_DF_iterate = STAMP_DF_noLung,
              STAMP_Fusion_iterate = STAMP_Fusion_noLung,
              STAMP_CNV_iterate = STAMP_CNV_noLung)
  
  outdir = paste(outdir.root, gsub("/$","",outdir.list[case_No]),"_noLung/", sep="")
  ## Generate match rate visualizations
  patient.list=patient.list.full
  source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/ClinicalTrial_MatchRate_Graphics.R",sep=""))
  
  closeAllConnections() 
  
  ## Extract fields of Case No. 4 for presentation
  if (isTRUE(case_No == "4")) {OnCore_Tally_fxn ()}
  
  remove(case_No)
}

closeAllConnections() 

#----------------------------------------------
## Replication of patients identified in NCI-MATCH Designated Lab Application: SNV/Indels
#----------------------------------------------
## Need to input "Patient ID information" and MOIs" tabs for parse ARMS
## Switch REF and ALT columns**
## Rename "gene" > "Gene Name" for "Inclusion Non-Hotspot Rules" in ARM-H, ARM-U
## Modifications to HistologicalDx_CTEP.csv: append "adenocarcinoma,lung,,,Non-small cell lung cancer,"

## Directories
#----------------------------------------------
outdir = paste(outdir.root,"LabApp_Replication/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 

# Specify output file
out.output = paste(outdir,Sys.Date(),".out.txt",sep="")
err.output = paste(outdir,Sys.Date(),".err.txt",sep="")

## Generate error file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

## Write output to file
#----------------------------------------------
sink(file = out.output, append = FALSE, split = FALSE)
options(max.print=999999)

## Specify parameters 
#----------------------------------------------
if (isTRUE(OnCore.file_LabRep == "FALSE")) {Internal_match = FALSE  
} else {Internal_match = TRUE}

if (isTRUE(NCI.file_LabRep == "FALSE")) {NCI_match = FALSE  
} else {NCI_match = TRUE}
NCI.ArmRemove <- NCI.ArmRemove_LabRep

if (isTRUE(NCI_match)) {
  PATIENT_VARIANT_REPORT <- suppressMessages(import_list(NCI.file_LabRep, setclass = "tbl"))
  Patient_Variant_Report_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1", sub(".*_", "", NCI.file_LabRep))), format= "%Y-%m-%d")
}

# Extract entries of interest
#----------------------------------------------
# Filter patients from 2017-07-01 to 2017-12-31, inclusive
start_date = "2017-07-01"
end_date = "2017-12-31"

STAMP_DF_QC$AssayReportDateReviewed <- as.Date(STAMP_DF_QC$AssayReportDateReviewed, format = "%Y-%m-%d")
STAMP_Fusion_QC$AssayReportDateReviewed <- as.Date(STAMP_Fusion_QC$AssayReportDateReviewed, format = "%Y-%m-%d")
STAMP_CNV_QC$AssayReportDateReviewed <- as.Date(STAMP_CNV_QC$AssayReportDateReviewed, format = "%Y-%m-%d")

STAMP_DF_LabRep <- STAMP_DF_QC[which(STAMP_DF_QC$AssayReportDateReviewed >= start_date &
                                       STAMP_DF_QC$AssayReportDateReviewed <= end_date),]
STAMP_Fusion_LabRep <- STAMP_Fusion_QC[which(STAMP_Fusion_QC$AssayReportDateReviewed >= start_date &
                                               STAMP_Fusion_QC$AssayReportDateReviewed <= end_date),]
STAMP_CNV_LabRep <- STAMP_CNV_QC[which(STAMP_CNV_QC$AssayReportDateReviewed >= start_date &
                                         STAMP_CNV_QC$AssayReportDateReviewed <= end_date),]

# Apply abbreviated primary tumor sites
shorthand_visual_fxn (DF = STAMP_DF_LabRep)
STAMP_DF_plot <- DF
shorthand_visual_fxn (DF = STAMP_Fusion_LabRep)
STAMP_Fusion_plot <- unique(DF[,!names(DF) %in% "Gene"])
shorthand_visual_fxn (DF = STAMP_CNV_LabRep)
STAMP_CNV_plot <- DF

# Data visualizations
cohort_id = "STAMPv2LabRep"
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Visualizations.R",sep=""))
remove(DF,STAMP_DF_plot,STAMP_Fusion_plot,STAMP_CNV_plot,cohort_id)

closeAllConnections()

# ITERATE through each case
#----------------------------------------------
# Specify dataframes to iterate
STAMP_DF_iterate = STAMP_DF_LabRep
STAMP_Fusion_iterate = STAMP_Fusion_LabRep
STAMP_CNV_iterate = STAMP_CNV_LabRep

sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

cat(paste("Outdirectory: ", outdir, sep=""),"\n")
adult.group_FILTER = FALSE
pathogenic_FILTER = TRUE
disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = TRUE

parameter_filter_fxn ()

outdir = paste(outdir,"Match_Results/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 

Iterate_Fxn(STAMP_DF_iterate = STAMP_DF_LabRep,
            STAMP_Fusion_iterate = STAMP_Fusion_LabRep,
            STAMP_CNV_iterate = STAMP_CNV_LabRep)

outdir = paste(outdir.root,"LabApp_Replication/",sep="")
## Generate match rate visualizations
patient.list=patient.list.full
source(paste(script.root,"SyapseExport_RetrospectiveAnalysis/ClinicalTrial_MatchRate_Graphics.R",sep=""))

LabRep_Summary_fxn ()

closeAllConnections() 

## Remove tempdir
#----------------------------------------------
# if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)}
