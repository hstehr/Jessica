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
outdir = paste(outdir.root,"/",sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
tempdir = paste(outdir.root,"/temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

# Specify output file
out.output = paste(outdir,Sys.Date(),".out",sep="")
err.output = paste(outdir,Sys.Date(),".err",sep="")

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
source("CustomPalettes.R")

## Default filters 
static.plots_FILTER = TRUE # Generate static visual graphics
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")
benign_removed <- c("Likely Benign")

## Customizeable filters
adult.group_FILTER = as.logical(adult_FILTER) # Adult = age >= 18+yo)

## Generate error file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

# ## Write output to file
# #----------------------------------------------
# sink(file = out.output, append = FALSE, split = FALSE)
# options(max.print=999999)

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
source("Syapse_Export_QC.R") # Includes manual edits for Syapse export
source("Syapse_VariantAnnotate.R")

# All STAMP version analyzed
#----------------------------------------------
cat("STAMP - Solid Tumor Actionable Mutation Panels","\n")
cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

# Remove benign mutations 
#----------------------------------------------
STAMP_DF <- STAMP_DF[!(STAMP_DF$VariantPathogenicityStatus %in% benign_removed), ]
# sort(unique(STAMP_DF$VariantPathogenicityStatus))

# Starting DF for each iteration
#----------------------------------------------
STAMP_DF_complete <- STAMP_DF

print(table(sort(STAMP_DF$AssayName)))

# Data visualizations
DF = STAMP_DF
cohort = "all"
source(paste(script.root,"Syapse_Visualizations.R",sep=""))

# FUNCTIONS
#----------------------------------------------
parameter_filter_fxn <- function() {
  ## Filters APPLIED
  if (isTRUE(adult.group_FILTER)) { ageName <- "AdultGroupON" } else { ageName <- "AdultGroupOFF" }
  if (isTRUE(pathogenic_FILTER)) { pathoName <- "pathogenicON" } else { pathoName <- "pathogenicOFF" }
  if (isTRUE(disease.group_FILTER)) { groupName <- "diseaseFILTER_groupON" } else { groupName <- "diseaseFILTER_groupOFF" }
  if (isTRUE(disease.group_FILTER & disease.site_FILTER)) { siteName <- "siteON" } else { siteName <- "siteOFF" }
  if (isTRUE(disease.code_FILTER)) { dxName <- "histologicaldxON" } else { dxName <- "histologicaldxOFF" }
  filterName <- paste(groupName, siteName, "_", dxName, "_", pathoName, "_",  ageName, sep="")
  
  assign("groupName", groupName, envir = .GlobalEnv)
  assign("siteName", siteName, envir = .GlobalEnv)
  assign("dxName", dxName, envir = .GlobalEnv)
  assign("pathoName", pathoName, envir = .GlobalEnv)
  assign("ageName", ageName, envir = .GlobalEnv)
  assign("filterName", filterName, envir = .GlobalEnv)
  
  ## Print parameters to output
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
}

Iterate_Fxn <- function() {
  patient.list <- sort(unique(STAMP_DF$PatientID))
  assign("patient.list", patient.list, envir = .GlobalEnv)
  
  if (length(patient.list) > 0) {
    
    if (isTRUE(Internal_match)) {
      source(paste(script.root,"Biomarker_Report_QC.R",sep=""))
      source(paste(script.root,"Biomarker_Report_Match.R",sep="")) 
    }
    
    if (isTRUE(NCI_match)) {
      source(paste(script.root,"Patient_Variant_Report_QC.R",sep=""))
      source(paste(script.root,"Patient_Variant_Report_InclusionMatch.R",sep=""))
      source(paste(script.root,"Patient_Variant_Report_NonHotspotMatch.R",sep=""))
    }
    
    if (isTRUE(static.plots_FILTER)) {source(paste(script.root,"ClinicalTrial_Graphics.R",sep=""))}
    # source(paste(script.root,"ClinicalTrial_Output_Details.R",sep=""))
    # source(paste(script.root,"ClinicalTrial_Output_tsv.R",sep=""))
    # source(paste(script.root,"ClinicalTrial_Output_Candidates.R",sep=""))
  }
}
  
closeAllConnections()

#----------------------------------------------
## Algorithm testing on retrospective SNV/Indel data
#----------------------------------------------
STAMP_DF <- STAMP_DF_complete

outdir = "~/Desktop/trials_Iterate/FullDataset_AllFilters"
if (!dir.exists(outdir)){dir.create(outdir)} 
setwd(outdir)

out.output.sub = paste(outdir,"/",Sys.Date(),".out",sep="")
sink(file = out.output.sub, append = FALSE, split = FALSE)
options(max.print=999999)

cat(paste("Outdirectory: ", outdir, sep=""),"\n")

pathogenic_FILTER = as.logical("TRUE")
disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = as.logical("TRUE")

parameter_filter_fxn ()

# Filter for adults
#----------------------------------------------
if (isTRUE(adult.group_FILTER)) {STAMP_DF <- STAMP_DF[which(STAMP_DF$PatientAge >= 18),]}

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

Iterate_Fxn()
  
closeAllConnections() 

#----------------------------------------------
## Algorithm testing on retrospective SNV/Indel data
#----------------------------------------------
STAMP_DF <- STAMP_DF_complete

# Filter for STAMP v2
#----------------------------------------------
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_DF$AssayName))

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# sort(unique(STAMP_DF$HistologicalDx))

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$PrimaryTumorSite),]
STAMP_DF <- STAMP_DF[!(STAMP_DF$PrimaryTumorSite %in% c("unknown","none","other primary site")), ]

# Collapse similar primary tumor site
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("testes","testis"))] <- "testis"
# sort(unique(STAMP_DF$PrimaryTumorSite))

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

# Data visualizations
outdir = paste(outdir.root,"/",sep="")
DF = STAMP_DF

# Abbreviate for visualization
DF$PrimaryTumorSite[which(DF$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
DF$PrimaryTumorSite[which(DF$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"

cohort = "all_QC"
source(paste(script.root,"Syapse_Visualizations.R",sep=""))

# Filter for adults
#----------------------------------------------
if (isTRUE(adult.group_FILTER)) {STAMP_DF <- STAMP_DF[which(STAMP_DF$PatientAge >= 18),]}

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

# ITERATE through each case
#----------------------------------------------
outdir.list <- c("CompleteDx_NoFilter","CompleteDx_PathoON","CompleteDx_DiseaseON","CompleteDx_AllFilter")
pathogenic_FILTER.list <- c("FALSE","TRUE","FALSE","TRUE")
disease_FILTER.list <- c("FALSE","FALSE","TRUE","TRUE")

for (case_No in 1:length(outdir.list)) {
  outdir = paste("~/Desktop/trials_Iterate/", outdir.list[case_No], sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  setwd(outdir)
  
  out.output.sub = paste(outdir,"/",Sys.Date(),".out",sep="")
  sink(file = out.output.sub, append = FALSE, split = FALSE)
  options(max.print=999999)
  
  cat(paste("Outdirectory: ", outdir, sep=""),"\n")
  cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")
  
  pathogenic_FILTER = as.logical(pathogenic_FILTER.list[case_No])
  disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = as.logical(disease_FILTER.list[case_No])
  
  parameter_filter_fxn ()
  Iterate_Fxn()

  closeAllConnections() 
  
  ## Extract fields of Case No. 4 for presentation
  #----------------------------------------------
  if (isTRUE(case_No == "4")) {
    OnCoreMatch <- read.csv(file=paste(tempdir,"OnCore_SNVIndel_Matched_2019-03_diseaseFILTER_groupONsiteON_pathogenicON_AdultGroupON.tsv",sep=""),sep="\t")
    
    # Breakdown by trial & gene
    DF_tabulate <- data.frame(OnCoreMatch %>% group_by(OnCore.No,VariantGene) %>% tally())
    colnames(DF_tabulate) <- c("OnCore.No","VariantGene","No.Orders")
    DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders,decreasing = TRUE),]
    
    write.table(DF_tabulate, file = paste(outdir, "/DF_tabulate.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    View(DF_tabulate)
    
    # Breakdown by trial 
    DF_tabulate <- data.frame(OnCoreMatch %>% group_by(OnCore.No) %>% tally())
    colnames(DF_tabulate) <- c("OnCore.No","No.Orders")
    DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders,decreasing = TRUE),]
    
    write.table(DF_tabulate, file = paste(outdir, "/DF_tabulate_trial.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    View(DF_tabulate)
    
    remove(OnCoreMatch,DF_tabulate)  
  }
  remove(case_No)
}

#----------------------------------------------
## Algorithm testing on retrospective SNV/Indel data - remove lung cases
#----------------------------------------------
# Filter lung primary tumor site
STAMP_DF <- STAMP_DF[which(tolower(STAMP_DF$PrimaryTumorSite) != "lung"),]

cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")

# ITERATE through each case
#----------------------------------------------
for (case_No in 1:length(outdir.list)) {
  outdir = paste("~/Desktop/trials_Iterate/", outdir.list[case_No],"_noLung", sep="")
  if (!dir.exists(outdir)){dir.create(outdir)} 
  setwd(outdir)
  
  out.output.sub = paste(outdir,"/",Sys.Date(),".out",sep="")
  sink(file = out.output.sub, append = FALSE, split = FALSE)
  options(max.print=999999)
  
  cat(paste("Outdirectory: ", outdir, sep=""),"\n")
  cat(paste(nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total patients", sep=""),"\n","\n")
  
  pathogenic_FILTER = as.logical(pathogenic_FILTER.list[case_No])
  disease.group_FILTER = disease.site_FILTER = disease.code_FILTER = as.logical(disease_FILTER.list[case_No])
  
  parameter_filter_fxn ()
  Iterate_Fxn()
  
  closeAllConnections() 
  
  ## Extract fields of Case No. 4 for presentation
  #----------------------------------------------
  if (isTRUE(case_No == "4")) {
    OnCoreMatch <- read.csv(file=paste(tempdir,"OnCore_SNVIndel_Matched_2019-03_diseaseFILTER_groupONsiteON_pathogenicON_AdultGroupON.tsv",sep=""),sep="\t")
    
    # Breakdown by trial & gene
    DF_tabulate <- data.frame(OnCoreMatch %>% group_by(OnCore.No,VariantGene) %>% tally())
    colnames(DF_tabulate) <- c("OnCore.No","VariantGene","No.Orders")
    DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders,decreasing = TRUE),]
    
    write.table(DF_tabulate, file = paste(outdir, "/DF_tabulate.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    View(DF_tabulate)
    
    # Breakdown by trial 
    DF_tabulate <- data.frame(OnCoreMatch %>% group_by(OnCore.No) %>% tally())
    colnames(DF_tabulate) <- c("OnCore.No","No.Orders")
    DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders,decreasing = TRUE),]
    
    write.table(DF_tabulate, file = paste(outdir, "/DF_tabulate_trial.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    View(DF_tabulate)
    
    remove(OnCoreMatch,DF_tabulate)
  }
}

## Remove tempdir
#----------------------------------------------
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)}
