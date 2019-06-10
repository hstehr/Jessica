suppressMessages(library("easypackages"))
suppressMessages(libraries("plyr","dplyr","Biobase","eeptools","splitstackshape",
                           "reshape","rio","stringr","openxlsx","tidyr"))

rm(list=ls())

# Identify entries to map = test orders 
#----------------------------------------------
# Load file
out.output = "~/Desktop/output.testing.txt"
tempdir = "~/Desktop/"
source("~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TestOrder_QC.R")
DF <- TRF_DF

# Extract fields of interest
UniqCombo <- data.frame(primaryTumorSite = TRF_DF$PrimaryTumorSite,
                        histologicalDiagnosis = TRF_DF$HistologicalDx,
                        stringsAsFactors = FALSE)
# Remove duplicate entries
UniqCombo <- unique(UniqCombo[,])
# Remove empty rows
UniqCombo <- UniqCombo[rowSums(is.na(UniqCombo)) != ncol(UniqCombo),]  

# Identify entries to map = current mappings
#----------------------------------------------
# Load file
setwd("~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/")
histoDx.key="SyapseExport_RetrospectiveAnalysis/HistologicalDx_CTEP.csv"
source("HistologicalDx_CTEP_Match.R")

DF_Mapped <- HistologicalDxCategory
# Remove duplicate entries
DF_Mapped <- unique(DF_Mapped[,c("histologicalDiagnosis","primaryTumorSite")])
# Remove empty rows
DF_Mapped <- DF_Mapped[rowSums(is.na(DF_Mapped)) != ncol(DF_Mapped),]  

# Extract entries missing mappings 
#----------------------------------------------
Entries2Map <- anti_join(UniqCombo,DF_Mapped,
                         by = c("primaryTumorSite", "histologicalDiagnosis"))
Entries2Map <- Entries2Map[order(Entries2Map$histologicalDiagnosis),]
Entries2Map <- Entries2Map[order(Entries2Map$primaryTumorSite),]

# Remove duplicate entries
Entries2Map <- unique(Entries2Map[,])
# Remove empty rows
Entries2Map <- Entries2Map[rowSums(is.na(Entries2Map)) != ncol(Entries2Map),]  

# sort(unique(Entries2Map$primaryTumorSite))
# sort(unique(Entries2Map$histologicalDiagnosis))

write.table(Entries2Map, file = "~/Desktop/Entries2Map.tsv",
            append = FALSE, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

remove(out.output,tempdir,TRF_DF,DF,patient_list,TRF.file,
       DF_Mapped,HistologicalDxCategory,UniqCombo,histoDx.key)

# Identify disease exclusion entries
#----------------------------------------------
# Upload patient variant reports that compromise of different arms
# All ARMs consistent from (2018-09-06 until 2019-05-30, inclusively)

setwd("~/Documents/ClinicalDataScience_Fellowship/")
NCI.ArmRemove="NULL"

NCI.file="PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2018-09-06.xlsx"
#----------------------------------------------
PATIENT_VARIANT_REPORT <- suppressMessages(import_list(NCI.file, setclass = "tbl"))
Patient_Variant_Report_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1", sub(".*_", "", NCI.file))), format= "%Y-%m-%d")

colname_disease <- c("CTEP.CATEGORY","CTEP.SUBCATEGORY","CTEP.TERM","SHORT.NAME","MedDRA.CODE")
DF_Histologic_Disease_Exclusion_Codes <- data.frame(matrix(NA, ncol = length(colname_disease) +1))

for (tab_No in which(names(PATIENT_VARIANT_REPORT) == "Disease Exclusion LOOK-UP Table")) {
  Disease_Exclusion_file <- PATIENT_VARIANT_REPORT[[tab_No]]
  
  ## Remove rows that are all empty
  Disease_Exclusion_file <- Disease_Exclusion_file[rowSums(is.na(Disease_Exclusion_file)) != 
                                                     ncol(Disease_Exclusion_file),]  
  
  ## Identify Arms based on iteration of specific string
  row_start_list = (which(Disease_Exclusion_file[[1]] == "Histologic Disease Exclusion Codes"))
  
  for (Histo_No in 1:length(row_start_list)) {
    
    # Assume entries begin 2 lines after specific string
    row_start = row_start_list[Histo_No] +2
    
    # Extract ARM_Name
    if (isTRUE(row_start_list[Histo_No] == 1)) { 
      # First Arm: top header line is read as column name
      Arm_Name <- colnames(Disease_Exclusion_file)[2]
      
    } else { 
      # Other Arms: Arm name is top header line relative to specific string
      Arm_Name <- as.character(Disease_Exclusion_file[row_start_list[Histo_No] -1,2])
    }
    
    # Abbreviate Arm name
    Arm_Name <- gsub("^(EAY131)(.*)", "ARM\\2", Arm_Name)
    
    if (isTRUE(Histo_No == length(row_start_list))) { 
      # Last Arm: last line is last line of file 
      row_end = as.numeric(nrow(Disease_Exclusion_file))
      
    } else { 
      # Other Arms: last line is before start of header of next iteration
      row_end = row_start_list[Histo_No +1] -2 
    }
    
    # Extract rows per ARM_Name
    DF_Histologic_Disease_Exclusion_pre <- Disease_Exclusion_file[c(row_start:row_end), 1:length(colname_disease)]
    DF_Histologic_Disease_Exclusion_pre$Arm_Name <- Arm_Name
    
    if (isTRUE(row_start_list[Histo_No] == 1)) { 
      # First Arm: start of dataframe
      DF_Histologic_Disease_Exclusion_Codes <- DF_Histologic_Disease_Exclusion_pre
      
    } else { 
      # Other Arms: append to dataframe
      DF_Histologic_Disease_Exclusion_Codes <- rbind(DF_Histologic_Disease_Exclusion_Codes,
                                                     DF_Histologic_Disease_Exclusion_pre)
    }
  }
  remove(DF_Histologic_Disease_Exclusion_pre, Disease_Exclusion_file,Arm_Name,Histo_No,row_end,row_start,row_start_list,tab_No)
}

colnames(DF_Histologic_Disease_Exclusion_Codes) <- c(colname_disease,"Arm_Name")
DF_Histologic_Disease_Exclusion_Codes <- DF_Histologic_Disease_Exclusion_Codes[,c("Arm_Name",colname_disease)]
DF_Histologic_Disease_Exclusion_Codes <- unique(DF_Histologic_Disease_Exclusion_Codes[,c("CTEP.CATEGORY","CTEP.SUBCATEGORY","CTEP.TERM","SHORT.NAME","MedDRA.CODE")])

DF_Histologic_Disease_Exclusion_Codes <- DF_Histologic_Disease_Exclusion_Codes[!is.na(DF_Histologic_Disease_Exclusion_Codes$MedDRA.CODE),]

write.table(DF_Histologic_Disease_Exclusion_Codes, file = "~/Desktop/Disease_Exclusion_Codes.tsv",
            append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

remove(Patient_Variant_Report_timestamp,PATIENT_VARIANT_REPORT,colname_disease,NCI.ArmRemove,NCI.file)
