sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP database QC PIPELINE: use to map test order volume
## Remove duplicates of entries based on multiple unique sys.date_created
## Structure patient DOB and calculate current age
#################################

## List of test requests from STAMP v2 (130 genes) regardless of mutations identified.
TRF.file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-05-03_Syapse_export_TestRequests.csv"
TRF_DF <- read.csv(file = TRF.file, header = TRUE, na.strings = c(""," ","<NA>","NA","None"), stringsAsFactors = FALSE, sep = ",")

# Shorten column names of DF
colnames(TRF_DF) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(TRF_DF))

## Remove duplicates of entries based on multiple unique sys.date_created
## Include cases of smpl.amendedString == "AMENDED"
#----------------------------------------------
TRF_DF_AMENDED <- data.frame()
patient_list <- unique(TRF_DF$sys.uniqueId)

for (id_num in 1:length(patient_list)) {
  patient_id <- patient_list[id_num]
  
  assay_list <- unique(TRF_DF$smpl.assayName[which(TRF_DF$sys.uniqueId == patient_id)])
  
  for (assay_num in 1:length(assay_list)) {
    assay_id <- assay_list[assay_num]
    
    if (length(unique(TRF_DF$sys.date_created[which(TRF_DF$sys.uniqueId == patient_id & TRF_DF$smpl.assayName == assay_id)])) > 1) {
      created_list <- unique(TRF_DF$sys.date_created[which(TRF_DF$sys.uniqueId == patient_id & TRF_DF$smpl.assayName == assay_id)])
      
      ## Retain most recent version if multiple sys.date_created.1 exist in report 
      TRF_DF_AMENDED_pre <- TRF_DF[which(TRF_DF$sys.uniqueId == patient_id & TRF_DF$smpl.assayName == assay_id),]
      TRF_DF_AMENDED_pre <- TRF_DF_AMENDED_pre[TRF_DF_AMENDED_pre$sys.date_created == max(created_list),]
      
      remove(created_list)
      
    } else {
      TRF_DF_AMENDED_pre <- TRF_DF[which(TRF_DF$sys.uniqueId == patient_id & TRF_DF$smpl.assayName == assay_id),]
    }
    
    TRF_DF_AMENDED <- rbind(TRF_DF_AMENDED, TRF_DF_AMENDED_pre)
    remove(assay_id,assay_num)
  }
  remove(patient_id,assay_list,id_num)
}

# Overwrite existing dataframe
TRF_DF <- TRF_DF_AMENDED
remove(TRF_DF_AMENDED,TRF_DF_AMENDED_pre)

# Use sys.date_created as alternative for missing smpl.dateReceived
date_NA <- which(is.na(TRF_DF$smpl.dateReceived))
for (i in 1:length(date_NA)) {
  row_id = date_NA[i]
  patient_id = TRF_DF$sys.uniqueId[row_id]
  TRF_DF$smpl.dateReceived[row_id] <- unique(TRF_DF$sys.date_created[TRF_DF$sys.uniqueId == patient_id])
  
  remove(i,row_id,patient_id)
}
remove(date_NA)

# Convert to date format
TRF_DF$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                         TRF_DF$smpl.dateReceived), "%Y-%m-%d")

## Structure patient DOB and input current age
#----------------------------------------------
# Extract components of DOB
TRF_DF$year = gsub("(^[[:digit:]]{4})(.*)", "\\1", TRF_DF$base.dob)
TRF_DF$month = gsub("(^[[:digit:]]{4}[-])([[:digit:]]{2})(.*)", "\\2", TRF_DF$base.dob)
TRF_DF$day = gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-])(.*)", "\\2", TRF_DF$base.dob)

for (row_No in 1:nrow(TRF_DF)) {
  TRF_DF$patient.dob[row_No] <- paste(TRF_DF$month[row_No], TRF_DF$day[row_No], TRF_DF$year[row_No], sep="/")
}
TRF_DF$patient.dob <- as.Date( TRF_DF$patient.dob, "%m/%d/%Y")

# Age rounded down to nearest integer -- relative to smpl.dateReceived
TRF_DF$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                          TRF_DF$smpl.dateReceived), "%Y-%m-%d")

for (row_No in 1:nrow(TRF_DF)) {
  if (isTRUE(is.na(TRF_DF$patient.dob[row_No]) | is.na(TRF_DF$smpl.dateReceived[row_No]))) {
    TRF_DF$Age[row_No] <- NA
  } else {
    TRF_DF$Age[row_No] <- as.numeric(floor(age_calc(dob = TRF_DF$patient.dob[row_No], 
                                            enddate = TRF_DF$smpl.dateReceived[row_No], 
                                            units = "years")))
  }
}

# Subset relevant columns
#----------------------------------------------
col_extract  <- c("sys.uniqueId","smpl.gender","base.dob","Age","smpl.primaryTumorSite",
                  "smpl.dateReceived","smpl.specimenType","smpl.percentTumor","smpl.histologicalDiagnosis",
                  "smpl.assayName","smpl.reportDateReviewed")
TRF_DF <- unique(TRF_DF[,col_extract])

colnames(TRF_DF) <- c("PatientID","PatientGender","PatientDOB","PatientAge","PrimaryTumorSite",
                      "AssayDateReceived","smpl.specimenType","smpl.percentTumor","HistologicalDx",
                      "AssayName","AssayReportDateReviewed")
remove(col_extract,row_No)

assign("TRF_DF", TRF_DF, envir = .GlobalEnv)

write.table(TRF_DF, file = paste(tempdir, "TRF_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
