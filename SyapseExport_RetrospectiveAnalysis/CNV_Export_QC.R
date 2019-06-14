## STAMP CNV database QC PIPELINE
#----------------------------------------------
## All CNVs marked as "AMP" (i.e. only high-confidence amplifications) for each TRF id from STAMP2. 
## If the TRF is not listed in the file, 
# then either there were no amplifications, or the sample could not be matched for some reason.  
## Existence of sites.addition.CNV not addressed 

# Fusion export for STAMP v2 (130 genes)
#----------------------------------------------
CNV.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-05-15_stamp_cnvs.txt"
STAMP_CNV <- read.csv(file = CNV.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
CNV_Export_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*/", "", CNV.file))), format= "%Y-%m-%d")

# Rename column names of DF
colnames(STAMP_CNV) <- c("PatientID","report_version","CNV_Gene","var.anno")

# Integrate with primary tumor site field from TRF_File
#----------------------------------------------
STAMP_CNV <- left_join(STAMP_CNV,
                       unique(TRF_DF[,c("PatientID","PatientGender","PatientDOB","PatientAge","PrimaryTumorSite",
                                        "AssayDateReceived","AssayReportDateReviewed",
                                        "smpl.specimenType","smpl.percentTumor","HistologicalDx","AssayName")]),
                       by = "PatientID")

# Collapse similar gene names
#----------------------------------------------
STAMP_CNV$CNV_Gene[which(tolower(STAMP_CNV$CNV_Gene) %in% c("nkx2", "nkx2-1"))] <- "NKX2-1"
# sort(unique(STAMP_CNV$CNV_Gene))

# General structurization: convert to lowercase 
#----------------------------------------------
STAMP_CNV$PrimaryTumorSite <- tolower(STAMP_CNV$PrimaryTumorSite)
STAMP_CNV$HistologicalDx <- tolower(STAMP_CNV$HistologicalDx)

# PrimaryTumorSite grouped according to medical specializations (PrimaryTumorSite.Category)
#----------------------------------------------
# Confirm PrimaryTumorSite is classified 
primaryTumorSite.STAMP <- sort(unique(STAMP_CNV$PrimaryTumorSite))
primaryTumorSite.key <- sort(unique(DiseaseGroupCategory_LongFormat$primaryTumorSite))

# for (elem_No in 1:length(primaryTumorSite.STAMP)) {
#   if (isTRUE(is.element(primaryTumorSite.STAMP[elem_No], primaryTumorSite.key) == FALSE)) {
#     cat("\n")
#     cat(paste(primaryTumorSite.STAMP[elem_No],
#               ": primaryTumorSite has not been classified into DiseaseGroupCategory", sep=""),"\n","\n")
#   }
# }

STAMP_CNV$PrimaryTumorSite.Category <- NA
for (row_No in 1:nrow(STAMP_CNV)) {
  DiseaseGroupCategory.name <- 
    DiseaseGroupCategory_LongFormat$Disease.Group.category[which(DiseaseGroupCategory_LongFormat$primaryTumorSite == 
                                                                   STAMP_CNV$PrimaryTumorSite[row_No])] 
  
  STAMP_CNV$PrimaryTumorSite.Category[row_No] <- paste(as.character(DiseaseGroupCategory.name), collapse=", ")
}

# If missing smpl.primaryTumorSite, primaryTumorSite.category = "unknown"
STAMP_CNV$PrimaryTumorSite.Category[which(STAMP_CNV$PrimaryTumorSite.Category == "")] <- "unknown"

## Remove entries with missing information = gender & DOB
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[which(!is.na(STAMP_CNV$PatientGender)),]
STAMP_CNV <- STAMP_CNV[!is.na(STAMP_CNV$PatientDOB),]

# Missing information
#----------------------------------------------
STAMP_CNV$VariantPathogenicityStatus <- "NULL"
STAMP_CNV$var.type <- "CNV"
STAMP_CNV$Locus <- "NULL"
STAMP_CNV$Tiles <- "NULL"
STAMP_CNV$mean.z <- "NULL"
STAMP_CNV$mcopies <- "NULL"

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_CNV", STAMP_CNV, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(STAMP_CNV, file = paste(tempdir, CNV_Export_timestamp, "_CNV_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Export entries per patient into "_SNVIndel.tsv" file
#----------------------------------------------
patient.list <- sort(unique(STAMP_CNV$PatientID))

for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Extract STAMP entries for individual patient
  DF_patient <- STAMP_CNV[which(STAMP_CNV$PatientID == patient_id),]
  
  # Order alphabetically by gene name
  DF_patient <- DF_patient[order(DF_patient$CNV_Gene, decreasing = FALSE),]
  
  ## Write to local computer
  write.table(DF_patient, file = paste(tempdir, patient_id, "_CNV.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
