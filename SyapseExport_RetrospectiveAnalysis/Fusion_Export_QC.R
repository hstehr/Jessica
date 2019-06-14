# Fusion export for STAMP v2 (130 genes)
#----------------------------------------------
Fusion.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-05-03_Syapse_export_STAMP2_Fusions.csv"
STAMP_Fusion <- read.csv(file = Fusion.file, header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")
Fusion_Export_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*/", "", Fusion.file))), format= "%Y-%m-%d")

# Shorten column names of DF
colnames(STAMP_Fusion) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(STAMP_Fusion))

## Remove duplicates of entries based on multiple unique sys.date_created
## Include cases of smpl.amendedString == "AMENDED"
#----------------------------------------------
STAMP_Fusion_AMENDED <- data.frame()
patient_list <- unique(STAMP_Fusion$sys.uniqueId)

for (id_num in 1:length(patient_list)) {
  patient_id <- patient_list[id_num]
  
  assay_list <- unique(STAMP_Fusion$smpl.assayName[which(STAMP_Fusion$sys.uniqueId == patient_id)])
  
  for (assay_num in 1:length(assay_list)) {
    assay_id <- assay_list[assay_num]
    
    if (length(unique(STAMP_Fusion$sys.date_created[which(STAMP_Fusion$sys.uniqueId == patient_id & STAMP_Fusion$smpl.assayName == assay_id)])) > 1) {
      created_list <- unique(STAMP_Fusion$sys.date_created[which(STAMP_Fusion$sys.uniqueId == patient_id & STAMP_Fusion$smpl.assayName == assay_id)])
      
      ## Retain most recent version if multiple sys.date_created.1 exist in report 
      STAMP_Fusion_AMENDED_pre <- STAMP_Fusion[which(STAMP_Fusion$sys.uniqueId == patient_id & STAMP_Fusion$smpl.assayName == assay_id),]
      STAMP_Fusion_AMENDED_pre <- STAMP_Fusion_AMENDED_pre[STAMP_Fusion_AMENDED_pre$sys.date_created == max(created_list),]
      
      remove(created_list)
      
    } else {
      STAMP_Fusion_AMENDED_pre <- STAMP_Fusion[which(STAMP_Fusion$sys.uniqueId == patient_id & STAMP_Fusion$smpl.assayName == assay_id),]
    }
    
    STAMP_Fusion_AMENDED <- rbind(STAMP_Fusion_AMENDED, STAMP_Fusion_AMENDED_pre)
    remove(assay_id,assay_num)
  }
  remove(patient_id,assay_list,id_num)
}

# Overwrite existing dataframe
STAMP_Fusion <- STAMP_Fusion_AMENDED
remove(STAMP_Fusion_AMENDED,STAMP_Fusion_AMENDED_pre,patient_list)

# Use sys.date_created as alternative for missing smpl.dateReceived
date_NA <- which(is.na(STAMP_Fusion$smpl.dateReceived))
for (i in 1:length(date_NA)) {
  row_id = date_NA[i]
  patient_id = STAMP_Fusion$sys.uniqueId[row_id]
  STAMP_Fusion$smpl.dateReceived[row_id] <- unique(STAMP_Fusion$sys.date_created[STAMP_Fusion$sys.uniqueId == patient_id])
  
  remove(i,row_id,patient_id)
}
remove(date_NA)

## Structure patient DOB and input current age
#----------------------------------------------
STAMP_Fusion$base.dob <- as.Date(STAMP_Fusion$base.dob, "%Y-%m-%d")

# Age rounded down to nearest integer -- relative to smpl.dateReceived
STAMP_Fusion$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                               STAMP_Fusion$smpl.dateReceived), "%Y-%m-%d")

STAMP_Fusion$Age <- as.numeric(floor(age_calc(dob = STAMP_Fusion$base.dob, 
                                              enddate = STAMP_Fusion$smpl.dateReceived, units = "years")))

# MANUAL EDITS
#----------------------------------------------
STAMP_Fusion$sys.label[which(STAMP_Fusion$sys.label == "TRK-FUSED GENE; TFG-ALK")] <- "TFG-ALK"
STAMP_Fusion$sys.label[which(STAMP_Fusion$sys.label == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6-RET")] <- "CCDC6-RET"
# sort(unique(STAMP_Fusion$sys.label))

STAMP_Fusion$smpl.fusionGene1[which(STAMP_Fusion$smpl.fusionGene1 == "TRK-FUSED GENE; TFG")] <- "TFG"
STAMP_Fusion$smpl.fusionGene1[which(STAMP_Fusion$smpl.fusionGene1 == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6")] <- "CCDC6"
# sort(unique(append(STAMP_Fusion$smpl.fusionGene1, STAMP_Fusion$smpl.fusionGene2)))

# Rename columns
#----------------------------------------------
colnames_extract <- c("sys.uniqueId","smpl.gender","base.dob","Age",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","smpl.fusionGene1","smpl.fusionGene2")
STAMP_Fusion <- STAMP_Fusion[,colnames_extract]

colnames_generic <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                      "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "AssayName","sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "Fusion_Detail","Gene1","Gene2")
colnames(STAMP_Fusion) <- colnames_generic
remove(colnames_extract,colnames_generic)

# Reformat to long format i.e. single gene per row 
#----------------------------------------------
colname_keep <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                  "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
                  "AssayName","sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
                  "smpl.amendedString","smpl.amendmentReason","Fusion_Detail")

## Convert to long format 
STAMP_Fusion <- melt(STAMP_Fusion, id.vars=colname_keep)
STAMP_Fusion <- STAMP_Fusion[,c(colname_keep,"value")]
colnames(STAMP_Fusion) <- c(colname_keep,"Gene")

# General structurization: convert to lowercase 
#----------------------------------------------
STAMP_Fusion$HistologicalDx <- tolower(STAMP_Fusion$HistologicalDx)
STAMP_Fusion$PrimaryTumorSite <- tolower(STAMP_Fusion$PrimaryTumorSite)

# PrimaryTumorSite grouped according to medical specializations (PrimaryTumorSite.Category)
#----------------------------------------------
# Confirm PrimaryTumorSite is classified 
primaryTumorSite.STAMP <- sort(unique(STAMP_Fusion$PrimaryTumorSite))
primaryTumorSite.key <- sort(unique(DiseaseGroupCategory_LongFormat$primaryTumorSite))

# for (elem_No in 1:length(primaryTumorSite.STAMP)) {
#   if (isTRUE(is.element(primaryTumorSite.STAMP[elem_No], primaryTumorSite.key) == FALSE)) {
#     cat("\n")
#     cat(paste(primaryTumorSite.STAMP[elem_No],
#               ": primaryTumorSite has not been classified into DiseaseGroupCategory", sep=""),"\n","\n")
#   }
# }

STAMP_Fusion$PrimaryTumorSite.Category <- NA
for (row_No in 1:nrow(STAMP_Fusion)) {
  DiseaseGroupCategory.name <- 
    DiseaseGroupCategory_LongFormat$Disease.Group.category[which(DiseaseGroupCategory_LongFormat$primaryTumorSite == 
                                                                   STAMP_Fusion$PrimaryTumorSite[row_No])] 
  
  STAMP_Fusion$PrimaryTumorSite.Category[row_No] <- paste(as.character(DiseaseGroupCategory.name), collapse=", ")
}

# If missing smpl.primaryTumorSite, primaryTumorSite.category = "unknown"
STAMP_Fusion$PrimaryTumorSite.Category[which(STAMP_Fusion$PrimaryTumorSite.Category == "")] <- "unknown"

## Remove entries with missing information = gender & DOB
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(!is.na(STAMP_Fusion$PatientGender)),]
STAMP_Fusion <- STAMP_Fusion[!is.na(STAMP_Fusion$PatientDOB),]

# Missing information
#----------------------------------------------
STAMP_Fusion$var.type <- "Fusion"
STAMP_Fusion$var.anno <- "Fusion"
STAMP_Fusion$VariantPathogenicityStatus <- "NULL"
STAMP_Fusion$Break <- "NULL"

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_Fusion", STAMP_Fusion, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(STAMP_Fusion, file = paste(tempdir, Fusion_Export_timestamp, "_Fusion_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Export entries per patient into "_SNVIndel.tsv" file
#----------------------------------------------
patient.list <- sort(unique(STAMP_Fusion$PatientID))

for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Extract STAMP entries for individual patient
  DF_patient <- STAMP_Fusion[which(STAMP_Fusion$PatientID == patient_id),]
  
  # Order alphabetically by gene name
  DF_patient <- DF_patient[order(DF_patient$Fusion_Detail, decreasing = FALSE),]
  
  ## Write to local computer
  write.table(DF_patient, file = paste(tempdir, patient_id, "_Fusion.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
