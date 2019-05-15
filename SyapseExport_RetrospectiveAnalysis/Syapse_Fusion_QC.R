Fusion.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-05-03_Syapse_export_STAMP2_Fusions.csv"
STAMP_Fusion <- read.csv(file = Fusion.file, header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")
Fusion_Export_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*/", "", Fusion.file))), format= "%Y-%m-%d")

## Output: "syapse_export_fusion_QC.csv"

cat(paste("Timestamp of Syapse_Fusions: ", Fusion_Export_timestamp, sep=""),"\n","\n")

# Load relevant file from global environment 
#----------------------------------------------
DF_Full <- STAMP_Fusion

# Shorten column names of DF
colnames(DF_Full) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_Full))

# Subset columns of interest
colnames_extract <- c("sys.uniqueId","base.dob","smpl.gender",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","smpl.hasOrderingPhysician",
                      "sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","smpl.fusionGene1","smpl.fusionGene2")

DF_Full <- DF_Full[,colnames_extract]

## Remove duplicates of entries based on multiple unique sys.date_created
## Include cases of smpl.amendedString == "AMENDED"
#----------------------------------------------
DF_Full_AMENDED <- data.frame()
patient_list <- unique(DF_Full$sys.uniqueId)

for (id_num in 1:length(patient_list)) {
  patient_id <- patient_list[id_num]
  
  assay_list <- unique(DF_Full$smpl.assayName[which(DF_Full$sys.uniqueId == patient_id)])
  
  for (assay_num in 1:length(assay_list)) {
    assay_id <- assay_list[assay_num]
    
    if (length(unique(DF_Full$sys.date_created[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id)])) > 1) {
      created_list <- unique(DF_Full$sys.date_created[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id)])
      
      ## Retain most recent version if multiple sys.date_created.1 exist in report 
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id),]
      DF_Full_AMENDED_pre <- DF_Full_AMENDED_pre[DF_Full_AMENDED_pre$sys.date_created == max(created_list),]
      
      remove(created_list)
      
    } else {
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id),]
    }
    
    DF_Full_AMENDED <- rbind(DF_Full_AMENDED, DF_Full_AMENDED_pre)
    remove(assay_id,assay_num)
  }
  remove(patient_id,assay_list,id_num)
}

# Overwrite existing dataframe
DF_Full <- DF_Full_AMENDED
remove(DF_Full_AMENDED,DF_Full_AMENDED_pre)

# Use sys.date_created as alternative for missing smpl.dateReceived
date_NA <- which(is.na(DF_Full$smpl.dateReceived))
for (i in 1:length(date_NA)) {
  row_id = date_NA[i]
  patient_id = DF_Full$sys.uniqueId[row_id]
  DF_Full$smpl.dateReceived[row_id] <- unique(DF_Full$sys.date_created[DF_Full$sys.uniqueId == patient_id])
  
  remove(i,row_id,patient_id)
}
remove(date_NA)

## Remove entries with missing information
#----------------------------------------------
## Remove entries without gender = 0 entries
DF_Full <- DF_Full[which(!is.na(DF_Full$smpl.gender)),]

# Remove patients without DOB = 0 entries
DF_Full <- DF_Full[!is.na(DF_Full$base.dob),]

## Structure patient DOB and input current age
#----------------------------------------------
DF_Full$base.dob <- as.Date(DF_Full$base.dob, "%Y-%m-%d")

# Age rounded down to nearest integer -- relative to smpl.dateReceived
DF_Full$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                          DF_Full$smpl.dateReceived), "%Y-%m-%d")

DF_Full$Age <- as.numeric(floor(age_calc(dob = DF_Full$base.dob, 
                                         enddate = DF_Full$smpl.dateReceived, units = "years")))

# Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
DF_Full <- DF_Full[which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]

## Format dataframe 
#----------------------------------------------
# Rename columns
colnames_extract <- c("sys.uniqueId","smpl.gender","base.dob","Age",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","smpl.hasOrderingPhysician",
                      "sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","smpl.fusionGene1","smpl.fusionGene2")
DF_Full <- DF_Full[,colnames_extract]

colnames_generic <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                      "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "AssayName","AssayOrderingPhysician",
                      "sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "Fusion_Detail","Gene1","Gene2")
colnames(DF_Full) <- colnames_generic

# General structurization
#----------------------------------------------
# Convert to lowercase 
DF_Full$HistologicalDx <- tolower(DF_Full$HistologicalDx)
DF_Full$PrimaryTumorSite <- tolower(DF_Full$PrimaryTumorSite)

# MANUAL EDITS
#----------------------------------------------
DF_Full$Fusion_Detail[which(DF_Full$Fusion_Detail == "TRK-FUSED GENE; TFG-ALK")] <- "TFG-ALK"
DF_Full$Fusion_Detail[which(DF_Full$Fusion_Detail == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6-RET")] <- "CCDC6-RET"
# sort(unique(DF_Full$Fusion_Detail))

DF_Full$Gene1[which(DF_Full$Gene1 == "TRK-FUSED GENE; TFG")] <- "TFG"
DF_Full$Gene1[which(DF_Full$Gene1 == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6")] <- "CCDC6"
# sort(unique(DF_Full$Gene1))

# sort(unique(DF_Full$Gene2))

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_Fusion", DF_Full, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(DF_Full, file = paste(tempdir, Fusion_Export_timestamp, "_Syapse_Fusion_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(patient_list,colnames_generic,colnames_extract,DF_Full)
