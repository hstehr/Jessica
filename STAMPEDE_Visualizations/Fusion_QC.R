sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP Fusion database QC PIPELINE
#################################

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

## Remove entries with missing information = gender & DOB
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(!is.na(STAMP_Fusion$smpl.gender)),]
STAMP_Fusion <- STAMP_Fusion[!is.na(STAMP_Fusion$base.dob),]

# General structurization: convert to lowercase 
#----------------------------------------------
STAMP_Fusion$smpl.histologicalDiagnosis <- tolower(STAMP_Fusion$smpl.histologicalDiagnosis)
STAMP_Fusion$smpl.primaryTumorSite <- tolower(STAMP_Fusion$smpl.primaryTumorSite)

# MANUAL EDITS
#----------------------------------------------
STAMP_Fusion$sys.label[which(STAMP_Fusion$sys.label == "TRK-FUSED GENE; TFG-ALK")] <- "TFG-ALK"
STAMP_Fusion$sys.label[which(STAMP_Fusion$sys.label == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6-RET")] <- "CCDC6-RET"
# sort(unique(STAMP_Fusion$sys.label))

STAMP_Fusion$smpl.fusionGene1[which(STAMP_Fusion$smpl.fusionGene1 == "TRK-FUSED GENE; TFG")] <- "TFG"
STAMP_Fusion$smpl.fusionGene1[which(STAMP_Fusion$smpl.fusionGene1 == "COILED-COIL DOMAIN-CONTAINING PROTEIN 6; CCDC6")] <- "CCDC6"
# sort(unique(append(STAMP_Fusion$smpl.fusionGene1, STAMP_Fusion$smpl.fusionGene2)))

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# sort(unique(STAMP_Fusion$smpl.histologicalDiagnosis))

# Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_Fusion$smpl.assayName))

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.primaryTumorSite),]
STAMP_Fusion <- STAMP_Fusion[!(STAMP_Fusion$smpl.primaryTumorSite %in% c("unknown","none","other primary site")), ]

# Collapse similar primary tumor site
STAMP_Fusion$smpl.primaryTumorSite[which(STAMP_Fusion$smpl.primaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_Fusion$smpl.primaryTumorSite[which(STAMP_Fusion$smpl.primaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"

# Abbreviate for STAMPEDE display
STAMP_Fusion$smpl.primaryTumorSite[which(STAMP_Fusion$smpl.primaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
STAMP_Fusion$smpl.primaryTumorSite[which(STAMP_Fusion$smpl.primaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"

STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.primaryTumorSite), ]
# sort(unique(STAMP_Fusion$smpl.primaryTumorSite))

# Filter entries for complete cases = percent tumor
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(!is.na(STAMP_Fusion$smpl.percentTumor)),]
STAMP_Fusion <- STAMP_Fusion[which(grepl("comment", STAMP_Fusion$smpl.percentTumor, ignore.case = TRUE) == FALSE),]
STAMP_Fusion$smpl.percentTumor <- gsub("%$", "", STAMP_Fusion$smpl.percentTumor)
STAMP_Fusion$smpl.percentTumor <- gsub("^>", "", STAMP_Fusion$smpl.percentTumor)
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$smpl.percentTumor != "None"), ]

STAMP_Fusion$smpl.percentTumor <- ceiling(as.numeric(STAMP_Fusion$smpl.percentTumor))
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.percentTumor),]
# sort(unique(STAMP_Fusion$smpl.percentTumor))

# Filter entries for complete cases = specimen type
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$smpl.specimenType != "other"),]
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.specimenType),]
# sort(unique(STAMP_Fusion$smpl.specimenType))

# # Filter for adults
# #----------------------------------------------
# STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$Age >= 18),]
# sort(as.numeric(unique(STAMP_Fusion$Age)))

# Rename columns
#----------------------------------------------
colnames_extract <- c("sys.uniqueId","smpl.gender","base.dob","Age",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","smpl.hasOrderingPhysician",
                      "sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","smpl.fusionGene1","smpl.fusionGene2")
STAMP_Fusion <- STAMP_Fusion[,colnames_extract]

colnames_generic <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                      "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "AssayName","AssayOrderingPhysician",
                      "sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "Fusion_Detail","Gene1","Gene2")
colnames(STAMP_Fusion) <- colnames_generic
remove(colnames_extract,colnames_generic)

# Examine number of missing fields
#----------------------------------------------
Fusion.list <- sort(unique(STAMP_Fusion$PatientID))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PatientID)))])
cat(paste("Number of Fusion entries without corresponding PatientID in POST-Filter STAMP_DF: ",
          Fusion.diff_No,sep=""),"\n","\n")

Fusion.diff_No <- length(which(is.na(STAMP_Fusion$PrimaryTumorSite)))
cat(paste("Number of Fusion entries missing PrimaryTumorSite: ",Fusion.diff_No,sep=""),"\n","\n")

Fusion.list <- sort(unique(STAMP_Fusion$PrimaryTumorSite))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))])
cat(paste("Number of Fusion entries without corresponding PrimaryTumorSite in POST-Filter STAMP_DF: ",
          Fusion.diff_No,sep=""),"\n",
    paste(unlist(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]),collapse=", "),"\n")

# Append elements not found in SNV/Indel DF from Fusion DF
Fusion.list <- Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]
if (length(Fusion.list) > 0) {
  sites.addition.Fusion = Fusion.list
} else {
  sites.addition.Fusion = NULL
}

Fusion.list <- sort(unique(append(STAMP_Fusion$Gene1,STAMP_Fusion$Gene2)))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))])
Fusion.list <- Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))]
cat(paste("Number of Fusion entries without corresponding Variant Gene in POST-Filter STAMP_DF: ",
          Fusion.diff_No, sep=""),"\n",
    paste(unlist(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))]),collapse=", "),"\n","\n")

# Append elements not found in SNV/Indel DF from Fusion DF
if (length(Fusion.list) > 0) {
  genes.addition.Fusion = Fusion.list
} else {
  genes.addition.Fusion = NULL
}

remove(Fusion.list,Fusion.diff_No,Fusion.file)

cat(paste("POST-Filter counts: ",nrow(STAMP_Fusion), " total entries and ", 
          length(unique(STAMP_Fusion[[1]])), " total test orders", sep=""),"\n","\n")

assign("STAMP_Fusion", STAMP_Fusion, envir = .GlobalEnv)
write.table(STAMP_Fusion, file = paste(tempdir, Fusion_Export_timestamp, "_Fusion_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
