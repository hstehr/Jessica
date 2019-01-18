setwd("~/Documents/ClinicalDataScience_Fellowship/")

## FUNCTIONS
#---------------------------------------------- 
Match_initial <- function(DF, var.type) {
  
  DF$NCI.aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF$Protein)
  DF$NCI.var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,4})(.*)", "\\2", DF$Protein)
  DF$NCI.aa.end <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4}[_]*)([[:alpha:]]{3})(.*)", "\\2", DF$Protein)
  
  ################################
  # Manual exclusion
  ################################
  # sort(unique(DF$NCI.aa.end))
  coding.change <- which(DF$NCI.aa.end %in% c("c.3082+1G>T","c.3082+1G>A","c.634+1G>C","c.635-1G>T",
                                              "c.3082_3082+26del","c.2888-18_2888-2del"))
  if (length(coding.change) > 0) {
    DF$NCI.aa.start[coding.change] <- NA
    DF$NCI.aa.end[coding.change] <- NA
    DF$NCI.var.position[coding.change] <- NA
  }

  DF$NCI.var.position <- as.numeric(DF$NCI.var.position)
  
  NA.DF <- DF[is.na(DF$NCI.aa.start) | is.na(DF$NCI.aa.end) | is.na(DF$NCI.var.position) |
                is.na(DF$aa.start) | is.na(DF$aa.end) | is.na(DF$var.position),]
  
  DF <- suppressMessages(anti_join(DF, NA.DF))
  DF <- DF[DF$aa.start == DF$NCI.aa.start,]
  DF <- DF[DF$aa.end == DF$NCI.aa.end,]
  DF <- DF[DF$var.position == DF$NCI.var.position,]
  
  assign("NA.DF", NA.DF, envir = .GlobalEnv)
  assign("DF", DF, envir = .GlobalEnv)
}

Arm_Exclusion <- function(patient_DF, matched_DF, exclusion_DF) {
  patient.list <- sort(unique(matched_DF$sys.uniqueId))
  
  ## Extract all STAMP entries for each patient
  patient_DF <- patient_DF[patient_DF$sys.uniqueId %in% patient.list,]
  
  for (pt_num in 1:length(patient.list)) {
    patient_id <- patient.list[pt_num]
    
    DF_patient_Full <- patient_DF[patient_DF$sys.uniqueId == patient_id,]
    DF_patient_Full$smpl.vaf <- as.character(DF_patient_Full$smpl.vaf)
    DF_patient_Full$smpl.vaf.1 <- as.character(DF_patient_Full$smpl.vaf.1)
    DF_patient_Full$sys.name.1 <- as.logical(DF_patient_Full$sys.name.1)
    
    DF_patient <- matched_DF[matched_DF$sys.uniqueId == patient_id,]
    DF_patient$smpl.vaf <- as.character(DF_patient$smpl.vaf)
    DF_patient$smpl.vaf.1 <- as.character(DF_patient$smpl.vaf.1)
    DF_patient$sys.name.1 <- as.logical(DF_patient$sys.name.1)
    
    DF_patient <- suppressMessages(full_join(DF_patient, DF_patient_Full))
    remove(DF_patient_Full)
    
    patient.gene.list <- sort(unique(DF_patient$base.gene))
    arm.list <- sort(unique(DF_patient$Arm_Name))
    
    for (arm_num in 1:length(arm.list)) {
      arm_id <- arm.list[arm_num]
      DF_patient_Arm <- DF_patient[which(is.na(DF_patient$Arm_Name) | DF_patient$Arm_Name == arm_id),]
      
      DF_Exclude_Arm <- DF_Exclusion_Variants[DF_Exclusion_Variants$Arm_Name == arm_id,]
      DF_Exclude_Arm_combo <- sort(unique(paste(DF_Exclude_Arm$Gene_Name,DF_Exclude_Arm$Protein,sep="_")))
      patient_Inclusion <- sort(unique(paste(DF_patient_Arm$base.gene,DF_patient_Arm$smpl.hgvsProtein,sep="_")))
      
      for (bio_num in 1:length(patient_Inclusion)) {
        if (patient_Inclusion[bio_num] %in% DF_Exclude_Arm_combo) {
          print(paste(patient_id, " disqualified from NCI-MATCH ", arm_id, " due to ", patient_Inclusion[bio_num], 
                      " mutation." , sep=""))
          
          matched_DF <- matched_DF[(matched_DF$sys.uniqueId == patient_id & matched_DF$Arm_Name == arm_id) == FALSE,]
        }
      }  
    }
  }
  assign("matched_DF", matched_DF, envir = .GlobalEnv)
}

summary_check <- function(DF, var.type) {
  print(paste("Summary details for variant type == ", var.type, sep=""))
  print("base.gene (patient entries)")
  print(sort(unique(DF$base.gene)))
  print("Gene_Name (NCI-MATCH trials)")
  print(sort(unique(DF$Gene_Name)))
  print("hgvs.Protein (patient entries)")
  print(sort(unique(DF$smpl.hgvsProtein)))
  print("Protein (NCI-MATCH trials)")
  print(sort(unique(DF$Protein)))
  cat("\n")
}

## Load relevant files
#---------------------------------------------- 
STAMP_all_variants_QC <- 
  read.csv(file = paste("ClinicalTrialMatching/", Syapse_Export_timestamp, "_syapse_export_DF_STAMP_VariantAnno.csv",sep=""),
           header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")

# Import DF_Exclusion_Variants
PATIENT_VARIANT_REPORT_QC <- 
  import_list(paste("ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, 
                    "_QC_.xlsx", sep=""), setclass = "tbl")
invisible(capture.output(lapply(names(PATIENT_VARIANT_REPORT_QC),
                                function(x) assign("DF_Exclusion_Variants", 
                                                   PATIENT_VARIANT_REPORT_QC[["DF_Exclusion_Variants"]],
                                                   envir = .GlobalEnv))))
remove(PATIENT_VARIANT_REPORT_QC)
## STAMP does not identify fusion mutations 
DF_Exclusion_Variants <- DF_Exclusion_Variants[DF_Exclusion_Variants$Variant_Type != "Fusion" &
                                                 DF_Exclusion_Variants$Variant_Type != "CNV",]

## STAMP entries extracted based on matching of Gene_Name and Variant_Type
int_file_01 = paste(getwd(), "/ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, 
                    "_QC_Matched_Variants", sep="")
DF_Output_Patient_Variant <- 
  read.csv(file = paste(int_file_01, ".csv",sep=""), header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## STAMP entries extracted based on matching of Gene_Name
int_file_02 = paste(getwd(), "/ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, 
                    "_QC_Matched_NonHotspot", sep="")
DF_Output_Patient_NonHotspot <- 
  read.csv(file = paste(int_file_02,".csv",sep=""), header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## Match STAMP entries based on Protein grouped by trial Variant_Type 
#----------------------------------------------
DF_Output_SNV <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "SNV",]
DF_Output_Ins <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Insertion",]
DF_Output_Del <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Deletion",]
DF_Output_Delins <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Delins",]

print(paste("All cases have been re-distributed based on variant type:",
            nrow(DF_Output_Patient_Variant) == (nrow(DF_Output_SNV) + nrow(DF_Output_Ins) + nrow(DF_Output_Del) + 
                                                  nrow(DF_Output_Delins))),sep="")
cat("\n")
remove(DF_Output_Patient_Variant)

## Parse through SNV matches 
#----------------------------------------------
Match_initial(DF = DF_Output_SNV, var.type = "SNV")
DF_Output_SNV <- DF
summary_check(DF = DF_Output_SNV, var.type = "SNV")

################################
# Manual matching = no match
################################
unique(NA.DF[,74:80])

# A <- unique(NA.DF[,c(50,56,57,64:67)])
# # Subset by each unique gene
# B <- data.frame()
# sort(unique(NA.DF$Gene_Name))
# for (row_No in 1:nrow(A)) {
#   if (isTRUE(grepl("MET", A$sys.label[row_No]))) {
#     B <- rbind(B, A[row_No,])}}
# View(B)
# remove(A,B)

## Parse through Deletion matches
#----------------------------------------------
Match_initial(DF = DF_Output_Del, var.type = "Del")
DF_Output_Del <- DF
summary_check(DF = DF_Output_Del, var.type = "Del")

################################
# Manual matching = no match
################################
unique(NA.DF[,74:80])

# A <- unique(NA.DF[,c(50,56,57,64:68)])
# # Subset by each unique gene
# B <- data.frame()
# sort(unique(NA.DF$Gene_Name))
# for (row_No in 1:nrow(A)) {
#   if (isTRUE(grepl("KIT", A$sys.label[row_No]))) {
#     B <- rbind(B, A[row_No,])}}
# View(B)
# remove(A,B)

## Parse through Insertion matches
#----------------------------------------------
Match_initial(DF = DF_Output_Ins, var.type = "Ins")
# No matches
print("Summary details for variant type == Ins: NO MATCHES")
cat("\n")
DF_Output_Ins <- DF
remove(DF_Output_Ins)

## Parse through Indel matches 
#----------------------------------------------
Match_initial(DF = DF_Output_Delins, var.type = "Delins")
DF_Output_Delins <- DF
DF_Output_Delins <- DF_Output_Delins[DF_Output_Delins$smpl.hgvsProtein == DF_Output_Delins$Protein,]
summary_check(DF = DF_Output_Delins, var.type = "Delins")

unique(NA.DF[,74:80])

## Merge matches 
#----------------------------------------------
DF_Output_Patient_Variant_Matched <- rbind(DF_Output_SNV[,1:86], DF_Output_Del[,1:86], DF_Output_Delins[,1:86])
remove(DF_Output_SNV,DF_Output_Del,DF_Output_Delins,DF,NA.DF)

## Parse through exclusion criteria per patient
#----------------------------------------------
print("Following patients (Variant_Matched) are disqualified from indicated NCI-MATCH trial(s) due to exclusion variant criteria indicated")
Arm_Exclusion(patient_DF = STAMP_all_variants_QC, 
              matched_DF = DF_Output_Patient_Variant_Matched, 
              exclusion_DF = DF_Exclusion_Variants)
cat("\n")
DF_Output_Patient_Variant_Matched <- matched_DF

## Parse through NonHotspot matches 
#----------------------------------------------
DF_Output_Patient_NonHotspot_Matched <- data.frame()
for (row_No in 1:nrow(DF_Output_Patient_NonHotspot)) {
  if (isTRUE(DF_Output_Patient_NonHotspot$Function[row_No] == "nonframeshiftDeletion")) {
    if (isTRUE(DF_Output_Patient_NonHotspot$var.type[row_No] == "Delins" | 
               DF_Output_Patient_NonHotspot$var.type[row_No] == "Deletion")) {
      DF_Output_Patient_NonHotspot_Matched <- rbind(DF_Output_Patient_NonHotspot_Matched, DF_Output_Patient_NonHotspot[row_No,])
    }
  } else if (isTRUE(DF_Output_Patient_NonHotspot$Function[row_No] == "nonframeshiftInsertion")) {
    if (isTRUE(DF_Output_Patient_NonHotspot$var.type[row_No] == "Delins" | 
               DF_Output_Patient_NonHotspot$var.type[row_No] == "Insertion")) {
      DF_Output_Patient_NonHotspot_Matched <- rbind(DF_Output_Patient_NonHotspot_Matched, DF_Output_Patient_NonHotspot[row_No,])
    }
  } else {
    if (DF_Output_Patient_NonHotspot$var.type[row_No] != "Synonymous") {
      DF_Output_Patient_NonHotspot_Matched <- rbind(DF_Output_Patient_NonHotspot_Matched, DF_Output_Patient_NonHotspot[row_No,])
    }
  }
}

## Parse through exclusion criteria per patient
#----------------------------------------------
print("Following patients (NonHotspot_Matched) are disqualified from indicated NCI-MATCH trial(s) due to exclusion variant criteria indicated")
Arm_Exclusion(patient_DF = STAMP_all_variants_QC, 
              matched_DF = DF_Output_Patient_NonHotspot_Matched, 
              exclusion_DF = DF_Exclusion_Variants)
DF_Output_Patient_NonHotspot_Matched <- matched_DF

remove(matched_DF,DF_Output_Patient_NonHotspot,STAMP_all_variants_QC,row_No,
       DF_Exclusion_Variants,Arm_Exclusion,Match_initial,summary_check)
cat("\n")

## Write to local computer
#----------------------------------------------
write.csv(DF_Output_Patient_Variant_Matched,
          file = paste(int_file_01, "_FINAL.csv", sep=""), na = "NA", row.names = FALSE)

write.csv(DF_Output_Patient_NonHotspot_Matched,
          file = paste(int_file_02, "_FINAL.csv", sep=""), na = "NA", row.names = FALSE)

# Delete intermediate file
if (isTRUE(deleteIntermediateFile)) {
  if (file.exists(paste(int_file_01,".csv",sep=""))){file.remove(paste(int_file_01,".csv",sep=""))}
  if (file.exists(paste(int_file_02,".csv",sep=""))){file.remove(paste(int_file_02,".csv",sep=""))}
}

remove(int_file_01,int_file_02)
