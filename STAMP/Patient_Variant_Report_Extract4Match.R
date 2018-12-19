## Clinical trial INPUT: Patient_Variant_Report
## Matched patient INPUT: ClinicalTrialMatching/Patient_Variant_Matched.csv
## Extract candidate NCI clinical trials based on "hgvs.Protein" components of STAMP entries
## Match STAMP entries based on "Protein" componenets of clinical trials
## Output: ClinicalTrialMatching/DF_Output_Patient_Variant_Matched.csv

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load Library
#----------------------------------------------
library("reshape")
library("tidyr")
library("splitstackshape")

## Summarycheck
#---------------------------------------------- 
summary_check <- function(DF) {
  print("base.gene")
  print(sort(unique(DF$base.gene)))
  cat("\n")
  print("Gene_Name")
  print(sort(unique(DF$Gene_Name)))
  cat("\n")
  print("hgvs.Protein")
  print(sort(unique(DF$smpl.hgvsProtein)))
  cat("\n")
  print("Protein")
  print(sort(unique(DF$Protein)))
  cat("\n")
}

## Load relevant files
#---------------------------------------------- 
## Active clinical trial files
DF_Exclusion_Variants <-
  read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_Exclusion_Variants.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## STAMP entries extracted based on matching of Gene_Name and Variant_Type
DF_Output_Patient_Variant <- 
  read.csv(file = "ClinicalTrialMatching/Patient_Variant_Matched.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## Match STAMP entries based on Protein grouped by Variant_Type 
#----------------------------------------------
# sort(unique(DF_Output_Patient_Variant$Variant_Type))
DF_Output_SNV <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "SNV",]
DF_Output_Indel <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Indel",]
DF_Output_Ins <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Insertion",]
DF_Output_Del <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Deletion",]
DF_Output_Dup <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$Variant_Type == "Duplication",]

remove(DF_Output_Patient_Variant)

## Parse through SNV matches 
#----------------------------------------------
DF_Output_SNV$NCI.aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF_Output_SNV$Protein)
DF_Output_SNV$NCI.var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,4})(.*)", "\\2", DF_Output_SNV$Protein)
DF_Output_SNV$NCI.aa.end <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4})([[:alpha:]]{3})(.*)", "\\2", DF_Output_SNV$Protein)

coding.change <- which(DF_Output_SNV$NCI.aa.end == "c.3082+1G>T")
DF_Output_SNV$NCI.aa.start[coding.change] <- NA
DF_Output_SNV$NCI.aa.end[coding.change] <- NA
DF_Output_SNV$NCI.var.position[coding.change] <- NA
DF_Output_SNV$NCI.var.position <- as.numeric(DF_Output_SNV$NCI.var.position)

NA.list <- which(is.na(DF_Output_SNV$NCI.aa.start))
NA.DF <- DF_Output_SNV[NA.list,c(50,53,56,57,64:68,71:77)]
sort(unique(NA.DF$Protein))
sort(unique(NA.DF$sys.label))

DF_Output_SNV <- DF_Output_SNV[DF_Output_SNV$aa.start == DF_Output_SNV$NCI.aa.start,]
DF_Output_SNV <- DF_Output_SNV[DF_Output_SNV$aa.end == DF_Output_SNV$NCI.aa.end,]
DF_Output_SNV <- DF_Output_SNV[DF_Output_SNV$var.position == DF_Output_SNV$NCI.var.position,]
DF_Output_SNV <- DF_Output_SNV[rowSums(is.na(DF_Output_SNV)) != ncol(DF_Output_SNV),]

summary_check(DF = DF_Output_SNV)

## Parse through Duplication matches = manual examination
#----------------------------------------------
summary_check(DF = DF_Output_Dup)

## Parse through Deletion matches
#----------------------------------------------
DF_Output_Del$NCI.aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Del$Protein)
DF_Output_Del$NCI.var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,3})(.*)", "\\2", DF_Output_Del$Protein)
DF_Output_Del$NCI.aa.end <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Del$Protein)

coding.change <- which(DF_Output_Del$Protein == "TSC2 Loss")
DF_Output_Del$NCI.aa.start[coding.change] <- NA
DF_Output_Del$NCI.aa.end[coding.change] <- NA
DF_Output_Del$NCI.var.position[coding.change] <- NA
DF_Output_Del$NCI.var.position <- as.numeric(DF_Output_Del$NCI.var.position)

NA.list <- which(is.na(DF_Output_Del$NCI.aa.start))
# NA.DF <- DF_Output_Del[NA.list,c(50,53,56,57,64:68,71:77)]
NA.DF <- DF_Output_Del[NA.list,]
sort(unique(NA.DF$Protein))
sort(unique(NA.DF$sys.label))

DF_Output_Del <- DF_Output_Del[DF_Output_Del$aa.start == DF_Output_Del$NCI.aa.start,]
DF_Output_Del <- DF_Output_Del[DF_Output_Del$aa.end == DF_Output_Del$NCI.aa.end,]
DF_Output_Del <- DF_Output_Del[DF_Output_Del$var.position == DF_Output_Del$NCI.var.position,]
DF_Output_Del <- DF_Output_Del[rowSums(is.na(DF_Output_Del)) != ncol(DF_Output_Del),]

DF_Output_Del <- rbind(DF_Output_Del, NA.DF)

summary_check(DF = DF_Output_Del)

## Parse through Insertion matches
#----------------------------------------------
DF_Output_Ins$NCI.aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Ins$Protein)
DF_Output_Ins$NCI.var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,3})(.*)", "\\2", DF_Output_Ins$Protein)
DF_Output_Ins$NCI.aa.end <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Ins$Protein)
DF_Output_Ins$NCI.var.position <- as.numeric(DF_Output_Ins$NCI.var.position)

# No matches
DF_Output_Ins <- DF_Output_Ins[DF_Output_Ins$aa.start == DF_Output_Ins$NCI.aa.start,]

## Parse through Indel matches 
#----------------------------------------------
DF_Output_Indel$NCI.aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Indel$Protein)
DF_Output_Indel$NCI.var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,3})(.*)", "\\2", DF_Output_Indel$Protein)
DF_Output_Indel$NCI.aa.end <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Indel$Protein)

DF_Output_Indel$NCI.var.position <- as.numeric(DF_Output_Indel$NCI.var.position)

DF_Output_Indel <- DF_Output_Indel[DF_Output_Indel$aa.start == DF_Output_Indel$NCI.aa.start,]
DF_Output_Indel <- DF_Output_Indel[DF_Output_Indel$aa.end == DF_Output_Indel$NCI.aa.end,]
DF_Output_Indel <- DF_Output_Indel[DF_Output_Indel$var.position == DF_Output_Indel$NCI.var.position,]

DF_Output_Indel$NCI.var.position2 <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*[[:alpha:]]{3})([[:digit:]]{3})(.*)", "\\2", DF_Output_Indel$Protein)
DF_Output_Indel$NCI.aa.end2 <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*[[:alpha:]]{3}[[:digit:]]{3}delins)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Indel$Protein)

DF_Output_Indel$NCI.var.position2 <- as.numeric(DF_Output_Indel$NCI.var.position2)

DF_Output_Indel$var.position2 <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*[[:alpha:]]{3})([[:digit:]]{3})(.*)", "\\2", DF_Output_Indel$smpl.hgvsProtein)
DF_Output_Indel$aa.end2 <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,3}[_]*[[:alpha:]]{3}[[:digit:]]{3}delins)([[:alpha:]]{3})(.*)", "\\2", DF_Output_Indel$smpl.hgvsProtein)
# DF_Output_Indel[which(DF_Output_Indel$smpl.hgvsProtein == "p.Leu747_Pro753insSer"),]
DF_Output_Indel$aa.end2[c(78,79)] <- "Ser"
# DF_Output_Indel[which(DF_Output_Indel$smpl.hgvsProtein == "p.Trp557_Lys558del"),]
DF_Output_Indel$aa.end2[75] <- "del"

DF_Output_Indel$var.position2 <- as.numeric(DF_Output_Indel$var.position2)

DF_Output_Indel <- DF_Output_Indel[DF_Output_Indel$aa.end2 == DF_Output_Indel$NCI.aa.end2,]
DF_Output_Indel <- DF_Output_Indel[DF_Output_Indel$var.position2 == DF_Output_Indel$NCI.var.position2,]
DF_Output_Indel <- DF_Output_Indel[rowSums(is.na(DF_Output_Indel)) != ncol(DF_Output_Indel),]

summary_check(DF = DF_Output_Indel)

## Merge matches 
#----------------------------------------------
DF_Output_Patient_Variant_Matched <- rbind(DF_Output_SNV[,1:86], DF_Output_Dup[,1:86],
                                           DF_Output_Del[,1:86], DF_Output_Ins[,1:86],
                                           DF_Output_Indel[,1:86])

remove(DF_Output_SNV,DF_Output_Dup,DF_Output_Del,DF_Output_Ins,DF_Output_Indel,
       coding.change,NA.list,NA.DF)

## Parse through exclusion criteria per patient
#----------------------------------------------
## STAMP does not identify fusion mutations 
DF_Exclusion_Variants <- DF_Exclusion_Variants[DF_Exclusion_Variants$Variant_Type != "Fusion",]

patient.list <- sort(unique(DF_Output_Patient_Variant_Matched$sys.uniqueId))

# sort(table(DF_Output_Patient_Variant_Matched$sys.uniqueId))
# for (pt_num in 1:length(patient.list)) {
#   pt_id <- patient.list[pt_num]
#   DF <- DF_Output_Patient_Variant_Matched[DF_Output_Patient_Variant_Matched$sys.uniqueId == pt_id,]
#   if (nrow(DF) >= 8) {
#     print(pt_id)
#     print(sort(table(DF$Arm_Name)))
#   }
# }

for (pt_num in 1:length(patient.list)) {
  patient_id <- patient.list[pt_num]
  
  DF_patient <- DF_Output_Patient_Variant_Matched[DF_Output_Patient_Variant_Matched$sys.uniqueId == patient_id,]
  
  patient.gene.list <- sort(unique(DF_patient$base.gene))
  arm.list <- sort(unique(DF_patient$Arm_Name))
  
  for (arm_num in 1:length(arm.list)) {
    arm_id <- arm.list[arm_num]
    DF_patient_Arm <- DF_patient[DF_patient$Arm_Name == arm_id,]
    
    DF_Exclude_Arm <- DF_Exclusion_Variants[DF_Exclusion_Variants$Arm_Name == arm_id,]
    DF_Exclude_Arm_combo <- sort(unique(paste(DF_Exclude_Arm$Gene_Name,DF_Exclude_Arm$Protein,sep="_")))
    patient_Inclusion <- sort(unique(paste(DF_patient_Arm$Gene_Name,DF_patient_Arm$Protein,sep="_")))
    
    for (bio_num in 1:length(patient_Inclusion)) {
      if (patient_Inclusion[bio_num] %in% DF_Exclude_Arm_combo) {
        print(paste("Mutations identified in patient ", patient_id, 
                    " may exclude him/her from qualifying for NCI clinical trial ", arm_id, sep=""))
        print("Matched biomarker information from Inclusion criteria:")
        print(patient_Inclusion)
        print("Biomarker information from Exclusion criteria:")
        print(DF_Exclude_Arm_combo)
      }
    }  
  }
}

remove(DF_Exclude_Arm,DF_patient,DF_patient_Arm,pt_num,
       arm_id,arm_num,arm.list,bio_num,DF_Exclude_Arm_combo,
       patient_id,patient_Inclusion,patient.gene.list,
       patient.list)

## Write to local computer
#----------------------------------------------
write.csv(DF_Output_Patient_Variant_Matched,
          file = "ClinicalTrialMatching/Patient_Variant_Matched_FINAL.csv",
          na = "NA",
          row.names = FALSE)

# ## How deal with following DFs??????
#----------------------------------------------
# DF_Inclusion_NonHotspot_Rules <-
#   read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_Inclusion_NonHotspot_Rules.csv",
#            header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
#
# DF_Exclusion_NonHotspot_Rules <-
#   read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_Exclusion_NonHotspot_Rules.csv",
#            header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
# 
# DF_IHC_Results <-
#   read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_IHC_Results.csv",
#            header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
# 
# DF_Comments <-
#   read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_Comments.csv",
#            header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
