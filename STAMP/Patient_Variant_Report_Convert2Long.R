## Clean up PATIENT_VARIANT_REPORT_TEMPLATE.xlsx
####################################################################################
## Biomarker.Description, Disease.Sites << XXX????
####################################################################################
## Output: "Patient_Variant_Report_Inclusion_NonHotspot_Rules.csv"
## Output: "Patient_Variant_Report_Exclusion_NonHotspot_Rules.csv"
## Output: "Patient_Variant_Report_Exclusion_Variants.csv"
## Output: "Patient_Variant_Report_Inclusion_Variants.csv"
## Output: "Patient_Variant_Report_IHC_Results.csv"
## Output: "Patient_Variant_Report_Comments.csv"

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/")

## Load Library
#----------------------------------------------
library("dplyr")
library("tidyr")
library("stringr")
library("rio")
library("tidyverse")

## Manual correction of fields
#----------------------------------------------
Field_Correction <- function(DF) {
  DF_add <- DF[DF$Variant_Type == "CNV" & 
                 grepl("Amplification", DF$Protein, ignore.case = TRUE) == TRUE,]
  if (nrow(DF_add) > 0) {
    DF_add$Variant_Type <- "Duplication"
  }
  
  DF_add_02 <- data.frame()
  DF_add_03 <- data.frame()
  if (length(which(DF$Variant_Type == "CNV" & 
                   grepl("CNV", DF$Protein, ignore.case = TRUE) == TRUE)) > 0) {
    row_No_CNV <- as.numeric(which(DF$Variant_Type == "CNV" & 
                                     grepl("CNV", DF$Protein, ignore.case = TRUE) == TRUE))
    
    DF_add_02 <- DF[row_No_CNV,]
    DF_add_03 <- DF_add_02
    
    if (nrow(DF_add_02) > 0) {
      DF_add_02$Variant_Type <- "Duplication"
      DF_add_03$Variant_Type <- "Insertion"
    }
    
    DF$Variant_Type[row_No_CNV] <- "Deletion"
    remove(row_No_CNV)
  }
  
  for (row_No in 1:nrow(DF)) {
    if (DF$Variant_Type[row_No] == "CNV" & 
        grepl("Loss", DF$Protein[row_No], ignore.case = TRUE) == TRUE) {
      DF$Variant_Type[row_No] <- "Deletion"
    } else if (DF$Variant_Type[row_No] == "CNV" & 
               grepl("Amplification", DF$Protein[row_No], ignore.case = TRUE) == TRUE) {
      DF$Variant_Type[row_No] <- "Insertion"
    } else if (DF$Variant_Type[row_No] == "SNV") {
      DF$Protein[row_No] <- gsub("(^p.[[:blank:]]*)(.*)", "\\2", DF$Protein[row_No])
      DF$Protein[row_No] <- paste("p.", DF$Protein[row_No], sep="")
    }
  }
  DF <- rbind(DF,DF_add,DF_add_02,DF_add_03)
  remove(DF_add,DF_add_02,DF_add_03,row_No)
  
  # Remove empty rows 
  DF <- DF[rowSums(is.na(DF)) != ncol(DF),]
  
  assign("DF", DF, envir = .GlobalEnv)
}

## Load relevant file
#----------------------------------------------
## Import excel file as single list
PATIENT_VARIANT_REPORT <- import_list("PATIENT_VARIANT_REPORT_TEMPLATE.xlsx", setclass = "tbl")
## Generate individual DF per worksheet in excel file
# invisible(capture.output(lapply(names(PATIENT_VARIANT_REPORT), 
#                                 function(x) assign(x, PATIENT_VARIANT_REPORT[[x]],
#                                                    envir = .GlobalEnv))))
# remove(`ARM-GENE LOOK-UP TABLE`,Version,`Disease Exclusion LOOK-UP Table`,
#        `GENE-STRAND Look-up table`,`Gene Ref Seq Look-up table`,`EXAMPLE DATA for Standard aMOIs`,
#        `Patient ID information`,MOIs,`lab list`,`Variant Type`)
# 
# ## Extract names of dataframes in global environment
# names(PATIENT_VARIANT_REPORT)
# # Total number of Arms
# length(unique(PATIENT_VARIANT_REPORT$`ARM-GENE LOOK-UP TABLE`$ARM))

## Specify parameters of output files 
#----------------------------------------------
DF_Inclusion_NonHotspot_Rules <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Inclusion_NonHotspot_Rules) <- c("Arm_Name", "Description", "oncominevariantclass",
                                             "Gene_Name","Exon","Function",	"Level_of_Evidence",
                                             "PRESENT", "Protein","VARIANT_TYPE","Chromosome",
                                             "position","Ref","Alt","ALLELE_FREQ","Variant_ID")

DF_Exclusion_NonHotspot_Rules <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Exclusion_NonHotspot_Rules) <- colnames(DF_Inclusion_NonHotspot_Rules)

DF_Exclusion_Variants <- data.frame(matrix(NA, ncol = 15))
colnames(DF_Exclusion_Variants) <- c("Arm_Name", "Gene_Name","Variant_ID","Variant_Type",
                                     "Protein","Level_of_Evidence","HGVS","Chromosome",
                                     "Position","Ref","Alt","PRESENT","ALLELE_FREQ",
                                     "CN_value","FUSION_READ_DEPTH")

DF_Inclusion_Variants <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Inclusion_Variants) <- append(colnames(DF_Exclusion_Variants), "Variant Source")

DF_IHC_Results <- data.frame(matrix(NA, ncol = 7))
colnames(DF_IHC_Results) <- c("Arm_Name","Gene","Status_POSITIVE_NEGATIVE_INDETERMINATE",
                           "Variant_PRESENT_NEGATIVE_EMPTY","Description","LOE","IHC_RESULT")

DF_Comments <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Comments) <- c("Arm_Name","Gene","X__1","X__2","X__3",
                           "X__4","X__5","X__6","X__7","X__8","X__9",
                           "X__10","X__11","X__12","X__13","X__14")


## Parse info across ARMs from original file into corresponding output files
#----------------------------------------------
#################################################################################
## What is the last worksheet preceding & following Arm info????
#################################################################################
arm_start = which(names(PATIENT_VARIANT_REPORT) == "Patient ID information") +1
arm_end = which(names(PATIENT_VARIANT_REPORT) == "MOIs") -1

for (Arm_No in arm_start:arm_end) {
  DF <- PATIENT_VARIANT_REPORT[[Arm_No]]
  Arm_Name <- colnames(DF)[1]
  ## Remove rows that are all empty
  DF <- DF[rowSums(is.na(DF)) != ncol(DF),]  

  ##################################################################
  ## is relevant tab always column 2??? How consistent is structure???
  ##################################################################
  row_end = (which(DF[[2]] == "Exclusion Variants")) -1
  
  if (isTRUE(row_end > 2)) {
    if (isTRUE(which(DF[[2]] == "Exclusion Non-Hotspot Rules") > 0)) {
      row_int_end = (which(DF[[2]] == "Exclusion Non-Hotspot Rules")) -1
      row_int_start = row_int_end +3
      
      DF_Inclusion_NonHotspot_Rules_pre <- DF[c(2:row_int_end),1:16]
      colnames(DF_Inclusion_NonHotspot_Rules_pre) <- colnames(DF_Inclusion_NonHotspot_Rules) 
      DF_Inclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Inclusion_NonHotspot_Rules <- rbind(DF_Inclusion_NonHotspot_Rules, DF_Inclusion_NonHotspot_Rules_pre)
      
      DF_Exclusion_NonHotspot_Rules_pre <- DF[c(row_int_start:row_end),1:16]
      colnames(DF_Exclusion_NonHotspot_Rules_pre) <- colnames(DF_Exclusion_NonHotspot_Rules) 
      DF_Exclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Exclusion_NonHotspot_Rules <- rbind(DF_Exclusion_NonHotspot_Rules, DF_Exclusion_NonHotspot_Rules_pre)
      
    } else {
      DF_Inclusion_NonHotspot_Rules_pre <- DF[c(2:row_end),1:16]
      colnames(DF_Inclusion_NonHotspot_Rules_pre) <- colnames(DF_Inclusion_NonHotspot_Rules) 
      DF_Inclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Inclusion_NonHotspot_Rules <- rbind(DF_Inclusion_NonHotspot_Rules, DF_Inclusion_NonHotspot_Rules_pre)
    }
  }
   
  row_start = (which(DF[[2]] == "Exclusion Variants")) +2
  row_end = (which(DF[[2]] == "Inclusion Variants")) -1
  DF_Exclusion_Variants_pre <- DF[c(row_start:row_end),1:15]
  colnames(DF_Exclusion_Variants_pre) <- colnames(DF_Exclusion_Variants) 
  DF_Exclusion_Variants_pre$Arm_Name <- Arm_Name
  DF_Exclusion_Variants <- rbind(DF_Exclusion_Variants, DF_Exclusion_Variants_pre)
  
  row_start = (which(DF[[2]] == "Inclusion Variants")) +2
  row_end = (which(DF[[2]] == "COMMENTS:")) -1
  
  if (isTRUE(which(DF[[2]] == "IHC Results") > 0)) {
    row_int_end = (which(DF[[2]] == "IHC Results")) -1
    row_int_start = row_int_end +3
    
    DF_Inclusion_Variants_pre <- DF[c(row_start:row_int_end),1:16]
    colnames(DF_Inclusion_Variants_pre) <- colnames(DF_Inclusion_Variants) 
    DF_Inclusion_Variants_pre$Arm_Name <- Arm_Name
    DF_Inclusion_Variants <- rbind(DF_Inclusion_Variants, DF_Inclusion_Variants_pre)

    DF_IHC_Results_pre <- DF[c(row_int_start:row_end),1:7]
    colnames(DF_IHC_Results_pre) <- colnames(DF_IHC_Results) 
    DF_IHC_Results_pre$Arm_Name <- Arm_Name
    DF_IHC_Results <- rbind(DF_IHC_Results, DF_IHC_Results_pre)
    
  } else {
    DF_Inclusion_Variants_pre <- DF[c(row_start:row_end),1:16]
    colnames(DF_Inclusion_Variants_pre) <- colnames(DF_Inclusion_Variants) 
    DF_Inclusion_Variants_pre$Arm_Name <- Arm_Name
    DF_Inclusion_Variants <- rbind(DF_Inclusion_Variants, DF_Inclusion_Variants_pre)
  }

  row_start = (which(DF[[2]] == "COMMENTS:")) +1
  row_end = nrow(DF)
  if (row_start < row_end) {
    DF_Comments_pre <- DF[c(row_start:row_end),1:16]
    colnames(DF_Comments_pre) <- colnames(DF_Comments) 
    DF_Comments_pre$Arm_Name <- Arm_Name
    DF_Comments<- rbind(DF_Comments, DF_Comments_pre)
  }

  if (exists("DF_Inclusion_NonHotspot_Rules_pre")) {remove(DF_Inclusion_NonHotspot_Rules_pre)}
  if (exists("DF_Exclusion_NonHotspot_Rules_pre")) {remove(DF_Exclusion_NonHotspot_Rules_pre)}
  if (exists("DF_Exclusion_Variants_pre")) {remove(DF_Exclusion_Variants_pre)}
  if (exists("DF_Inclusion_Variants_pre")) {remove(DF_Inclusion_Variants_pre)}
  if (exists("DF_IHC_Results_pre")) {remove(DF_IHC_Results_pre)}
  if (exists("DF_Comments_pre")) {remove(DF_Comments_pre)}
  remove(Arm_Name,row_start,row_end)
  }

## Remove rows that are all empty or where column 2 == c("none", "None")
DF_Inclusion_NonHotspot_Rules <- 
  DF_Inclusion_NonHotspot_Rules[rowSums(is.na(DF_Inclusion_NonHotspot_Rules)) != ncol(DF_Inclusion_NonHotspot_Rules),]

DF_Exclusion_NonHotspot_Rules <- 
  DF_Exclusion_NonHotspot_Rules[rowSums(is.na(DF_Exclusion_NonHotspot_Rules)) != 
                                  ncol(DF_Exclusion_NonHotspot_Rules),]

DF_Exclusion_Variants <- 
  DF_Exclusion_Variants[rowSums(is.na(DF_Exclusion_Variants)) != ncol(DF_Exclusion_Variants) &
                          DF_Exclusion_Variants[[2]] != "None", ]  

DF_Inclusion_Variants <- 
  DF_Inclusion_Variants[rowSums(is.na(DF_Inclusion_Variants)) != 
                                  ncol(DF_Inclusion_Variants),]  

DF_IHC_Results <- DF_IHC_Results[rowSums(is.na(DF_IHC_Results)) != ncol(DF_IHC_Results) &
                                   DF_IHC_Results[[2]] != "none",]  

DF_Comments <- DF_Comments[rowSums(is.na(DF_Comments)) != ncol(DF_Comments),]

remove(PATIENT_VARIANT_REPORT, DF)

#################################################################
## How deal with notes under COMMENT section? What are the headers??? 
#################################################################
## Manual edit of DF_Comments
DF_Comments$Note[max(which(DF_Comments$Arm_Name == "ARM-T"))] <- 
  paste(DF_Comments$Gene[1], DF_Comments$Gene[2], DF_Comments$Gene[3], sep=" ")
DF_Comments <- DF_Comments[DF_Comments$Arm_Name == "ARM-T" & 
                      !is.na(DF_Comments$Note),]
DF_Comments <- DF_Comments[,c(1:11,ncol(DF_Comments))]
colnames(DF_Comments)[1:11] <- colnames(DF_Exclusion_Variants)[1:11]

# ## Simple review of data
# #----------------------------------------------
# table(sort(DF_Inclusion_Variants$Gene_Name))
# table(sort(DF_Inclusion_Variants$Variant_Type))
# table(sort(DF_Exclusion_Variants$Gene_Name))
# table(sort(DF_Exclusion_Variants$Variant_Type))
# table(sort(DF_Inclusion_NonHotspot_Rules$Gene_Name))
# table(sort(DF_Exclusion_NonHotspot_Rules$Gene_Name))
# table(sort(DF_IHC_Results$Gene))
# table(sort(DF_Comments$Gene_Name))

remove(arm_end,arm_start,Arm_No,row_int_end,row_int_start)

## Manual correction of fields = DF_Inclusion_Variants
#----------------------------------------------
Field_Correction(DF = DF_Inclusion_Variants)
DF_Inclusion_Variants <- DF

# DF_Inclusion_Variants[which(DF_Inclusion_Variants$Variant_Type == "MNV"), c(1:11)]
DF_Inclusion_Variants$Variant_Type[c(7,72)] <- "SNV"
DF_Inclusion_Variants$Protein[c(7,72)] <- "p.Gly719Cys"

# DF_Inclusion_Variants[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082G>C"), c(1:11)]
DF_Inclusion_Variants$Protein[65] <- "p.Asp1028His"

# Input coding info for intronic mutations
# DF_Inclusion_Variants[which(DF_Inclusion_Variants$HGVS == "NM_001127500.2:c.3082+1G>T"), c(1:11)]
DF_Inclusion_Variants$Protein[63] <- "c.3082+1G>T"

# DF_Inclusion_Variants$HGVS == "NM_001127500.2:c.3082_3082+26del"
DF_Inclusion_Variants$Protein[64] <- "c.3082_3082+26del"

# DF_Inclusion_Variants[which(DF_Inclusion_Variants$Protein == "p.G983_G1054del"),]
DF_Inclusion_Variants$Protein[c(60,80)] <- "p.Gly983_Gly1054del"

# Annotate for consistency  
for (row_No in 1:nrow(DF_Inclusion_Variants)) {
  if (DF_Inclusion_Variants$Variant_Type[row_No] == "Indel") {
    if (grepl("del.*ins", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Indel"
    } else if (grepl("ins", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Insertion"
    } else if (grepl("del", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Deletion"
    }
  }
}

## Manual correction of fields = DF_Exclusion_Variants
#----------------------------------------------
Field_Correction(DF = DF_Exclusion_Variants)
DF_Exclusion_Variants <- DF

remove(DF,row_No)

## Write to local computer
#----------------------------------------------
write.csv(DF_Inclusion_Variants,
          file = "Patient_Variant_Report_Inclusion_Variants.csv",
          na = "NA",
          row.names = FALSE)

write.csv(DF_Exclusion_Variants,
          file = "Patient_Variant_Report_Exclusion_Variants.csv",
          na = "NA",
          row.names = FALSE)

# write.csv(DF_Inclusion_NonHotspot_Rules,
#           file = "Patient_Variant_Report_Inclusion_NonHotspot_Rules.csv",
#           na = "NA",
#           row.names = FALSE)
# 
# write.csv(DF_Exclusion_NonHotspot_Rules,
#           file = "Patient_Variant_Report_Exclusion_NonHotspot_Rules.csv",
#           na = "NA",
#           row.names = FALSE)
# 
# write.csv(DF_IHC_Results,
#           file = "Patient_Variant_Report_IHC_Results.csv",
#           na = "NA",
#           row.names = FALSE)
# 
# write.csv(DF_Comments,
#           file = "Patient_Variant_Report_Comments.csv",
#           na = "NA",
#           row.names = FALSE)
