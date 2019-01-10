setwd("~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/")

## Load relevant file
#----------------------------------------------
OnCore_Biomarker_Report <- read.csv(file = paste("Biomarker Report_", OnCore_Biomarker_Report_timestamp, ".csv", sep=""),
                                    header = TRUE, na.strings = c("NA", ""),
                                    stringsAsFactors = FALSE, sep = ",")

print(paste("Timestamp of OnCore_Biomarker_Report: ", OnCore_Biomarker_Report_timestamp, sep=""))

## Confirm Current.Status are OPEN TO ACCRUAL
#----------------------------------------------
if (length(OnCore_Biomarker_Report$Current.Status[!grepl("OPEN TO ACCRUAL", OnCore_Biomarker_Report$Current.Status)]) > 0) {
  print("Current status of clinical trials from Biomarker Report that are not OPEN TO ACCRUAL:")
  print(paste("NOTE: corresponding rows (n=", length(which(OnCore_Biomarker_Report$Current.Status != "OPEN TO ACCRUAL")),
  ") have been removed from downstream processing", sep=""))
  cat("\n")
  print(OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$Current.Status != "OPEN TO ACCRUAL"), ], 
        row.names = FALSE)
  
  ## Remove clinical trials not OPEN TO ACCRUAL from spreadsheet
  rowkeep <- which(OnCore_Biomarker_Report$Current.Status == "OPEN TO ACCRUAL")
  OnCore_Biomarker_Report <- OnCore_Biomarker_Report[rowkeep, ]
} else {
  print("Current.Status of all clinical trials from Biomarker Report are OPEN TO ACCRUAL")
}

## Confirm Protocol.Type are Treatment
#----------------------------------------------
if (length(OnCore_Biomarker_Report$Protocol.Type[!grepl("Treatment", OnCore_Biomarker_Report$Protocol.Type)]) > 0) {
  print("Protocol type of clinical trials from Biomarker Report that are not TREATMENT:")
  print(paste("NOTE: corresponding rows (n=", length(which(OnCore_Biomarker_Report$Protocol.Type != "Treatment")),
              ") have been removed from downstream processing", sep=""))
  cat("\n")
  print(OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$Protocol.Type != "Treatment"), ], 
        row.names = FALSE)

  ## Remove clinical trials not OPEN TO ACCRUAL from spreadsheet
  rowkeep <- which(OnCore_Biomarker_Report$Protocol.Type == "Treatment")
  OnCore_Biomarker_Report <- OnCore_Biomarker_Report[rowkeep, ]
} else {
  print("Protocol.Type of all clinical trials from Biomarker Report are TREATMENT")
}


## Match Age.Group column with patient DOB
#----------------------------------------------
for (x in 1:nrow(OnCore_Biomarker_Report)) {
  if (OnCore_Biomarker_Report$Age.Group[x] == "A") {
    OnCore_Biomarker_Report$Age.Group[x] <- "Adult"
  }
}
print(paste("Age.Group of clinical trials from Biomarker Report: ", 
            unique(OnCore_Biomarker_Report$Age.Group), sep=""))

# Convert Biomarker.Description to long format
#----------------------------------------------
colname_keep <- colnames(OnCore_Biomarker_Report)

## Columns appended to end of dataframe
OnCore_Biomarker_Report$Biomarker <- OnCore_Biomarker_Report$Biomarker.Description
OnCore_Biomarker_Report <- cSplit(OnCore_Biomarker_Report, "Biomarker", ",", 
                                  stripWhite = TRUE, type.convert="as.character",
                                  drop = FALSE)

## Convert to long format 
OnCore_Biomarker_Report <- melt(OnCore_Biomarker_Report, 
                                id.vars=colname_keep)
OnCore_Biomarker_Report <- 
  OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$variable != "Biomarker" &
                                  !is.na(OnCore_Biomarker_Report$value)), ]

## Strip Biomarker.Description into components
#----------------------------------------------
OnCore_Biomarker_Report <- cSplit(OnCore_Biomarker_Report, "value", "|", 
                                  stripWhite = TRUE, type.convert="as.character")
colnames(OnCore_Biomarker_Report)[17:19] <- 
  c("Biomarker_GeneName","Biomarker_Condition","Biomarker_Detail")
OnCore_Biomarker_Report$Biomarker_Comment <- NA
colname_keep <- append(colname_keep, colnames(OnCore_Biomarker_Report)[17:20])
OnCore_Biomarker_Report <- subset(OnCore_Biomarker_Report, select = colname_keep)

## Confirm Biomarker Gene alterations are within permitted terms
#----------------------------------------------
alterations_okay <- c("AMPLIFICATION", "DELETION","FUSION","MUTATION","All alterations")
rowkeep <- which(OnCore_Biomarker_Report$Biomarker_Condition %in% alterations_okay)
                   
if (nrow(OnCore_Biomarker_Report) != length(rowkeep)) {
  rowNA <- which(!OnCore_Biomarker_Report$Biomarker_Condition %in% alterations_okay)

  print("Biomarker Gene alterations from Biomarker Report that are not defined using the following terms:")
  print(alterations_okay)
  print(paste("NOTE: corresponding rows (n=", 
              length(rowNA), ") have been removed from downstream processing", sep=""))
  cat("\n")
  print(OnCore_Biomarker_Report[rowNA,c(1:15)], row.names = FALSE)
  
  ## Remove respective rows from spreadsheet
  OnCore_Biomarker_Report <- OnCore_Biomarker_Report[rowkeep, ]
} else {
  print("Biomarker Gene alterations from Biomarker Report are defined using the following term(s):")
  print(alterations_okay)
}
remove(rowkeep)

## Extract Biomarker.Description comments
#----------------------------------------------
for (x in 1:nrow(OnCore_Biomarker_Report)) {
  if (isTRUE(grepl("-", OnCore_Biomarker_Report$Biomarker_Detail[x]))) {
    OnCore_Biomarker_Report$Biomarker_Comment[x] <- 
      sapply(strsplit(OnCore_Biomarker_Report$Biomarker_Detail[x], "-"), "[", 2)
    OnCore_Biomarker_Report$Biomarker_Detail[x] <- 
      sapply(strsplit(OnCore_Biomarker_Report$Biomarker_Detail[x], "-"), "[", 1)
  }
}

OnCore_Biomarker_Report$Biomarker_Detail <- 
  gsub("^\\*[[:blank:]]*$", "All mutations accepted", OnCore_Biomarker_Report$Biomarker_Detail)

################################
## Manual correction: Amino acid change
################################
OnCore_Biomarker_Report$Biomarker_Detail[OnCore_Biomarker_Report$Biomarker_Detail == "V600*"] <- "Val600*"

# Disease.Group grouped according to medical specialties  
# Most general Disease.Site is assigned priority 
################################
# Manual assignment
################################
OnCore_Biomarker_Report$Disease.Group.category <- NA
for (row_No in 1:nrow(OnCore_Biomarker_Report)) {
  if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Cutaneous Oncology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Dermatology"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Developmental Therapeutics" |
             OnCore_Biomarker_Report$Disease.Group[row_No] == "Non-CRG Specific") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Any Site"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Gastrointestinal Oncology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Gastroenterology"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Genitourinary Oncology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Genitourinary"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Head & Neck Oncology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Otolaryngology"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Thoracic Oncology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Pulmonology"
  } else if (OnCore_Biomarker_Report$Disease.Group[row_No] == "Hematology") {
    OnCore_Biomarker_Report$Disease.Group.category[row_No] <- "Hematology"
  }
}

OnCore_Biomarker_Report$Disease.Site.category <- OnCore_Biomarker_Report$Disease.Sites
OnCore_Biomarker_Report$Disease.Site.category[grep("Any Site", OnCore_Biomarker_Report$Disease.Sites)] <- "Any Site"

################################
# Manual correction for consistency
################################
OnCore_Biomarker_Report$Disease.Site.category [OnCore_Biomarker_Report$Disease.Site.category == 
                                                 "Leukemia, other; Other Hematopoietic"] <- "Other Leukemia; Other Hematopoietic"
  
# Convert Disease.Site.category to long format
#----------------------------------------------
colname_keep <- append(colname_keep, "Disease.Group.category")
## Strip Disease.Site.category into components
OnCore_Biomarker_Report <- cSplit(OnCore_Biomarker_Report, "Disease.Site.category", ";", 
                                  stripWhite = TRUE, type.convert="as.character",
                                  drop = FALSE)

## Convert to long format 
OnCore_Biomarker_Report <- melt(OnCore_Biomarker_Report, 
                                id.vars=append(colname_keep, "Disease.Site.category"))
OnCore_Biomarker_Report <- OnCore_Biomarker_Report[which(!is.na(OnCore_Biomarker_Report$value)), ]

colnames(OnCore_Biomarker_Report)[as.numeric(ncol(OnCore_Biomarker_Report))] <- c("Disease.Site")
colname_keep <- append(colname_keep, "Disease.Site")
OnCore_Biomarker_Report <- subset(OnCore_Biomarker_Report, select = colname_keep)

# ## Data Review
# #----------------------------------------------
# print(paste("Total number of unique OnCore.No ", length(unique(OnCore_Biomarker_Report$OnCore.No)), sep=""))
# print(paste("Total number of unique disease groups: ", length(sort(unique(OnCore_Biomarker_Report$Disease.Group))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Disease.Group))
# print(paste("Total number of unique disease group categories: ", length(sort(unique(OnCore_Biomarker_Report$Disease.Group.category))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Disease.Group.category))
# print(paste("Total number of unique disease sites: ", length(sort(unique(OnCore_Biomarker_Report$Disease.Site))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Disease.Site))
# print(paste("Total number of unique biomarker genes: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName))
# print(paste("Total number of unique biomarker conditions: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_Condition))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_Condition))
# print(paste("Total number of unique Biomarker details: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_Detail))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_Detail))

# Confirm all entries have a Disease.Group.category and Disease.Site
# which(is.na(OnCore_Biomarker_Report$Disease.Group.category))
# sort(table(OnCore_Biomarker_Report$Disease.Group.category))
# which(is.na(OnCore_Biomarker_Report$Disease.Site))
# sort(table(OnCore_Biomarker_Report$Disease.Site))

remove(alterations_okay,colname_keep,x,row_No)
cat("\n")

## Write to local computer
#----------------------------------------------
write.csv(OnCore_Biomarker_Report,
          file = paste("Biomarker_Report_LongFormat_", OnCore_Biomarker_Report_timestamp, ".csv", sep=""),
          na = "NA", row.names = FALSE)
