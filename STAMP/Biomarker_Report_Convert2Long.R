## Clean up Biomarker Report_2018-10_OnCore_Biomarker_Report.csv
## Specific parameters converted to long format 
## Biomarker.Description, Disease.Sites
## Output: "Biomarker_Report_LongFormat.csv"

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/")

## Load Library
#----------------------------------------------
library("splitstackshape")
library("reshape")

## Load relevant file
#----------------------------------------------
OnCore_Biomarker_Report <- read.csv(file = "Biomarker Report_2018-10.csv",
                                    header = TRUE,
                                    na.strings = c("NA", ""),
                                    stringsAsFactors = FALSE,
                                    sep = ",")

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


## Should match Age.Group column with patient DOB
#----------------------------------------------
for (x in 1:nrow(OnCore_Biomarker_Report)) {
  if (OnCore_Biomarker_Report$Age.Group[x] == "A") {
    OnCore_Biomarker_Report$Age.Group[x] <- "Adult"
  }
}

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
  print("Biomarker Gene alterations from Biomarker Report are all defined using the following terms:")
  print(alterations_okay)
}
remove(rowkeep)

## Extract Biomarker.Description comments
#######################
for (x in 1:nrow(OnCore_Biomarker_Report)) {
  if (isTRUE(grepl("-", OnCore_Biomarker_Report$Biomarker_Detail[x]))) {
    OnCore_Biomarker_Report$Biomarker_Comment[x] <- 
      sapply(strsplit(OnCore_Biomarker_Report$Biomarker_Detail[x], "-"), "[", 2)
    OnCore_Biomarker_Report$Biomarker_Detail[x] <- 
      sapply(strsplit(OnCore_Biomarker_Report$Biomarker_Detail[x], "-"), "[", 1)
  }
}

###############################################################
## What does "V600*" mean? Input as separate field?????????
###############################################################
OnCore_Biomarker_Report$Biomarker_Detail <- 
  gsub("^\\*[[:blank:]]*$", "All mutations accepted", OnCore_Biomarker_Report$Biomarker_Detail)
OnCore_Biomarker_Report$Biomarker_Detail <- 
  gsub("\\*$", " (All mutations accepted)", OnCore_Biomarker_Report$Biomarker_Detail)


# Convert Disease.Sites to long format
#----------------------------------------------
## Strip Disease.Sites into components
OnCore_Biomarker_Report <- cSplit(OnCore_Biomarker_Report, "Disease.Sites", ";", 
                                  stripWhite = TRUE, type.convert="as.character",
                                  drop = FALSE)

## Convert to long format 
OnCore_Biomarker_Report <- melt(OnCore_Biomarker_Report, 
                                id.vars=colname_keep)
OnCore_Biomarker_Report <- OnCore_Biomarker_Report[which(!is.na(OnCore_Biomarker_Report$value)), ]

colnames(OnCore_Biomarker_Report)[21] <- c("Disease.Site")
colname_keep <- append(colname_keep, "Disease.Site")
OnCore_Biomarker_Report <- subset(OnCore_Biomarker_Report, select = colname_keep)

# ## Simple review of data
# #----------------------------------------------
# print(paste("Total number of unique OnCore.No ", length(unique(OnCore_Biomarker_Report$OnCore.No)), sep=""))
# print(paste("Total number of unique disease groups: ", length(sort(unique(OnCore_Biomarker_Report$Disease.Group))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Disease.Group))
# print(paste("Total number of unique disease sites: ", length(sort(unique(OnCore_Biomarker_Report$Disease.Site))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Disease.Site))
# print(paste("Total number of unique biomarker genes: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName))
# print(paste("Total number of unique biomarker conditions: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_Condition))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_Condition))
# print(paste("Total number of unique Biomarker details: ", length(sort(unique(OnCore_Biomarker_Report$Biomarker_Detail))), sep=""))
# sort(unique(OnCore_Biomarker_Report$Biomarker_Detail))

## Write to local computer
#----------------------------------------------
write.csv(OnCore_Biomarker_Report,
          file = "Biomarker_Report_LongFormat.csv",
          na = "NA",
          row.names = FALSE)
