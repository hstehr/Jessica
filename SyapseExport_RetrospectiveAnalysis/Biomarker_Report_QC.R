## Convert to long format: Biomarker.Description, Disease.Sites
## Classification #1 (Disease.Site): most general from Disease.Sites
## Classification #2 (Disease.Group.category): based on Disease.Group
## Output: "Biomarker_Report_QC.tsv"

cat(paste("Timestamp of OnCore_Biomarker_Report processing START: ", Sys.time(), sep=""),"\n","\n")
OnCore_Biomarker_QC <- OnCore_Biomarker_Report

#################################
## CONFIRMATION
#################################

## CONFIRMATION: Current.Status are OPEN TO ACCRUAL
if (length(OnCore_Biomarker_QC$Current.Status[!grepl("OPEN TO ACCRUAL", OnCore_Biomarker_QC$Current.Status)]) > 0) {
  cat("Current status of clinical trials from Biomarker Report that are not OPEN TO ACCRUAL:","\n")
  cat(paste("NOTE: corresponding rows (n=", length(which(OnCore_Biomarker_QC$Current.Status != "OPEN TO ACCRUAL")),
            ") have been removed from downstream processing", sep=""),"\n")
  print(OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$Current.Status != "OPEN TO ACCRUAL"), ], 
        row.names = FALSE)
  
  ## Remove clinical trials not OPEN TO ACCRUAL from spreadsheet
  rowkeep <- which(OnCore_Biomarker_QC$Current.Status == "OPEN TO ACCRUAL")
  OnCore_Biomarker_QC <- OnCore_Biomarker_QC[rowkeep, ]
  
} else {
  cat("Current.Status of all clinical trials from Biomarker Report are OPEN TO ACCRUAL","\n")
}

## CONFIRMATION: Protocol.Type are Treatment
if (length(OnCore_Biomarker_QC$Protocol.Type[!grepl("Treatment", OnCore_Biomarker_QC$Protocol.Type)]) > 0) {
  cat("Protocol type of clinical trials from Biomarker Report that are not TREATMENT:","\n")
  cat(paste("NOTE: corresponding rows (n=", length(which(OnCore_Biomarker_QC$Protocol.Type != "Treatment")),
            ") have been removed from downstream processing", sep=""),"\n")
  print(OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$Protocol.Type != "Treatment"), ], 
        row.names = FALSE)
  
  ## Remove clinical trials not OPEN TO ACCRUAL from spreadsheet
  rowkeep <- which(OnCore_Biomarker_QC$Protocol.Type == "Treatment")
  OnCore_Biomarker_QC <- OnCore_Biomarker_QC[rowkeep, ]
  
} else {
  cat("Protocol.Type of all clinical trials from Biomarker Report are TREATMENT","\n")
}

## CONFIRMATION: Age.Group is Adult
OnCore_Biomarker_QC$Age.Group[which(OnCore_Biomarker_QC$Age.Group == "A")] <- "Adult"
cat(paste("Age.Group of clinical trials from Biomarker Report: ", unique(OnCore_Biomarker_QC$Age.Group), sep=""),"\n")

#################################
## Biomarker.Description
#################################

# Convert Biomarker.Description to long format
#----------------------------------------------
colname_keep <- colnames(OnCore_Biomarker_QC)
col_No <- as.numeric(ncol(OnCore_Biomarker_QC))

## Columns appended to end of dataframe
OnCore_Biomarker_QC$Biomarker <- OnCore_Biomarker_QC$Biomarker.Description
OnCore_Biomarker_QC <- cSplit(OnCore_Biomarker_QC, "Biomarker", ",", 
                              stripWhite = TRUE, type.convert="as.character",
                              drop = FALSE)

## Convert to long format 
OnCore_Biomarker_QC <- melt(OnCore_Biomarker_QC, id.vars=colname_keep)
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$variable != "Biomarker" &
                                                   !is.na(OnCore_Biomarker_QC$value)), ]

## Strip Biomarker.Description into components
#----------------------------------------------
OnCore_Biomarker_QC <- cSplit(OnCore_Biomarker_QC, "value", "|", 
                              stripWhite = TRUE, type.convert="as.character")
colnames(OnCore_Biomarker_QC)[c((col_No +2):(col_No +4))] <- c("Biomarker_GeneName","Biomarker_Condition","Biomarker_Detail")
OnCore_Biomarker_QC$Biomarker_Comment <- NA
colname_keep <- append(colname_keep, colnames(OnCore_Biomarker_QC)[c((col_No +2):(col_No +5))])
OnCore_Biomarker_QC <- subset(OnCore_Biomarker_QC, select = colname_keep)

## Extract Biomarker.Description comments
#----------------------------------------------
for (x in 1:nrow(OnCore_Biomarker_QC)) {
  if (isTRUE(grepl("-", OnCore_Biomarker_QC$Biomarker_Detail[x]))) {
    OnCore_Biomarker_QC$Biomarker_Comment[x] <- 
      sapply(strsplit(OnCore_Biomarker_QC$Biomarker_Detail[x], "-"), "[", 2)
    OnCore_Biomarker_QC$Biomarker_Detail[x] <- 
      sapply(strsplit(OnCore_Biomarker_QC$Biomarker_Detail[x], "-"), "[", 1)
  }
}

## Annotate Biomarker_Detail
#----------------------------------------------
OnCore_Biomarker_QC$Biomarker_Detail <- 
  gsub("^\\*[[:blank:]]*$", "All mutations accepted", OnCore_Biomarker_QC$Biomarker_Detail)

################################
## Manual correction: Amino acid change
################################
OnCore_Biomarker_QC$Biomarker_Detail[OnCore_Biomarker_QC$Biomarker_Detail == "V600*"] <- "Val600*"

## CONFIRMATION: Biomarker Gene alterations are within permitted terms
#----------------------------------------------
alterations_okay <- c("AMPLIFICATION", "DELETION","FUSION","MUTATION","All alterations")
rowkeep <- which(OnCore_Biomarker_QC$Biomarker_Condition %in% alterations_okay)

if (nrow(OnCore_Biomarker_QC) != length(rowkeep)) {
  rowNA <- which(!OnCore_Biomarker_QC$Biomarker_Condition %in% alterations_okay)
  
  cat("Biomarker Gene alterations from Biomarker Report that are not defined using the following terms:","\n")
  print(alterations_okay)
  cat(paste("NOTE: corresponding rows (n=", 
            length(rowNA), ") have been removed from downstream processing", sep=""),"\n")
  print(OnCore_Biomarker_QC[rowNA,c(1:15)], row.names = FALSE)
  
  ## Remove respective rows from spreadsheet
  OnCore_Biomarker_QC <- OnCore_Biomarker_QC[rowkeep, ]
  
} else {
  cat("Biomarker Gene alterations from Biomarker Report are defined using the following term(s):","\n")
  print(alterations_okay)
}

#################################
## Disease.Group
#################################
# Convert to lowercase 
OnCore_Biomarker_QC$Disease.Group <- tolower(OnCore_Biomarker_QC$Disease.Group)

# Confirm Disease.Group is classified 
Disease.Group.Report <- sort(unique(OnCore_Biomarker_QC$Disease.Group))
Disease.Group.key <- sort(unique(DiseaseGroupCategory_LongFormat$Disease.Group))

for (elem_No in 1:length(Disease.Group.Report)) {
  if (isTRUE(is.element(Disease.Group.Report[elem_No], Disease.Group.key) == FALSE)) {
    cat(paste(Disease.Group.Report[elem_No], 
              ": Disease.Group has not been classified into Disease.Group.category", sep=""),"\n")
    
  }
}

# Assign Disease.Group.category based on Disease.Group
OnCore_Biomarker_QC$Disease.Group.category <- NA
for (row_No in 1:nrow(OnCore_Biomarker_QC)) {
  DiseaseGroupCategory.name <- 
    unique(DiseaseGroupCategory_LongFormat$Disease.Group.category[which(DiseaseGroupCategory_LongFormat$Disease.Group == 
                                                                          OnCore_Biomarker_QC$Disease.Group[row_No])]) 
  
  OnCore_Biomarker_QC$Disease.Group.category[row_No] <- paste(as.character(DiseaseGroupCategory.name), collapse=", ")
}

# Convert Disease.Group.category to long format
#----------------------------------------------
## Strip Disease.Group.category into components
OnCore_Biomarker_QC <- cSplit(OnCore_Biomarker_QC, "Disease.Group.category", ",", 
                              stripWhite = TRUE, type.convert="as.character",
                              drop = FALSE)

## Convert to long format 
OnCore_Biomarker_QC <- melt(OnCore_Biomarker_QC, id.vars=colname_keep)
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$variable != "Disease.Group.category" &
                                                   !is.na(OnCore_Biomarker_QC$value)), ]
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[,c(colname_keep, "value")]
colnames(OnCore_Biomarker_QC)[as.numeric(ncol(OnCore_Biomarker_QC))] <- c("Disease.Group.category")
colname_keep <- append(colname_keep, "Disease.Group.category")

#################################
## Disease.Sites
#################################
# Convert to lowercase 
OnCore_Biomarker_QC$Disease.Sites <- tolower(OnCore_Biomarker_QC$Disease.Sites)

# Create new column
OnCore_Biomarker_QC$Disease.Site <- OnCore_Biomarker_QC$Disease.Sites

################################
# Manual correction for consistency & ontology
################################
OnCore_Biomarker_QC$Disease.Site[OnCore_Biomarker_QC$Disease.Site == 
                                   "leukemia, other; other hematopoietic"] <- "other leukemia; other hematopoietic"

OnCore_Biomarker_QC$Disease.Site <- 
  gsub("other leukemia; other hematopoietic","blood",OnCore_Biomarker_QC$Disease.Site)
OnCore_Biomarker_QC$Disease.Site <- 
  gsub("other digestive organ",
       DiseaseGroupCategory$primaryTumorSite[DiseaseGroupCategory$Disease.Group.category == "gastroenterology"],
       OnCore_Biomarker_QC$Disease.Site)

# Most general Disease.Site is assigned priority 
OnCore_Biomarker_QC$Disease.Site[grep("any site", OnCore_Biomarker_QC$Disease.Site)] <- "any site"

# Standarize notation for new element
OnCore_Biomarker_QC$Disease.Site <- gsub(",|and",";", OnCore_Biomarker_QC$Disease.Site)

# Convert Disease.Site.category to long format
#----------------------------------------------
## Strip Disease.Site.category into components
OnCore_Biomarker_QC <- cSplit(OnCore_Biomarker_QC, "Disease.Site", ";", 
                              stripWhite = TRUE, type.convert="as.character",
                              drop = FALSE)

## Convert to long format 
OnCore_Biomarker_QC <- melt(OnCore_Biomarker_QC, id.vars=colname_keep)
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$variable != "Disease.Site" &
                                                   !is.na(OnCore_Biomarker_QC$value)), ]
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[,c(colname_keep, "value")]
colnames(OnCore_Biomarker_QC)[as.numeric(ncol(OnCore_Biomarker_QC))] <- c("Disease.Site")
colname_keep <- append(colname_keep, "Disease.Site")

################################
## Format dataframe 
################################
# Select columns of interest
colname_keep <- c("OnCore.No","Age.Group","Disease.Group.category","Disease.Site",
                  "Biomarker_GeneName","Biomarker_Condition","Biomarker_Detail")
OnCore_Biomarker_QC <- OnCore_Biomarker_QC[order(OnCore_Biomarker_QC$OnCore.No),colname_keep]

## Write variable in global environment
#----------------------------------------------
assign("OnCore_Biomarker_QC", OnCore_Biomarker_QC, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(OnCore_Biomarker_QC, paste(tempdir,OnCore_Biomarker_Report_timestamp, "_Biomarker_Report_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(alterations_okay,col_No,colname_keep,Disease.Group.key,Disease.Group.Report,
       DiseaseGroupCategory.name,elem_No,row_No,rowkeep,x)
cat("\n")
cat(paste("Timestamp of OnCore_Biomarker_Report processing FINISH: ", Sys.time(), sep=""),"\n","\n")
