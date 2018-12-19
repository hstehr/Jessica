## Clinical trial INPUT: ClinicalTrialMatching/Biomarker_Report_LongFormat.csv
## Matched patient INPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched.csv
## Extract candidate clinical trials based on gene_id and disease_site of STAMP entries
## Match STAMP entries based on Biomarker.Description of clinical trials
## Output: ClinicalTrialMatching/OnCore_Biomarker_Matched_merged.csv

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load Library
#----------------------------------------------
library("reshape")
library("tidyr")
library("splitstackshape")
library("plyr")

## Load relevant active clinical trial files
#----------------------------------------------
OnCore_Biomarker_Report <- 
  read.csv(file = "ClinicalTrialMatching/Biomarker_Report_LongFormat.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
for (row_No in 1:nrow(OnCore_Biomarker_Report)) {
  if (is.na(OnCore_Biomarker_Report$Biomarker_Comment[row_No])) {
    OnCore_Biomarker_Report$Biomarker_Comment[row_No] <- "None"
  }
}

## STAMP entries extracted based on matching of gene name and disease site
DF_Output_OnCore_Biomarker <- 
  read.csv(file = "ClinicalTrialMatching/OnCore_Biomarker_Matched.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
ncol_No2Keep = as.numeric(which(colnames(DF_Output_OnCore_Biomarker) == "current.age") +1)
colnames(DF_Output_OnCore_Biomarker)[[ncol_No2Keep]] <- "disease.site_matched"

## Remove duplicated clinical trials (extended melt procedure)
row_total <- as.numeric(ncol(DF_Output_OnCore_Biomarker))
for (row_No in 1:nrow(DF_Output_OnCore_Biomarker)) {
  DF <- DF_Output_OnCore_Biomarker[row_No,(ncol_No2Keep +1):row_total]
  DF <- DF[!duplicated(as.list(DF))]
  nexist = as.numeric(ncol(DF))
  ntotal = row_total - ncol_No2Keep
  DF[,nexist+1:(ntotal - nexist)] <- NA
  
  DF_Output_OnCore_Biomarker[row_No,(ncol_No2Keep +1):row_total] <- DF
}
remove(row_total,ntotal,nexist,DF)

# Convert matched clinical trial to long format
#----------------------------------------------
colname_keep <- colnames(DF_Output_OnCore_Biomarker)[1:ncol_No2Keep]
DF_Output_OnCore_Biomarker <- melt(DF_Output_OnCore_Biomarker, 
                                   id.vars=colname_keep)
DF_Output_OnCore_Biomarker <- 
  DF_Output_OnCore_Biomarker[which(!is.na(DF_Output_OnCore_Biomarker$value)), 
                             c(colname_keep,"value")]
DF_Output_OnCore_Biomarker$value <- as.character(DF_Output_OnCore_Biomarker$value)

# Input associated clinical trial info from OnCore_Biomarker_Report
#----------------------------------------------
# Format structure of DF_Output_OnCore_Biomarker
ncol_original <- as.numeric(ncol((DF_Output_OnCore_Biomarker)))
DF_Output_OnCore_Biomarker <- cbind(DF_Output_OnCore_Biomarker,
                                    data.frame(matrix(NA, ncol = ncol(OnCore_Biomarker_Report), 
                                                      nrow = nrow(DF_Output_OnCore_Biomarker))))
colnames(DF_Output_OnCore_Biomarker) <- c(colnames(DF_Output_OnCore_Biomarker)[1:ncol_original],
                                          colnames(OnCore_Biomarker_Report))

## Specify parameters of output files 
#----------------------------------------------
DF_Output_Biomarker_Matched <- data.frame(matrix(NA, ncol = ncol(DF_Output_OnCore_Biomarker)))
colnames(DF_Output_Biomarker_Matched) <- colnames(DF_Output_OnCore_Biomarker) 

# Extract corresponding row based on OnCore.No
for (match_num in 1:nrow(DF_Output_OnCore_Biomarker)) {
  match_id <- DF_Output_OnCore_Biomarker$value[match_num]
  gene_id <- DF_Output_OnCore_Biomarker$base.gene[match_num]
  site_id <- DF_Output_OnCore_Biomarker$disease.site_matched[match_num]
  
  match_info <- OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$OnCore.No == match_id &
                        OnCore_Biomarker_Report$Biomarker_GeneName == gene_id &
                          OnCore_Biomarker_Report$Disease.Site == site_id),]
  
  if (nrow(match_info) > 1) {

    ## reformat structure of match_info
    col_extract <- c("Biomarker_GeneName", "Biomarker_Condition","Biomarker_Detail",
                     "Biomarker_Comment","Disease.Site")
    
    ncol_original_match = as.numeric(ncol(match_info))
    col_new = (nrow(match_info) -1)*length(col_extract) +ncol_original_match
    match_info[,c((ncol_original_match+1):col_new)] <- NA

    row_end = ncol(DF_Output_OnCore_Biomarker) + (col_new -ncol_original_match)

    for (Criteria_No in 2:nrow(match_info)) {
      row_end_int = (Criteria_No -1)*length(col_extract) +ncol_original_match
      row_start_int = row_end_int -(length(col_extract) -1)

      match_info[1,c(row_start_int:row_end_int)] <- match_info[Criteria_No,col_extract]
    }
    DF_Output_Pre <- cbind(DF_Output_OnCore_Biomarker[match_num,c(1:ncol_original)],match_info[1,])

  } else {
    DF_Output_Pre <- cbind(DF_Output_OnCore_Biomarker[match_num,c(1:ncol_original)],match_info)
  }
  
  DF_Output_Biomarker_Matched <- rbind.fill(DF_Output_Biomarker_Matched, DF_Output_Pre)
}

DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[rowSums(is.na(DF_Output_Biomarker_Matched)) != ncol(DF_Output_Biomarker_Matched),]

# which(DF_Output_Biomarker_Matched$value != DF_Output_Biomarker_Matched$OnCore.No)
DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[,!colnames(DF_Output_Biomarker_Matched) == "value"]

remove(DF_Output_OnCore_Biomarker,DF_Output_Pre,
       match_info,OnCore_Biomarker_Report,col_new,
       Criteria_No,gene_id,match_id,match_num,ncol_original,
       ncol_original_match,row_end,row_end_int,row_No,
       row_start_int,site_id)

# Convert Biomarker.Description to long format
#----------------------------------------------
ncol_No2Keep = as.numeric(which(colnames(DF_Output_Biomarker_Matched) == "Title"))
colname_keep <- colnames(DF_Output_Biomarker_Matched)[1:ncol_No2Keep]

# Merge individual Biomarker.Description = MANUAL INPUT
col_start <- as.numeric(which(colnames(DF_Output_Biomarker_Matched) == "Title"))
biomarker_iter = (ncol(DF_Output_Biomarker_Matched) - col_start) %/% (length(col_extract))

for (n_iter in biomarker_iter:1) {
  col_int_start = col_start +1 +(length(col_extract)*(n_iter -1))
  col_int_end = col_int_start +length(col_extract) -1
  
  print(paste("For n_iter=", n_iter, " col_start=", col_int_start, 
              " and col_end=", col_int_end, sep=""))
}

# For n_iter=3 col_start=97 and col_end=101
DF_Output_Biomarker_Matched <- unite(DF_Output_Biomarker_Matched, 
                                     newcol.3, c(97:101), sep = "_", remove = TRUE)
# For n_iter=2 col_start=92 and col_end=96
DF_Output_Biomarker_Matched <- unite(DF_Output_Biomarker_Matched, 
                                     newcol.2, c(92:96), sep = "_", remove = TRUE)
# For n_iter=1 col_start=87 and col_end=91
DF_Output_Biomarker_Matched <- unite(DF_Output_Biomarker_Matched, 
                                     newcol.1, c(87:91), sep = "_", remove = TRUE)

# Convert Biomarker.Description to long format
DF_Output_Biomarker_Matched <- melt(DF_Output_Biomarker_Matched, 
                                    id.vars=colname_keep)

DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[which(DF_Output_Biomarker_Matched$value != "NA_NA_NA_NA_NA"), 
                             c(colname_keep,"value")]
DF_Output_Biomarker_Matched$value <- as.character(DF_Output_Biomarker_Matched$value)

# Strip Biomarker.Description into components
DF_Output_Biomarker_Matched <- cSplit(DF_Output_Biomarker_Matched, "value", "_", 
                                  stripWhite = TRUE, type.convert="as.character")
colnames(DF_Output_Biomarker_Matched)[(col_start +1):ncol(DF_Output_Biomarker_Matched)] <- col_extract 

# Match var.anno from patient with Biomarker_Condition from clinical trial reports
DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[which(DF_Output_Biomarker_Matched$var.anno == 
                                      DF_Output_Biomarker_Matched$Biomarker_Condition),]

###############################################################
# Match with Biomarker_Detail from clinical trial reports??????
# "All mutations accepted", "Amino acid change", "Translocation partner gene"
###############################################################
remove(biomarker_iter,col_extract,col_int_start,col_int_end,
       col_start,colname_keep,n_iter,ncol_No2Keep)

## Write file on local computer
#----------------------------------------------
write.csv(DF_Output_Biomarker_Matched,
          file = "ClinicalTrialMatching/OnCore_Biomarker_Matched_FINAL.csv",
          na = "NA",
          row.names = FALSE)
