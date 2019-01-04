setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load relevant active clinical trial files
#----------------------------------------------
OnCore_Biomarker_Report <- 
  read.csv(file = paste("ClinicalTrialMatching/Biomarker_Report_LongFormat_", OnCore_Biomarker_Report_timestamp, ".csv", sep=""),
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
for (row_No in 1:nrow(OnCore_Biomarker_Report)) {
  if (is.na(OnCore_Biomarker_Report$Biomarker_Comment[row_No])) {
    OnCore_Biomarker_Report$Biomarker_Comment[row_No] <- "None"
  }
}

## STAMP entries extracted based on matching of gene name and disease group
int_file_01 = paste(getwd(), "/ClinicalTrialMatching/OnCore_Biomarker_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    filterName_initial, sep="")
DF_Output_OnCore_Biomarker <- read.csv(file = paste(int_file_01,".csv",sep=""), 
                                       header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")
ncol_No2Keep = as.numeric(which(colnames(DF_Output_OnCore_Biomarker) == "primaryTumorSite.category"))

# ## Remove duplicated clinical trials (extended melt procedure)
# row_total <- as.numeric(ncol(DF_Output_OnCore_Biomarker))
# for (row_No in 1:nrow(DF_Output_OnCore_Biomarker)) {
#   DF <- DF_Output_OnCore_Biomarker[row_No,(ncol_No2Keep +1):row_total]
#   DF <- DF[!duplicated(as.list(DF))]
#   nexist = as.numeric(ncol(DF))
#   ntotal = row_total - ncol_No2Keep
#   DF[,nexist+1:(ntotal - nexist)] <- NA
#   
#   DF_Output_OnCore_Biomarker[row_No,(ncol_No2Keep +1):row_total] <- DF
# }
# remove(row_total,ntotal,nexist,DF)

# Convert matched clinical trial to long format
#----------------------------------------------
colname_keep <- colnames(DF_Output_OnCore_Biomarker)[1:ncol_No2Keep]

# Merge individual matched description
col_extract <- c("Disease.Group.category" ,"OnCore.No")
col_start <- as.numeric(which(colnames(DF_Output_OnCore_Biomarker) == "primaryTumorSite.category"))
biomarker_iter = (ncol(DF_Output_OnCore_Biomarker) - col_start) %/% (length(col_extract))

for (n_iter in biomarker_iter:1) {
  col_int_start = col_start +1 +(length(col_extract)*(n_iter -1))
  col_int_end = col_int_start +length(col_extract) -1
  
  print(paste("For n_iter=", n_iter, " col_start=", col_int_start, 
              " and col_end=", col_int_end, sep=""))
  }

################################
# Manual edit
################################
if (isTRUE(filterName_initial == "diseaseFILTER_groupOFF")) {
  
  # For n_iter=3 col_start=78 and col_end=79
  DF_Output_OnCore_Biomarker <- unite(DF_Output_OnCore_Biomarker, 
                                      newcol.3, c(78:79), sep = "_", remove = TRUE)
  # For n_iter=2 col_start=76 and col_end=77
  DF_Output_OnCore_Biomarker <- unite(DF_Output_OnCore_Biomarker, 
                                      newcol.2, c(76:77), sep = "_", remove = TRUE)
  # For n_iter=1 col_start=74 and col_end=75
  DF_Output_OnCore_Biomarker <- unite(DF_Output_OnCore_Biomarker, 
                                      newcol.1, c(74:75), sep = "_", remove = TRUE)
} else if (isTRUE(filterName_initial == "diseaseFILTER_groupON")) {
  
  # For n_iter=1 col_start=74 and col_end=75
  DF_Output_OnCore_Biomarker <- unite(DF_Output_OnCore_Biomarker, 
                                      newcol.1, c(74:75), sep = "_", remove = TRUE)
}

# Convert matched description to long format
DF_Output_OnCore_Biomarker <- melt(DF_Output_OnCore_Biomarker, 
                                   id.vars=colname_keep)
DF_Output_OnCore_Biomarker <- 
  DF_Output_OnCore_Biomarker[which(DF_Output_OnCore_Biomarker$value != "NA_NA"), 
                              c(colname_keep,"value")]
DF_Output_OnCore_Biomarker$value <- as.character(DF_Output_OnCore_Biomarker$value)

# Strip matched description into components
DF_Output_OnCore_Biomarker <- cSplit(DF_Output_OnCore_Biomarker, "value", "_", 
                                      stripWhite = TRUE, type.convert="as.character")

# Input associated clinical trial info from OnCore_Biomarker_Report
#----------------------------------------------
# Format structure of DF_Output_OnCore_Biomarker
ncol_original <- as.numeric(ncol((DF_Output_OnCore_Biomarker)))
colnames(DF_Output_OnCore_Biomarker)[[ncol_original -1]] <- "Disease.Group.category"
colnames(DF_Output_OnCore_Biomarker)[[ncol_original]] <- "value"

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
  site_id <- DF_Output_OnCore_Biomarker$Disease.Group.category[match_num]
  var_id <- DF_Output_OnCore_Biomarker$var.anno[match_num]
  
  match_info <- OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$OnCore.No == match_id &
                        OnCore_Biomarker_Report$Biomarker_GeneName == gene_id &
                          OnCore_Biomarker_Report$Disease.Group.category == site_id &
                          OnCore_Biomarker_Report$Biomarker_Condition == var_id),]
  
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
    DF_Output_Pre <- cbind(DF_Output_OnCore_Biomarker[match_num,1:ncol_original],match_info[1,])

  } else {
    DF_Output_Pre <- cbind(DF_Output_OnCore_Biomarker[match_num,1:ncol_original],match_info)
  }
  
  DF_Output_Biomarker_Matched <- rbind.fill(DF_Output_Biomarker_Matched, DF_Output_Pre)
}

DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[rowSums(is.na(DF_Output_Biomarker_Matched)) != ncol(DF_Output_Biomarker_Matched),]

# which(DF_Output_Biomarker_Matched$value != DF_Output_Biomarker_Matched$OnCore.No)
DF_Output_Biomarker_Matched <- 
  DF_Output_Biomarker_Matched[,!colnames(DF_Output_Biomarker_Matched) == "value"]

remove(DF_Output_OnCore_Biomarker,DF_Output_Pre,match_info,OnCore_Biomarker_Report,
       gene_id,match_id,match_num,ncol_original,site_id,var_id,col_new,Criteria_No,
       ncol_original_match,row_end,row_end_int,row_start_int,row_No)

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

################################
# Manual edit
################################
# For n_iter=2 col_start=95 and col_end=99
DF_Output_Biomarker_Matched <- unite(DF_Output_Biomarker_Matched,
                                     newcol.3, c(95:99), sep = "_", remove = TRUE)
# For n_iter=1 col_start=90 and col_end=94
DF_Output_Biomarker_Matched <- unite(DF_Output_Biomarker_Matched,
                                     newcol.2, c(90:94), sep = "_", remove = TRUE)

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
remove(biomarker_iter,col_extract,col_int_start,col_int_end,col_start,colname_keep,n_iter,ncol_No2Keep)

## Match Disease.Site from clinical trial
#----------------------------------------------

# ################################
# # Manual Confirmation of 
# ################################
# Disease.Site.trial <- unique(DF_Output_Biomarker_Matched$Disease.Site)
# primaryTumorSite.patient <- sort(unique(DF_Output_Biomarker_Matched$smpl.primaryTumorSite))
# 
# # Strip site_id into components
# site_id_unlist <- list()
# site_id_component <- strsplit(Disease.Site.trial, split = ",")
# for (elem_No in 1:length(site_id_component)) {
#   for (elem_No_sub in 1:length(site_id_component[[elem_No]])) {
#     site_id_component_pre <- strsplit(site_id_component[[elem_No]][elem_No_sub], split = "and")
#     site_id_unlist <- append(site_id_unlist, site_id_component_pre)
#   }
# }
# 
# site_id_component <- list()
# for (elem_No in 1:length(unlist(site_id_unlist))) {
#   # Remove leading and/or trailing whitespace from character strings
#   site_id_component <- append(site_id_component, trimws(unlist(site_id_unlist)[elem_No]))
# }
# site_id_component <- unlist(site_id_component)
# 
# print(site_id_component)
# print(primaryTumorSite.patient)
# 
# print("Disease site matches with primaryTumorSite(s):")
# for (elem_No in 1:length(site_id_component)) {
#   if (length(which(grepl(site_id_component[elem_No], primaryTumorSite.patient, ignore.case = TRUE))) > 0) {
#     cat(paste(site_id_component[elem_No], " (disease site) matches with following primaryTumorSite(s): ", sep=""), "\n")
#     print(primaryTumorSite.patient[which(grepl(site_id_component[elem_No], primaryTumorSite.patient, ignore.case = TRUE))])
#     cat("\n")  
#   }
# }
# remove(Disease.Site.trial,primaryTumorSite.patient,site_id_unlist,site_id_component,
#        elem_No,elem_No_sub,site_id_component_pre)
  
DF_Output_Biomarker_Matched_FINAL <- data.frame(matrix(NA, ncol = ncol(DF_Output_Biomarker_Matched)))
colnames(DF_Output_Biomarker_Matched_FINAL) <- colnames(DF_Output_Biomarker_Matched) 

patient.list <- unique(DF_Output_Biomarker_Matched$sys.uniqueId)

for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  OnCore.patient <- sort(unique(DF_Output_Biomarker_Matched$OnCore.No
                                [DF_Output_Biomarker_Matched$sys.uniqueId == patient.list[patient_num]]))
  
  for (trial_num in 1:length(OnCore.patient)) {
    OnCore_id <- OnCore.patient[trial_num]
    
    DF_patient <- DF_Output_Biomarker_Matched[which(DF_Output_Biomarker_Matched$sys.uniqueId == patient_id &
                                                      DF_Output_Biomarker_Matched$OnCore.No == OnCore_id),]
    
    Disease.Site.trial <- unique(DF_patient$Disease.Site)
    primaryTumorSite.patient <- unique(DF_patient$smpl.primaryTumorSite)
    # smpl.specimenSite is not an appropriate proxy because it can be the site of a metastasis
    disease_site_match <- NA
    
    ## Modify based on whether disease.site_FILTER == TRUE
    if (isTRUE(disease.site_FILTER & disease.group_FILTER)) {
      
      # Remove corresponding trial INFO from "DF_patient" if no match
      if (!("Any Site" %in% Disease.Site.trial)) {
        
        for (site_num in 1:length(Disease.Site.trial)) {
          site_id <- Disease.Site.trial[site_num]
          
          # Strip site_id into components
          site_id_unlist <- list()
          site_id_component <- strsplit(site_id, split = ",")
          for (elem_No in 1:listLen(site_id_component)) {
            site_id_component_pre <- strsplit(site_id_component[[1]][elem_No], split = "and")
            site_id_unlist <- append(site_id_unlist, site_id_component_pre)
          }
          
          site_id_component <- list()
          for (elem_No in 1:length(unlist(site_id_unlist))) {
            # Remove leading and/or trailing whitespace from character strings
            site_id_component <- append(site_id_component, trimws(unlist(site_id_unlist)[elem_No]))
          }
          site_id_component <- unlist(site_id_component)
          
          for (elem_No in 1:length(site_id_component)) {
            if (isTRUE(grepl(site_id_component[elem_No], primaryTumorSite.patient, ignore.case = TRUE))) {
              disease_site_match <- append(disease_site_match, "TRUE")
            } else {
              disease_site_match <- append(disease_site_match, "FALSE")
            }
          }
          
          if ("TRUE" %in% disease_site_match) {
            disease_site_match <- TRUE
          } else {
            disease_site_match <- FALSE
          }
          
          if (isTRUE(disease_site_match == FALSE)) {
            DF_patient <- DF_patient[DF_patient$Disease.Site != site_id,]
          }
        }
        remove(site_num,site_id,site_id_unlist,site_id_component_pre,site_id_component,elem_No)
      }
    }
    
    if (nrow(DF_patient) > 0) {
      DF_Output_Biomarker_Matched_FINAL <- rbind(DF_Output_Biomarker_Matched_FINAL, DF_patient)
    }
  }
}

DF_Output_Biomarker_Matched_FINAL <- 
  DF_Output_Biomarker_Matched_FINAL[rowSums(is.na(DF_Output_Biomarker_Matched_FINAL)) != ncol(DF_Output_Biomarker_Matched_FINAL),]

remove(Disease.Site.trial,OnCore_id,OnCore.patient,patient_id,patient_num,patient.list,disease_site_match,
       primaryTumorSite.patient,trial_num,DF_patient,DF_Output_Biomarker_Matched)

# ## Data Review
# #----------------------------------------------
# print(paste("Total number of unique OnCore.No ", length(unique(DF_Output_Biomarker_Matched_FINAL$OnCore.No)), sep=""))
# print(paste("Total number of unique disease group categories: ", length(sort(unique(DF_Output_Biomarker_Matched_FINAL$Disease.Group.category))), sep=""))
# sort(unique(DF_Output_Biomarker_Matched_FINAL$Disease.Group.category))
# print(paste("Total number of unique disease sites: ", length(sort(unique(DF_Output_Biomarker_Matched_FINAL$Disease.Site))), sep=""))
# sort(unique(DF_Output_Biomarker_Matched_FINAL$Disease.Site))
# print(paste("Total number of unique biomarker genes: ", length(sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_GeneName))), sep=""))
# sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_GeneName))
# print(paste("Total number of unique biomarker conditions: ", length(sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_Condition))), sep=""))
# sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_Condition))
# print(paste("Total number of unique Biomarker details: ", length(sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_Detail))), sep=""))
# sort(unique(DF_Output_Biomarker_Matched_FINAL$Biomarker_Detail))

# Match with Biomarker_Detail
# Options: "All mutations accepted", "Amino acid change", "Translocation partner gene"
#----------------------------------------------
# Extract Amino acid change criteria
DF_Output_Biomarker_Matched_FINAL_edit <- 
  DF_Output_Biomarker_Matched_FINAL[DF_Output_Biomarker_Matched_FINAL$Biomarker_Detail != "All mutations accepted",]
DF_Output_Biomarker_Matched_FINAL <- 
  DF_Output_Biomarker_Matched_FINAL[DF_Output_Biomarker_Matched_FINAL$Biomarker_Detail == "All mutations accepted",]

################################
# Manual parse
################################
# Remove corresponding trial INFO if no match in "Amino acid change"
for (row_No in 1:nrow(DF_Output_Biomarker_Matched_FINAL_edit)) {
  
  # sort(unique(DF_Output_Biomarker_Matched_FINAL_edit$Biomarker_Detail))
  # gsub function removes asterisk 
  if (isTRUE(grepl(gsub("(^[[:alpha:]]{3}[[:digit:]]{,4})(.*)", "p.\\1", DF_Output_Biomarker_Matched_FINAL_edit$Biomarker_Detail[row_No]), 
            DF_Output_Biomarker_Matched_FINAL_edit$smpl.hgvsProtein[row_No], ignore.case = TRUE) == FALSE)) {
    DF_Output_Biomarker_Matched_FINAL_edit <- 
      DF_Output_Biomarker_Matched_FINAL_edit[DF_Output_Biomarker_Matched_FINAL_edit$smpl.hgvsProtein != 
                                               DF_Output_Biomarker_Matched_FINAL_edit$smpl.hgvsProtein[row_No],]
  }
}

DF_Output_Biomarker_Matched_FINAL <- rbind(DF_Output_Biomarker_Matched_FINAL, DF_Output_Biomarker_Matched_FINAL_edit)
remove(DF_Output_Biomarker_Matched_FINAL_edit,row_No)

## Write to local computer
#----------------------------------------------
write.csv(DF_Output_Biomarker_Matched_FINAL, file = paste(int_file_01, filterName_int, "_FINAL.csv", sep=""),
          na = "NA", row.names = FALSE)

# Delete intermediate file
if (isTRUE(deleteIntermediateFile)) {
  if (file.exists(paste(int_file_01,".csv",sep=""))){file.remove(paste(int_file_01,".csv",sep=""))}
}

remove(int_file_01)
