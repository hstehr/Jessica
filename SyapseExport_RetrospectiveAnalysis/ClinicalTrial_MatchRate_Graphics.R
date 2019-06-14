# Generate graphics illustrating match rate

Output_FINAL <- data.frame()
plotTitle_label <- "Trial Matching"
plot_filename <- paste("_", filterName, sep="")
ymax <- c()

colnames_extract <- c("PatientID","Trial.No","Match_Detail")
colnames_extract.OnCore <- c("PatientID","OnCore.No","Match_Detail")
colnames_extract.MATCH <- c("PatientID","Arm_Name","Match_Detail")

if (isTRUE(Internal_match)) {
  plotTitle_label <- paste("OnCore ",plotTitle_label,sep="")
  plot_filename <- paste("OnCore_", OnCore_Biomarker_Report_timestamp, plot_filename, sep="")
  
  # Import candidate matches and specify empty dataframes
  #---------------------------------------------- 
  # SNV Indels
  file_name = paste(tempdir,"OnCore_SNVIndel_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_OnCore) > 0) {continue_SNVIndel <- as.logical("TRUE")
    } else {continue_SNVIndel <- as.logical("FALSE")}
    
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  # CNVs
  file_name = paste(tempdir,"OnCore_CNV_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_OnCore) > 0) {continue_CNV <- as.logical("TRUE")
    } else {continue_CNV <- as.logical("FALSE")}
    
  } else {continue_CNV <- as.logical("FALSE")}
  
  # Fusions
  file_name = paste(tempdir,"OnCore_Fusion_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_Fusion_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_Fusion_OnCore) > 0) {continue_Fusion <- as.logical("TRUE")
    } else {continue_Fusion <- as.logical("FALSE")}
    
  } else {continue_Fusion <- as.logical("FALSE")}
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_OnCore_FINAL <- data.frame(matrix(ncol = length(colnames_extract.OnCore), nrow = 0))
  colnames(Output_OnCore_FINAL) <- colnames_extract.OnCore
  
  if (isTRUE(continue_SNVIndel)) {
    Output_SNVIndel_OnCore$Match_Detail <- "SNVIndel_OnCore"
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_SNVIndel_OnCore[,colnames_extract.OnCore])
    
    remove(Output_SNVIndel_OnCore,continue_SNVIndel)
  }
  
  if (isTRUE(continue_CNV)) {
    Output_CNV_OnCore$Match_Detail <- "CNV_OnCore"
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_CNV_OnCore[,colnames_extract.OnCore])
    
    remove(Output_CNV_OnCore,continue_CNV)
  }
  
  if (isTRUE(continue_Fusion)) {
    Output_Fusion_OnCore$Match_Detail <- "Fusion_OnCore"
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_Fusion_OnCore[,colnames_extract.OnCore])
    
    remove(Output_Fusion_OnCore,continue_Fusion)
  }
  
  Output_OnCore_FINAL <- unique(Output_OnCore_FINAL[,])
  colnames(Output_OnCore_FINAL) <- colnames_extract
  
  Output_OnCore_FINAL$Trial.Type <- "OnCore"
  Output_FINAL <- rbind(Output_FINAL,Output_OnCore_FINAL)
  remove(Output_OnCore_FINAL)
}

if (isTRUE(NCI_match)) {
  plotTitle_label <- paste("NCI-MATCH ",plotTitle_label,sep="")
  plot_filename <- paste("_NCI-MATCH_", Patient_Variant_Report_timestamp, plot_filename, sep="")
  
  # Import candidate matches and specify empty dataframes
  #---------------------------------------------- 
  # SNVIndel_Variant
  file_name = paste(tempdir,"NCI_SNVIndel_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_Variant) > 0) {continue_SNVIndel_Variant <- as.logical("TRUE")
    } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
    
  } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
  
  # CNV_Variant
  file_name = paste(tempdir,"NCIMatch_CNV_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_Variant) > 0) {continue_CNV_Variant <- as.logical("TRUE")
    } else {continue_CNV_Variant <- as.logical("FALSE")}
    
  } else {continue_CNV_Variant <- as.logical("FALSE")}
  
  # Fusion_Variant
  file_name = paste(tempdir,"NCI_Fusion_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_Fusion_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_Fusion_Variant) > 0) {continue_Fusion_Variant <- as.logical("TRUE")
    } else {continue_Fusion_Variant <- as.logical("FALSE")}
    
  } else {continue_Fusion_Variant <- as.logical("FALSE")}
  
  # SNVIndel_NonHotspot
  file_name = paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_NonHotspot <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_NonHotspot) > 0) {continue_SNVIndel_NonHotspot <- as.logical("TRUE")
    } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
    
  } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
  
  # CNV_NonHotspot
  file_name = paste(tempdir,"NCI_CNV_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_NonHotspot <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_NonHotspot) > 0) {continue_CNV_NonHotspot <- as.logical("TRUE")
    } else {continue_CNV_NonHotspot <- as.logical("FALSE")}
    
  } else {continue_CNV_NonHotspot <- as.logical("FALSE")}
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_NCI_FINAL <- data.frame(matrix(ncol = length(colnames_extract.MATCH), nrow = 0))
  colnames(Output_NCI_FINAL) <- colnames_extract.MATCH
  
  if (isTRUE(continue_SNVIndel_Variant)) {
    Output_SNVIndel_Variant$Match_Detail <- "SNVIndel_Variant"
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_Variant[,colnames_extract.MATCH])
    
    remove(Output_SNVIndel_Variant,continue_SNVIndel_Variant)
  }
  
  if (isTRUE(continue_CNV_Variant)) {
    Output_CNV_Variant$Match_Detail <- "CNV_Variant"
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_Variant[,colnames_extract.MATCH])
    
    remove(Output_CNV_Variant,continue_CNV_Variant)
  }
  
  if (isTRUE(continue_Fusion_Variant)) {
    Output_Fusion_Variant$Match_Detail <- "Fusion_Variant"
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_Fusion_Variant[,colnames_extract.MATCH])
    
    remove(Output_Fusion_Variant,continue_Fusion_Variant)
  }
  
  if (isTRUE(continue_SNVIndel_NonHotspot)) {
    Output_SNVIndel_NonHotspot$Match_Detail <- "SNVIndel_NonHotspot"
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_NonHotspot[,colnames_extract.MATCH])
    
    remove(Output_SNVIndel_NonHotspot,continue_SNVIndel_NonHotspot)
  }
  
  if (isTRUE(continue_CNV_NonHotspot)) {
    Output_CNV_NonHotspot$Match_Detail <- "CNV_NonHotspot"
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_NonHotspot[,colnames_extract.MATCH])
    
    remove(Output_CNV_NonHotspot,continue_CNV_NonHotspot)
  }
  
  Output_NCI_FINAL <- unique(Output_NCI_FINAL[,])
  colnames(Output_NCI_FINAL) <- colnames_extract
  
  Output_NCI_FINAL$Trial.Type <- "NCI.MATCH"
  Output_FINAL <- rbind(Output_FINAL,Output_NCI_FINAL)
  remove(Output_NCI_FINAL)
}

plotTitle_label <- gsub("NCI-MATCH OnCore", "NCI-MATCH & OnCore", plotTitle_label)
plot_filename <- paste("Total_MatchDistribution", plot_filename, sep="")

# Remove duplicated patient_id-Arm_Name combinations
Output_FINAL_unique <- unique(Output_FINAL[,c("PatientID","Trial.No","Trial.Type")])

if (isTRUE(Internal_match)) {
  # OnCore: tally trial matches from each source
  #---------------------------------------------- 
  Trial_Count_OnCore <- 
    Output_FINAL_unique[which(Output_FINAL_unique$Trial.Type == "OnCore"),] %>% group_by(PatientID,Trial.Type) %>% tally()
  
  max_ct <- max(Trial_Count_OnCore$n)
  Trial_Count_OnCore_tally <- data.frame(OnCore.Trials.Matched = seq(1,max_ct,1),
                                         No.Orders = NA,
                                         stringsAsFactors = FALSE)
  for (i in 1:max_ct) {
    Trial_Count_OnCore_tally$No.Orders[i] = length(which(Trial_Count_OnCore$n == i))
  }
  Trial_Count_OnCore_tally <- rbind(Trial_Count_OnCore_tally,
                                    data.frame(OnCore.Trials.Matched = 0,
                                               No.Orders = length(patient.list) - sum(Trial_Count_OnCore_tally$No.Orders)))
  
  Trial_Count_OnCore_tally <- Trial_Count_OnCore_tally[order(Trial_Count_OnCore_tally$OnCore.Trials.Matched),]
  
  Trial_Count_OnCore_tally$No.Orders <- as.numeric(Trial_Count_OnCore_tally$No.Orders)
  Trial_Count_OnCore_tally$OnCore.Trials.Matched <- factor(Trial_Count_OnCore_tally$OnCore.Trials.Matched,
                                                           levels = seq(0,max_ct,1))
  
  ymax <- append(ymax, max(Trial_Count_OnCore_tally$No.Orders))
}

if (isTRUE(NCI_match)) {
  # NCI.MATCH: tally trial matches from each source
  #---------------------------------------------- 
  Trial_Count_MATCH <- 
    Output_FINAL_unique[which(Output_FINAL_unique$Trial.Type == "NCI.MATCH"),] %>% group_by(PatientID,Trial.Type) %>% tally()
  
  max_ct <- max(Trial_Count_MATCH$n)
  Trial_Count_MATCH_tally <- data.frame(MATCH.Trials.Matched = seq(1,max_ct,1),
                                        No.Orders = NA,
                                        stringsAsFactors = FALSE)
  for (i in 1:max_ct) {
    Trial_Count_MATCH_tally$No.Orders[i] = length(which(Trial_Count_MATCH$n == i))
  }
  Trial_Count_MATCH_tally <- rbind(Trial_Count_MATCH_tally,
                                   data.frame(MATCH.Trials.Matched = 0,
                                              No.Orders = length(patient.list) - sum(Trial_Count_MATCH_tally$No.Orders)))
  Trial_Count_MATCH_tally <- Trial_Count_MATCH_tally[order(Trial_Count_MATCH_tally$MATCH.Trials.Matched),]
  
  Trial_Count_MATCH_tally$No.Orders <- as.numeric(Trial_Count_MATCH_tally$No.Orders)
  Trial_Count_MATCH_tally$MATCH.Trials.Matched <- factor(Trial_Count_MATCH_tally$MATCH.Trials.Matched,
                                                         levels = seq(0,max_ct,1))
  
  ymax <- append(ymax, max(Trial_Count_MATCH_tally$No.Orders))
}

# Total: tally trial matches from each source
#---------------------------------------------- 
Trial_Count_Total <- 
  unique(Output_FINAL_unique[,c("PatientID","Trial.No")]) %>% group_by(PatientID) %>% tally()

max_ct <- max(Trial_Count_Total$n)
Trial_Count_Total_tally <- data.frame(Total.Trials.Matched = seq(1,max_ct,1),
                                      No.Orders = NA,
                                      stringsAsFactors = FALSE)
for (i in 1:max_ct) {
  Trial_Count_Total_tally$No.Orders[i] = length(which(Trial_Count_Total$n == i))
}
Trial_Count_Total_tally <- rbind(Trial_Count_Total_tally,
                                 data.frame(Total.Trials.Matched = 0,
                                            No.Orders = length(patient.list) - sum(Trial_Count_Total_tally$No.Orders)))
Trial_Count_Total_tally <- Trial_Count_Total_tally[order(Trial_Count_Total_tally$Total.Trials.Matched),]

Trial_Count_Total_tally$No.Orders <- as.numeric(Trial_Count_Total_tally$No.Orders)
Trial_Count_Total_tally$Total.Trials.Matched <- factor(Trial_Count_Total_tally$Total.Trials.Matched,
                                                       levels = seq(0,max_ct,1))

ymax <- append(ymax, max(Trial_Count_Total_tally$No.Orders))

# Generate visualizations 
#---------------------------------------------- 
ymax <- ceiling(max(ymax)/10)*10

plot_Trial_Count(DF = Trial_Count_Total_tally, 
                 ymax = ymax,
                 trialColumn = "Total.Trials.Matched",
                 fileName_pre = plot_filename,
                 plotTitle = plotTitle_label,
                 outdir = outdir) 

remove(Trial_Count_Total,Trial_Count_Total_tally)

if (isTRUE(Internal_match)) {
  plot_Trial_Count(DF = Trial_Count_OnCore_tally, 
                   ymax = ymax,
                   trialColumn = "OnCore.Trials.Matched", 
                   fileName_pre = paste("OnCore_MatchDistribution_", OnCore_Biomarker_Report_timestamp, "_", 
                                        groupName,siteName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "OnCore Clinical Trial Matching",
                   outdir = outdir) 
  
  remove(Trial_Count_OnCore,Trial_Count_OnCore_tally)
}

if (isTRUE(NCI_match)) {
  plot_Trial_Count(DF = Trial_Count_MATCH_tally, 
                   ymax = ymax,
                   trialColumn = "MATCH.Trials.Matched", 
                   fileName_pre = paste("MATCH_Total_MatchDistribution_",
                                        Patient_Variant_Report_timestamp, "_", 
                                        dxName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "NCI-MATCH Trial Matching",
                   outdir = outdir) 
  
  remove(Trial_Count_MATCH,Trial_Count_MATCH_tally)
}

remove(plotTitle_label,plot_filename,colnames_extract,colnames_extract.OnCore,
       colnames_extract.MATCH,file_name,Output_FINAL,Output_FINAL_unique,max_ct,ymax)
