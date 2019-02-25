# Generate graphics illustrating matching rate

## ILLUSTRATION PLOTs: Matching rate
#---------------------------------------------- 
# Determine which patient_id-Arm_Name combo are duplicated
# Compare results from Inclusion_Variants and Inclusion_NonHotspot
if (isTRUE(nrow(DF_Output_Patient_Variant) > 0 & nrow(DF_Output_Patient_NonHotspot) > 0)) {
  DF_NCI_Matched_Dup <- rbind(unique(data.frame(id = DF_Output_Patient_Variant$PatientID,
                                                arm = DF_Output_Patient_Variant$Arm_Name,
                                                trial = "NCI.Inclusion.Variants", stringsAsFactors = FALSE)),
                              unique(data.frame(id = DF_Output_Patient_NonHotspot$PatientID,
                                                arm = DF_Output_Patient_NonHotspot$Arm_Name,
                                                trial = "NCI.Inclusion.NonHotspot", stringsAsFactors = FALSE)))
  DF_NCI_Matched_Dup <- DF_NCI_Matched_Dup$id[which(duplicated(DF_NCI_Matched_Dup[,1:2]))]
}

# Tally No. trial matches from each source
Trial_Count <- data.frame(patient_id = STAMP_DF$PatientID,
                          patient_age = as.numeric(STAMP_DF$PatientAge),
                          stringsAsFactors = FALSE)

# Remove duplicate entries 
Trial_Count <- unique(Trial_Count)

# Subset only ADULT patients (18+yo)
if (isTRUE(pathogenic_FILTER)) {
  Trial_Count <- Trial_Count[Trial_Count$patient_age >= 18,]
}

Trial_Count <- cbind(Trial_Count, 
                     data.frame(No.Internal.Trials = NA, No.NCI.Trials = NA, 
                                No.NCI.NonHotspot = NA, stringsAsFactors = FALSE))

for (row_No in 1:nrow(Trial_Count)) {
  pt_id <- Trial_Count$patient_id[row_No]
  
  if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
    Trial_Count$No.Internal.Trials[row_No] <- 
      as.numeric(length(unique(DF_Output_OnCore_Biomarker$OnCore.No
                               [which(DF_Output_OnCore_Biomarker$PatientID == pt_id)])))  
  } else {
    Trial_Count$No.Internal.Trials[row_No] <- 0
  }
  
  if (isTRUE(exists("DF_Output_Patient_Variant"))) {
    Trial_Count$No.NCI.Trials[row_No] <- 
      as.numeric(length(unique(DF_Output_Patient_Variant$Arm_Name
                               [which(DF_Output_Patient_Variant$PatientID == pt_id)])))
  } else {
    Trial_Count$No.NCI.Trials[row_No] <- 0
  }
  
  if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
    Trial_Count$No.NCI.NonHotspot[row_No] <- 
      as.numeric(length(unique(DF_Output_Patient_NonHotspot$Arm_Name
                               [which(DF_Output_Patient_NonHotspot$PatientID == pt_id)])))
  } else {
    Trial_Count$No.NCI.NonHotspot[row_No] <- 0
  }
  
  # Remove duplicated patient_id-Arm_Name combo for NCI-MATCH trials (Variant and NonHotspot matches)
  if (isTRUE(nrow(DF_Output_Patient_Variant) > 0 & nrow(DF_Output_Patient_NonHotspot) > 0)) {
    if (Trial_Count$patient_id[row_No] %in% DF_NCI_Matched_Dup) {
      Trial_Count$No.NCI.NonHotspot[row_No] <- Trial_Count$No.NCI.NonHotspot[row_No] -1
    }
  }
}
Trial_Count$No.Total.Trials <- rowSums(Trial_Count[,3:5])

# Calculate y-axis of plots 
max_list <- length(which(Trial_Count$No.Total.Trials == 0))

if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
  max_list <- append(max_list, length(which(Trial_Count$No.Internal.Trials == 0)))
}
if (isTRUE(exists("DF_Output_Patient_Variant"))) {
  max_list <- append(max_list, length(which(Trial_Count$No.NCI.Trials == 0)))
}
if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
  max_list <- append(max_list, length(which(Trial_Count$No.NCI.NonHotspot == 0)))
}

y_max = round(max(max_list)/100, digits = 1) * 100

## FUNCTION: save histogram plot in local computer
plot_Trial_Count <- function(DF, trialColumn, fileName_pre, plotTitle, y_max) {
  
  # Plot parameters 
  total_matched <- as.numeric(length(which(DF[[trialColumn]] != 0)))
  x_max = as.numeric(max(DF[[trialColumn]]))
  if (x_max > 2) { image_width = x_max +2
  } else { image_width = x_max +3
  }
  
  tiff(filename = paste(outdir, fileName_pre, "_", Sys.Date(), ".tiff", sep=""),
       width = image_width, height = 7, units = "in", res = 150)
  
  plot <- ggplot(DF, aes(x = DF[[trialColumn]])) +
    geom_histogram(bins = x_max +1, col="gray") +
    
    geom_text(stat='count', aes(label=..count..), vjust=-1) +
    labs(title = plotTitle,
         subtitle = paste("N = ", total_matched, " / ", nrow(DF), " adult patients matched", sep="")) +
    
    scale_x_continuous(name="No. Clinical Trials", breaks = seq(0, x_max, 1)) +
    scale_y_continuous(name="No. Individuals", limits = c(0,y_max), breaks = seq(0,y_max, by = 100)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=12),
          plot.subtitle = element_text(hjust=1, face="bold",size=12),
          legend.text=element_text(size=10),
          axis.text=element_text(size=15), 
          axis.title=element_text(size=15,face="bold"))
  
  print(plot)
  dev.off()
}

plot_Trial_Count(DF = Trial_Count, y_max = y_max,
                 trialColumn = "No.Total.Trials", 
                 fileName_pre = paste("Total_MatchDistribution_Internal_", OnCore_Biomarker_Report_timestamp, 
                                      "_NCI-MATCH_", Patient_Variant_Report_timestamp, "_", filterName, sep=""), 
                 plotTitle = "NCI-MATCH & Stanford Internal Clinical Trial Matching") 

if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
  plot_Trial_Count(DF = Trial_Count, y_max = y_max,
                   trialColumn = "No.Internal.Trials", 
                   fileName_pre = paste("OnCore_Biomarker_MatchDistribution_", OnCore_Biomarker_Report_timestamp, "_", 
                                        groupName,siteName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "Stanford Internal Clinical Trial Matching") 
}

if (isTRUE(exists("DF_Output_Patient_Variant"))) {
  plot_Trial_Count(DF = Trial_Count, y_max = y_max,
                   trialColumn = "No.NCI.Trials", 
                   fileName_pre = paste("Patient_Variant_MatchDistribution_",
                                        Patient_Variant_Report_timestamp, "_", 
                                        dxName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "NCI-MATCH Trial Matching (Variants)") 
}

if (isTRUE(exists("DF_Output_Patient_NonHotspot")) & nrow(DF_Output_Patient_NonHotspot) > 0) {
  plot_Trial_Count(DF = Trial_Count, y_max = y_max,
                   trialColumn = "No.NCI.NonHotspot", 
                   fileName_pre = paste("Patient_NonHotspot_MatchDistribution_",
                                        Patient_Variant_Report_timestamp, "_", 
                                        dxName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "NCI-MATCH Trial Matching (NonHotspot)") 
}

if (isTRUE(exists("DF_NCI_Matched_Dup"))) { remove(DF_NCI_Matched_Dup)}
remove(Trial_Count,pt_id,row_No,plot_Trial_Count,max_list)
