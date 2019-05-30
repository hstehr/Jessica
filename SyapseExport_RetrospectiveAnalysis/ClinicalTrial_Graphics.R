# Generate graphics illustrating matching rate

plotTitle_label <- "Trial Matching"
plot_filename <- paste("_", filterName, sep="")

## Load dataframes
if (isTRUE(Internal_match)) {
  plotTitle_label <- paste("OnCore ",plotTitle_label,sep="")
  plot_filename <- paste("OnCore_", OnCore_Biomarker_Report_timestamp, plot_filename, sep="")
  
  DF_Output_OnCore_Biomarker <-
    read.csv(file = paste(tempdir,"OnCore_SNVIndel_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                          groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
}

if (isTRUE(NCI_match)) {
  plotTitle_label <- paste("NCI-MATCH ",plotTitle_label,sep="")
  plot_filename <- paste("_NCI-MATCH_", Patient_Variant_Report_timestamp, plot_filename, sep="")
  
  DF_Output_Patient_Variant <- 
    read.csv(file = paste(tempdir,"NCI_SNVIndel_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                          dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  DF_Output_Patient_NonHotspot <- 
    read.csv(file = paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                          dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
}

plotTitle_label <- gsub("NCI-MATCH OnCore", "NCI-MATCH & OnCore", plotTitle_label)
plot_filename <- paste("Total_MatchDistribution", plot_filename, sep="")

## ILLUSTRATION PLOTs: Matching rate
#---------------------------------------------- 
# Determine which patient_id-Arm_Name combo are duplicated
# Compare results from Inclusion_Variants and Inclusion_NonHotspot
if (isTRUE(NCI_match & nrow(DF_Output_Patient_Variant) > 0 & nrow(DF_Output_Patient_NonHotspot) > 0)) {
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
cat(paste(nrow(Trial_Count), " total patients plotted", sep=""),"\n","\n")

if (isTRUE(Internal_match)) {
Trial_Count <- cbind(Trial_Count, 
                     data.frame(No.Internal.Trials = NA, stringsAsFactors = FALSE))
}

if (isTRUE(NCI_match)) {
Trial_Count <- cbind(Trial_Count, 
                     data.frame(No.NCI.Trials = NA, 
                                No.NCI.NonHotspot = NA, stringsAsFactors = FALSE))
}

for (row_No in 1:nrow(Trial_Count)) {
  pt_id <- Trial_Count$patient_id[row_No]
  
  if (isTRUE(Internal_match)) {
    if (nrow(DF_Output_OnCore_Biomarker) > 0) {
      Trial_Count$No.Internal.Trials[row_No] <- 
        as.numeric(length(unique(DF_Output_OnCore_Biomarker$OnCore.No
                                 [which(DF_Output_OnCore_Biomarker$PatientID == pt_id)])))
    } else {
      Trial_Count$No.Internal.Trials[row_No] <- 0
    }
  }
  
  if (isTRUE(NCI_match)) {
    if (nrow(DF_Output_Patient_Variant) > 0) {
      Trial_Count$No.NCI.Trials[row_No] <- 
        as.numeric(length(unique(DF_Output_Patient_Variant$Arm_Name
                                 [which(DF_Output_Patient_Variant$PatientID == pt_id)])))
    } else {
      Trial_Count$No.NCI.Trials[row_No] <- 0
    }
    
    if (nrow(DF_Output_Patient_NonHotspot) > 0) {
      Trial_Count$No.NCI.NonHotspot[row_No] <- 
        as.numeric(length(unique(DF_Output_Patient_NonHotspot$Arm_Name
                                 [which(DF_Output_Patient_NonHotspot$PatientID == pt_id)])))
    } else {
      Trial_Count$No.NCI.NonHotspot[row_No] <- 0
    }
  
  # Remove duplicated patient_id-Arm_Name combo for NCI-MATCH trials (Variant and NonHotspot matches)
  if (nrow(DF_Output_Patient_Variant) > 0 & nrow(DF_Output_Patient_NonHotspot) > 0) {
    if (Trial_Count$patient_id[row_No] %in% DF_NCI_Matched_Dup) {
      Trial_Count$No.NCI.NonHotspot[row_No] <- Trial_Count$No.NCI.NonHotspot[row_No] -1
    }
  }
  }
}
Trial_Count$No.Total.Trials <- rowSums(Trial_Count[,3:ncol(Trial_Count)])

# Calculate y-axis of plots 
max_list <- length(which(Trial_Count$No.Total.Trials == 0))

if (isTRUE(Internal_match)) {
  if (nrow(DF_Output_OnCore_Biomarker) > 0) {
  max_list <- append(max_list, length(which(Trial_Count$No.Internal.Trials == 0)))
  }
}

if (isTRUE(NCI_match)) {
  if (nrow(DF_Output_Patient_Variant) > 0) {
    max_list <- append(max_list, length(which(Trial_Count$No.NCI.Trials == 0)))
  }
  if (nrow(DF_Output_Patient_NonHotspot) > 0) {
    max_list <- append(max_list, length(which(Trial_Count$No.NCI.NonHotspot == 0)))
  }
}

y_max = ceiling(max(max_list)/100) * 100
if (isTRUE(y_max < 1000)) {
  if (isTRUE(y_max <= 20)) {y_increment = 1
  } else if (isTRUE(y_max <= 30)) {y_increment = 2
  } else if (isTRUE(y_max <= 100)) {y_increment = 5
  } else if (isTRUE(y_max <= 200)) {y_increment = 10
  } else if (isTRUE(y_max <= 500)) {y_increment = 50
  } else if (isTRUE(y_max <= 1000)) {y_increment = 100
  } else {y_increment = 10
  }
} else {y_increment = 250
}

## FUNCTION: save histogram plot in local computer
plot_Trial_Count <- function(DF, trialColumn, fileName_pre, plotTitle, y_max) {
  
  # Plot parameters 
  total_matched <- as.numeric(length(which(DF[[trialColumn]] != 0)))
  x_max = as.numeric(max(DF[[trialColumn]]))
  if (x_max > 2) { image_width = x_max +2
  } else { image_width = x_max +3
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF$Tabulate <- DF[[trialColumn]]
  
  DF_tabulate <- data.frame(DF %>% group_by(Tabulate) %>% tally())
  colnames(DF_tabulate) <- c("No.Trials","No.Orders")
  
  NoTrial.missing <- setdiff(seq(0,x_max), unique(DF_tabulate$No.Trials))
  if (length(NoTrial.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                              data.frame(No.Trials=NoTrial.missing,No.Orders=0, stringsAsFactors = FALSE))
  }
  DF_tabulate$No.Trials <- factor(DF_tabulate$No.Trials,
                                             levels=c(seq(0,x_max,1)))
  DF_tabulate$No.Orders <- as.numeric(DF_tabulate$No.Orders)
  
  # Specify palette 
  x_rows = nrow(DF_tabulate)
  custom.palette = get(paste("custom.hues.",x_rows,sep=""))
  
  plot <- ggplot(DF_tabulate, aes(x=No.Trials, y=No.Orders, fill=No.Trials)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=No.Orders), vjust=-0.75, size=4.0)+
    
    labs(title = plotTitle,
         subtitle = paste("Match rate: ",format(round(as.numeric(100*total_matched/nrow(DF)),1),nsmall = 1), "% (n = ", 
                          total_matched, "/", nrow(DF), ")", sep="")) +
    
    scale_x_discrete(name="No. Clinical Trials Matched") +
    scale_y_continuous(name="No. Test Orders", breaks = seq(0,ceiling(y_max/y_increment)*y_increment,y_increment), 
                       limits=c(0,y_max)) +
    
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
          plot.subtitle = element_text(hjust=1, face="bold",size=12),
          legend.text=element_text(size=10),
          axis.text=element_text(size=12,face="bold"), 
          axis.title=element_text(size=12,face="bold"),
          
          legend.position="none") +
    
    scale_fill_manual(values = custom.palette)
  
  tiff(filename = paste(outdir, "/", fileName_pre, "_", Sys.Date(), ".tiff", sep=""),
       width = image_width, height = 7, units = "in", res = 350)
  print(plot)
  dev.off()
}

plot_Trial_Count(DF = Trial_Count, 
                 y_max = y_max,
                 trialColumn = "No.Total.Trials",
                 fileName_pre = plot_filename,
                 plotTitle = plotTitle_label) 

if (isTRUE(Internal_match)) {
  if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
    plot_Trial_Count(DF = Trial_Count, 
                     y_max = y_max,
                     trialColumn = "No.Internal.Trials", 
                     fileName_pre = paste("OnCore_MatchDistribution_", OnCore_Biomarker_Report_timestamp, "_", 
                                          groupName,siteName, "_", pathoName, "_",  ageName, sep=""), 
                     plotTitle = "OnCore Clinical Trial Matching") 
  }
}

if (isTRUE(NCI_match)) {
  for (row_No in 1:nrow(Trial_Count)) {
    Trial_Count$No.NCI.Total[row_No] <- 
      as.numeric(Trial_Count$No.NCI.Trials[row_No] + Trial_Count$No.NCI.NonHotspot[row_No])
  }
  
  plot_Trial_Count(DF = Trial_Count, y_max = y_max,
                   trialColumn = "No.NCI.Total", 
                   fileName_pre = paste("MATCH_Total_MatchDistribution_",
                                        Patient_Variant_Report_timestamp, "_", 
                                        dxName, "_", pathoName, "_",  ageName, sep=""), 
                   plotTitle = "NCI-MATCH Trial Matching") 
}

if (isTRUE(exists("DF_NCI_Matched_Dup"))) {remove(DF_NCI_Matched_Dup)}
if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {remove(DF_Output_OnCore_Biomarker)}
if (isTRUE(exists("DF_Output_Patient_Variant"))) {remove(DF_Output_Patient_Variant)}
if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {remove(DF_Output_Patient_NonHotspot)}

remove(Trial_Count,pt_id,row_No,plot_Trial_Count,max_list)
