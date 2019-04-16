# Generate notification for positive candidate trials
# PositiveMatch.txt

email.comment <- paste("This email was generated on ", Sys.time(), ".",sep="")

if (isTRUE(Internal_match)) {
  email.comment <- paste("OnCore Biomarker Report updated on ", OnCore_Biomarker_Report_timestamp, ". ",
                         email.comment, sep="")
  
  # Import candidate matches
  #---------------------------------------------- 
  Output_SNVIndel_OnCore <- 
    read.csv(paste(tempdir,"OnCore_SNVIndel_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                   groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_CNV_OnCore <- 
    read.csv(paste(tempdir,"OnCore_CNV_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                   groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_Fusion_OnCore <- 
    read.csv(paste(tempdir,"OnCore_Fusion_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                   groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Specify empty dataframes
  #---------------------------------------------- 
  if (nrow(Output_SNVIndel_OnCore) > 0) {continue_SNVIndel <- as.logical("TRUE")
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  if (nrow(Output_CNV_OnCore) > 0) {continue_CNV <- as.logical("TRUE")
  } else {continue_CNV <- as.logical("FALSE")}
  
  if (nrow(Output_Fusion_OnCore) > 0) {continue_Fusion <- as.logical("TRUE")
  } else {continue_Fusion <- as.logical("FALSE")}
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.OnCore <- c("OnCore.No","Title","Biomarker.Description","Disease.Group","Disease.Sites","PI","PI.Email","Primary.CRC","Primary.CRC.Email")
  
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Type","Variant_Detail","VariantHGVSProtein")
  colnames.CNV.original <- c("PatientID","CNV_Gene","Variant_Type","var.anno","VariantHGVSProtein")
  colnames.Fusion.original <- c("PatientID","Gene","Variant_Type","Fusion_Detail","VariantHGVSProtein")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Type","Variant_Detail","VariantHGVSProtein")
  
  # Extract relevant columns from candidate matches
  #---------------------------------------------- 
  if (isTRUE(continue_SNVIndel)) {
    # Input missing columns 
    Output_SNVIndel_OnCore$Variant_Type = "SNV/Indel"
    Output_SNVIndel_OnCore$Variant_Detail = NA
    
    # Extract relevant columns
    Output_SNVIndel_OnCore <- Output_SNVIndel_OnCore[,c(colnames.SNVIndel.original,colnames.OnCore)]
    
    # Rename columns
    colnames(Output_SNVIndel_OnCore) <- c(colnames.output,colnames.OnCore)
  }
  
  if (isTRUE(continue_CNV)) {
    # Format back to original terminologies
    Output_CNV_OnCore$var.anno <- gsub("amplification","AMP",Output_CNV_OnCore$var.anno)
    Output_CNV_OnCore$var.anno <- gsub("deletion","DEL",Output_CNV_OnCore$var.anno)
    
    # Input missing columns 
    Output_CNV_OnCore$Variant_Type = "CNV"
    Output_CNV_OnCore$VariantHGVSProtein = NA
    
    # Extract relevant columns
    Output_CNV_OnCore <- Output_CNV_OnCore[,c(colnames.CNV.original,colnames.OnCore)]
    
    # Rename columns
    colnames(Output_CNV_OnCore) <- c(colnames.output,colnames.OnCore)
  }
  
  if (isTRUE(continue_Fusion)) {
    # Input missing columns 
    Output_Fusion_OnCore$Variant_Type = "Fusion"
    Output_Fusion_OnCore$VariantHGVSProtein = NA
    
    # Extract relevant columns
    Output_Fusion_OnCore <- Output_Fusion_OnCore[,c(colnames.Fusion.original,colnames.OnCore)]
    
    # Rename columns
    colnames(Output_Fusion_OnCore) <- c(colnames.output,colnames.OnCore)
  }
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_OnCore_FINAL <- data.frame(matrix(ncol = (length(colnames.output)+length(colnames.OnCore)), nrow = 0))
  colnames(Output_OnCore_FINAL) <- c(colnames.output,colnames.OnCore)
  
  if (isTRUE(continue_SNVIndel)) {
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_SNVIndel_OnCore)
  }
  
  if (isTRUE(continue_CNV)) {
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_CNV_OnCore)
  }
  
  if (isTRUE(continue_Fusion)) {
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_Fusion_OnCore)
  }
  
  # Remove duplicates
  Output_OnCore_FINAL <- unique(Output_OnCore_FINAL[order(Output_OnCore_FINAL$OnCore.No, decreasing = FALSE),])
  
  # Extract unique patient_id
  patient.oncore.matched <- sort(unique(Output_OnCore_FINAL$PatientID))
  
  remove(Output_SNVIndel_OnCore,Output_CNV_OnCore,Output_Fusion_OnCore,
         continue_SNVIndel,continue_CNV,continue_Fusion,colnames.output,
         colnames.OnCore,colnames.SNVIndel.original,colnames.CNV.original,colnames.Fusion.original)
  
} else {patient.oncore.matched <- c()}

if (isTRUE(NCI_match)) {
  email.comment <- paste("Patient Variant Report updated on ", Patient_Variant_Report_timestamp, ". ",
                         email.comment, sep="")
  
  # Import candidate matches
  #---------------------------------------------- 
  Output_SNVIndel_Variant <- 
    read.csv(paste(tempdir,"NCI_SNVIndel_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_CNV_Variant <- 
    read.csv(paste(tempdir,"NCIMatch_CNV_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_Fusion_Variant <- 
    read.csv(paste(tempdir,"NCI_Fusion_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_SNVIndel_NonHotspot <- 
    read.csv(paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  Output_CNV_NonHotspot <- 
    read.csv(paste(tempdir,"NCI_CNV_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Specify empty dataframes
  #---------------------------------------------- 
  if (nrow(Output_SNVIndel_Variant) > 0) {continue_SNVIndel_Variant <- as.logical("TRUE")
  } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
  
  if (nrow(Output_CNV_Variant) > 0) {continue_CNV_Variant <- as.logical("TRUE")
  } else {continue_CNV_Variant <- as.logical("FALSE")}
  
  if (nrow(Output_Fusion_Variant) > 0) {continue_Fusion_Variant <- as.logical("TRUE")
  } else {continue_Fusion_Variant <- as.logical("FALSE")}
  
  if (nrow(Output_SNVIndel_NonHotspot) > 0) {continue_SNVIndel_NonHotspot <- as.logical("TRUE")
  } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
  
  if (nrow(Output_CNV_NonHotspot) > 0) {continue_CNV_NonHotspot <- as.logical("TRUE")
  } else {continue_CNV_NonHotspot <- as.logical("FALSE")}
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.NCI.Variant <- c("Arm_Name","Gene_Name","Exon_Number","Protein","NCI_Type")
  colnames.NCI.NonHotspot <- c("Arm_Name","Gene_Name","Exon","oncominevariantclass_Function","NCI_Type")
  
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Detail","VariantHGVSProtein","Exon_Number","var.anno")
  colnames.CNV.original <- c("PatientID","CNV_Gene","Variant_Detail","VariantHGVSProtein","Exon_Number","var.anno")
  colnames.Fusion.original <- c("PatientID","Gene","Fusion_Detail","VariantHGVSProtein","Exon_Number","var.anno")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Detail","VariantHGVSProtein","Exon_Number","var.anno")
  colnames.NCI.output <- c("Arm_Name","Gene_Name","Exon","Biomarker_Detail","NCI_Type")
  
  # Extract relevant columns from candidate matches
  #---------------------------------------------- 
  if (isTRUE(continue_SNVIndel_Variant)) {
    # Input missing columns 
    Output_SNVIndel_Variant$Variant_Detail = NA
    Output_SNVIndel_Variant$NCI_Type = "Variant_Inclusion"
    
    # Extract relevant columns
    Output_SNVIndel_Variant <- Output_SNVIndel_Variant[,c(colnames.SNVIndel.original,colnames.NCI.Variant)]
    
    # Rename columns
    colnames(Output_SNVIndel_Variant) <- c(colnames.output,colnames.NCI.output)
  }
  
  if (isTRUE(continue_CNV_Variant)) {
    # Format back to original terminologies
    Output_CNV_Variant$var.anno <- gsub("amplification","AMP",Output_CNV_Variant$var.anno)
    Output_CNV_Variant$var.anno <- gsub("deletion","DEL",Output_CNV_Variant$var.anno)
    
    # Input missing columns 
    Output_CNV_Variant$Variant_Detail = NA
    Output_CNV_Variant$VariantHGVSProtein = NA
    Output_CNV_Variant$Exon_Number = NA
    Output_CNV_Variant$NCI_Type = "Variant_Inclusion"
    
    # Extract relevant columns
    Output_CNV_Variant <- Output_CNV_Variant[,c(colnames.CNV.original,colnames.NCI.Variant)]
    
    # Rename columns
    colnames(Output_CNV_Variant) <- c(colnames.output,colnames.NCI.output)
  }
  
  if (isTRUE(continue_Fusion_Variant)) {
    # Input missing columns 
    Output_Fusion_Variant$VariantHGVSProtein = NA
    Output_Fusion_Variant$Exon_Number = NA
    Output_Fusion_Variant$NCI_Type = "Variant_Inclusion"
    
    # Extract relevant columns
    Output_Fusion_Variant <- Output_Fusion_Variant[,c(colnames.Fusion.original,colnames.NCI.Variant)]
    
    # Rename columns
    colnames(Output_Fusion_Variant) <- c(colnames.output,colnames.NCI.output)
  }
  
  if (isTRUE(continue_SNVIndel_NonHotspot)) {
    # Input missing columns 
    Output_SNVIndel_NonHotspot$Variant_Detail = NA
    Output_SNVIndel_NonHotspot$oncominevariantclass_Function <- 
      tolower(gsub("NA","",paste(Output_SNVIndel_NonHotspot$oncominevariantclass,Output_SNVIndel_NonHotspot$Function)))
    Output_SNVIndel_NonHotspot$NCI_Type = "NonHotspot_Inclusion"
    
    # Extract relevant columns
    Output_SNVIndel_NonHotspot <- Output_SNVIndel_NonHotspot[,c(colnames.SNVIndel.original,colnames.NCI.NonHotspot)]
    
    # Rename columns
    colnames(Output_SNVIndel_NonHotspot) <- c(colnames.output,colnames.NCI.output)
  }
  
  if (isTRUE(continue_CNV_NonHotspot)) {
    # Format back to original terminologies
    Output_CNV_NonHotspot$var.anno <- gsub("amplification","AMP",Output_CNV_NonHotspot$var.anno)
    Output_CNV_NonHotspot$var.anno <- gsub("deletion","DEL",Output_CNV_NonHotspot$var.anno)
    
    # Input missing columns 
    Output_CNV_NonHotspot$VariantHGVSProtein = NA
    Output_CNV_NonHotspot$Variant_Detail = NA
    Output_CNV_NonHotspot$Exon_Number = NA
    Output_CNV_NonHotspot$oncominevariantclass_Function <- 
      tolower(gsub("NA","",paste(Output_CNV_NonHotspot$oncominevariantclass,Output_CNV_NonHotspot$Function)))
    Output_CNV_NonHotspot$NCI_Type = "NonHotspot_Inclusion"
    
    # Extract relevant columns
    Output_CNV_NonHotspot <- Output_CNV_NonHotspot[,c(colnames.CNV.original,colnames.NCI.NonHotspot)]
    
    # Rename columns
    colnames(Output_CNV_NonHotspot) <- c(colnames.output,colnames.NCI.output)
  }
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_NCI_FINAL <- data.frame(matrix(ncol = (length(colnames.output)+length(colnames.NCI.output)), nrow = 0))
  colnames(Output_NCI_FINAL) <- c(colnames.output,colnames.NCI.output)
  
  if (isTRUE(continue_SNVIndel_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_Variant)
  }
  
  if (isTRUE(continue_CNV_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_Variant)
  }
  
  if (isTRUE(continue_Fusion_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_Fusion_Variant)
  }
  
  if (isTRUE(continue_SNVIndel_NonHotspot)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_NonHotspot)
  }
  
  if (isTRUE(continue_CNV_NonHotspot)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_NonHotspot)
  }
  
  # Remove duplicates
  Output_NCI_FINAL <- unique(Output_NCI_FINAL[order(Output_NCI_FINAL$Arm_Name, decreasing = FALSE),])
  
  # Extract unique patient_id
  patient.nci.matched <- sort(unique(Output_NCI_FINAL$PatientID))
  
  remove(Output_SNVIndel_Variant,Output_CNV_Variant,Output_Fusion_Variant,Output_SNVIndel_NonHotspot,Output_CNV_NonHotspot,
         continue_SNVIndel_Variant,continue_CNV_Variant,continue_Fusion_Variant,continue_SNVIndel_NonHotspot,continue_CNV_NonHotspot,
         colnames.NCI.Variant,colnames.NCI.NonHotspot,colnames.SNVIndel.original,colnames.CNV.original,colnames.Fusion.original,
         colnames.output,colnames.NCI.output)
  
} else {patient.nci.matched <- c()}

# Compile complete list of patient ids
patient.list.matched <- unique(c(patient.oncore.matched,patient.nci.matched))

# Iterate through each patient_id of patient.list.matched
#---------------------------------------------- 
if (isTRUE(length(patient.list.matched) > 0)) {
  
  for (patient_num  in 1:length(patient.list.matched)) {
    patient_id <- patient.list.matched[patient_num]
    
    # Extract candidate trials based on patient_id
    #----------------------------------------------
    if (length(patient.oncore.matched) > 0) {Output_OnCore_FINAL <- Output_OnCore_FINAL[which(Output_OnCore_FINAL$PatientID == patient_id), ]}
    if (length(patient.nci.matched) > 0) {Output_NCI_FINAL <- Output_NCI_FINAL[which(Output_NCI_FINAL$PatientID == patient_id), ]}
    
    ## Write output to file (detailed match results per variant)
    #----------------------------------------------
    sink(file = paste(outdir, patient_id, ".PositiveMatch.txt", sep=""),
         append = FALSE, split = FALSE)
    options(max.print=999999)
    
    # Output patient bio
    cat(paste(patient_id, " may qualify for the following clinical trial(s) due to mutations identified in STAMP assay.",
              sep=""), "\n","\n")
    
    # Output trial INFO from Matched_Internal
    #----------------------------------------------
    if (isTRUE(Internal_match & nrow(Output_OnCore_FINAL) > 0 )) {
      
      # Extract list of matched OnCore.No
      OnCoreNo.list <- sort(unique(Output_OnCore_FINAL$OnCore.No))
      
      # Iterate through each matched trial_id per patient
      #----------------------------------------------
      for (num_internal in 1:length(OnCoreNo.list)) {
        trial_id <- OnCoreNo.list[num_internal]
        
        # Extract trial INFO based on trial_id match in Matched_Internal
        DF_candidate_trial <- 
          unique(Output_OnCore_FINAL[which(Output_OnCore_FINAL$OnCore.No == trial_id), ])
        
        # Remove trailing whitespace
        DF_candidate_trial$Biomarker.Description <- gsub("\n$", "", DF_candidate_trial$Biomarker.Description)
        
        # Output notification
        #----------------------------------------------
        cat(paste("No.", num_internal, ": Clinical Trial #", trial_id, ": ", unique(DF_candidate_trial$Title), sep=""),"\n",
            "----------------------------------------------------------------------", "\n")
        cat(paste("Principal Investigator: ", unique(DF_candidate_trial$PI), " (",
                  unique(DF_candidate_trial$PI.Email), ")", sep=""),"\n")
        cat(paste("Primary Clinical Research Coordinator: ", unique(DF_candidate_trial$Primary.CRC),
                  " (", unique(DF_candidate_trial$Primary.CRC.Email), ")", sep=""),"\n","\n")
        cat("Biomarker criteria:", "\n", "\t", gsub("\n","",unique(DF_candidate_trial$Biomarker.Description)),"\n","\n")
        cat("Relevant mutations identified in patient:", "\n")
        
        for (entry_num in 1:nrow(DF_candidate_trial)) {
          if (DF_candidate_trial$Variant_Type[entry_num] == "SNV/Indel") {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], ": ", DF_candidate_trial$VariantHGVSProtein[entry_num], sep="")
          } else if (DF_candidate_trial$Variant_Type[entry_num] == "Fusion") {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], " (", DF_candidate_trial$Variant_Detail[entry_num], " Fusion)", sep="")
          } else if (DF_candidate_trial$Variant_Type[entry_num] == "CNV") {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], " (", DF_candidate_trial$Variant_Detail[entry_num], " CNV)", sep="")
          } else {
            DF_candidate_trial$Output[entry_num] <- 
              paste("Potential ERROR in pipeline - manual assessment is needed")
          }
        }
        
        Output.patient <-  unique(DF_candidate_trial$Output)
        for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n")}
        cat("\n")
        
        cat("Trial criteria to be manually assessed: ", "\n")
        cat("\t",paste("Disease group criteria: ", unique(DF_candidate_trial$Disease.Group), sep=""),"\n")
        cat("\t",paste("Disease site criteria: ", unique(DF_candidate_trial$Disease.Sites), sep=""),"\n","\n","\n")
        
        remove(trial_id,DF_candidate_trial,Output.patient,entry_num)
      }
      remove(Output_OnCore_FINAL,OnCoreNo.list)
      
    } else {num_internal = 0}
    
    # Output trial INFO from Matched_NCI.Variants
    #----------------------------------------------
    if (isTRUE(NCI_match & nrow(Output_NCI_FINAL) > 0)) {
      
      # Extract list of matched Arm_No
      NCIArm.Variant.list <- sort(unique(Output_NCI_FINAL$Arm_Name))
      
      # Iterate through each matched trial_id per patient
      #----------------------------------------------
      for (num_variant in 1:length(NCIArm.Variant.list)) {
        trial_id <- NCIArm.Variant.list[num_variant]
        
        # Extract trial INFO based on trial_id match in Matched_NCI.Variants
        DF_candidate_trial <- 
          unique(Output_NCI_FINAL[which(Output_NCI_FINAL$Arm_Name == trial_id), ])
        
        # Output notification
        #----------------------------------------------
        cat(paste("No.", (num_internal + num_variant), ": NCI-MATCH Trial Treatment ", trial_id, sep=""),"\n",
            "----------------------------------------------------------------------", "\n")
        
        for (entry_num in 1:nrow(DF_candidate_trial)) {
          if (DF_candidate_trial$NCI_Type[entry_num] == "Variant_Inclusion") {
            
            if (DF_candidate_trial$var.anno[entry_num] == "Fusion") {
              DF_candidate_trial$Output[entry_num] <- 
                paste(DF_candidate_trial$Gene_Name[entry_num], " Fusion (refer to trial criteria for details)", sep="")  
              
            } else if (DF_candidate_trial$var.anno[entry_num] %in% c("AMP","DEL")) {
              DF_candidate_trial$Output[entry_num] <- 
                paste(DF_candidate_trial$Gene_Name[entry_num], " ", DF_candidate_trial$Biomarker_Detail[entry_num], sep="")
              
            } else {
              DF_candidate_trial$Output[entry_num] <- 
                paste(DF_candidate_trial$Gene_Name[entry_num], ": ", DF_candidate_trial$Biomarker_Detail[entry_num], sep="")  
            }
            
          } else if (DF_candidate_trial$NCI_Type[entry_num] == "NonHotspot_Inclusion") {
            DF_candidate_trial$Output[entry_num] <- 
              gsub("Exon NA ","", paste(DF_candidate_trial$Gene_Name[entry_num], " Exon ", DF_candidate_trial$Exon[entry_num], " (",
                    gsub("[[:blank:]]$","",gsub("^[[:blank:]]+","",DF_candidate_trial$Biomarker_Detail[entry_num])), ")", sep=""))
            
          } else {
            DF_candidate_trial$Output[entry_num] <- 
              paste("Potential ERROR in pipeline - manual assessment is needed")
          }
        }
        
        Output.patient <-  unique(DF_candidate_trial$Output[which(DF_candidate_trial$NCI_Type == "Variant_Inclusion")])
        if (length(Output.patient) > 0) {
          cat("Variant inclusion criteria: ", "\n")
          for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n")}  
        }
        
        Output.patient <-  unique(DF_candidate_trial$Output[which(DF_candidate_trial$NCI_Type == "NonHotspot_Inclusion")])
        if (length(Output.patient) > 0) {
          cat("Nonhotspot rule inclusion criteria: ", "\n")
          for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n")}
        }
        
        cat("\n")
        cat("Relevant mutations identified in patient:", "\n")
        
        for (entry_num in 1:nrow(DF_candidate_trial)) {
          if (DF_candidate_trial$var.anno[entry_num] == "MUTATION") {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], " Exon ", DF_candidate_trial$Exon_Number[entry_num], ": ", 
                    DF_candidate_trial$VariantHGVSProtein[entry_num], sep="")
            
            
          } else if (DF_candidate_trial$var.anno[entry_num] == "Fusion") {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], " (", DF_candidate_trial$Variant_Detail[entry_num], " Fusion)", sep="")
          } else if (DF_candidate_trial$var.anno[entry_num] %in% c("AMP","DEL")) {
            DF_candidate_trial$Output[entry_num] <- 
              paste(DF_candidate_trial$Variant_Gene[entry_num], " (", DF_candidate_trial$var.anno[entry_num], " CNV)", sep="")
          } else {
            DF_candidate_trial$Output[entry_num] <- 
              paste("Potential ERROR in pipeline - manual assessment is needed")
          }
        }
        
        Output.patient <-  unique(DF_candidate_trial$Output)
        for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n")}
        cat("\n")
        
        # Output additional ARM criteria
        #---------------------------------------------- 
        cat("Trial criteria to be manually assessed: ", "\n")
        
        assess_comment = ""
        
        # Exclusion_NonHotspot_Rules matched by trial_id
        #----------------------------------------------
        Exclusion_NonHotspot_Output <-
          Exclusion_NonHotspot_Rules[Exclusion_NonHotspot_Rules$Arm_Name == trial_id,]
        
        if (nrow(Exclusion_NonHotspot_Output) > 0) {
          
          ## Output applicable Exclusion NonHotspot rules
          #----------------------------------------------
          cat("Exclusion NonHotspot Rules:", "\n")
          
          Exclusion_NonHotspot_Output$Output <- paste(Exclusion_NonHotspot_Output$Gene_Name, " Exon ", Exclusion_NonHotspot_Output$Exon, " (",
                                                      Exclusion_NonHotspot_Output$Function, tolower(Exclusion_NonHotspot_Output$oncominevariantclass), ")", sep="")
          Exclusion_NonHotspot_Output$Output <- gsub("\\(NA","\\(", gsub("NA\\)$","\\)",gsub(" Exon NA ", " ",Exclusion_NonHotspot_Output$Output)))
          
          Output.trial <-  unique(Exclusion_NonHotspot_Output$Output)
          for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n")}
          cat("\n")
          
        } else {
          assess_comment <- paste(assess_comment, "There are no exclusion nonhotspot rules indicated for ",trial_id,". ",sep="")
        }
        
        # IHC_Results matched by trial_id
        #----------------------------------------------
        IHC_Output <- IHC_Results[IHC_Results$Arm_Name == trial_id,]
        
        if (nrow(IHC_Output) > 0) {
          
          ## Output applicable IHC_Results criteria
          #----------------------------------------------
          cat("IHC results: ", "\n")
          Output.trial <- unique(paste(IHC_Output$Gene,
                                       " (Status: ", IHC_Output$Status_POSITIVE_NEGATIVE_INDETERMINATE,
                                       "; Variant: ", IHC_Output$Variant_PRESENT_NEGATIVE_EMPTY, ")", sep=""))
          for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
          cat("\n")
          
        } else {
          assess_comment <- paste(assess_comment, "There are no immunohistochemical (IHC) assay criteria indicated for ",trial_id,". ",sep="")
        }
        
        # # Comments matched by trial_id
        # #----------------------------------------------
        # Comments_Output <- Comments[Comments$Arm_Name == trial_id,]
        # 
        # if (nrow(Comments_Output) > 0) {
        #   
        #   ## Output applicable Comments criteria
        #   #----------------------------------------------
        #   cat("Comments:", "\n")
        #   print(Comments_Output, row.names = FALSE)
        #   cat("\n")
        #   
        # } else {
        #   assess_comment <- paste(assess_comment, "There are no comments indicated for ",trial_id,". ",sep="")
        # }
        
        # Disease_Exclusion_Codes matched by trial_id
        #----------------------------------------------
        Disease_Exclusion_Output <- Disease_Exclusion_Codes[Disease_Exclusion_Codes$Arm_Name == trial_id,]
        
        if (nrow(Disease_Exclusion_Output) > 0) {
          
          ## Output applicable Disease_Exclusion_Codes criteria
          #----------------------------------------------
          cat("Histological disease exclusion codes:","\n")
          Output.trial <- sort(unique(Disease_Exclusion_Output$SHORT.NAME))
          for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
        
        } else {
          assess_comment <- paste(assess_comment, "There are no histological disease exclusion codes indicated for ",trial_id,".",sep="")
        }
        
        if (nchar(assess_comment) > 0) {
          cat("\n",assess_comment,"\n") 
        }
        cat("\n","\n")
        
        remove(trial_id,DF_candidate_trial,entry_num,Output.patient,Exclusion_NonHotspot_Output,
               IHC_Output,Disease_Exclusion_Output)
        if (isTRUE(exists("Output.trial"))){remove(Output.trial)}
      }
      remove(Output_NCI_FINAL,NCIArm.Variant.list,num_variant)
    }
    
    cat("NOTES:","\n")
    if (isTRUE(Internal_match)) {
      cat("[Stanford Internal Trials] Trial criteria indicated by age group, variant pathogenicity, disease group, disease site and comments need to be manually assessed.","\n")
    }
    if (isTRUE(NCI_match)) {
      cat("[NCI-MATCH Trials] Trial criteria indicated by age group, variant pathogenicity, exclusion nonhotspot rules, IHC results, comments, and disease exclusions need to be manually assessed.","\n")
      cat("[NCI-MATCH Trials] If candidate trials have been matched based on inclusion nonhotspot rules, the pathogenicity status and variant type of the listed variants need to be manually assessed.","\n")
    }
    cat("\n")
    cat(email.comment,"\n")
    
    sink()
    
    remove(patient_id,num_internal)
  }
  remove(patient_num)
}

remove(email.comment,patient.oncore.matched,patient.nci.matched,patient.list.matched)

cat(paste("Timestamp of algorithm matching output FINISH: ", Sys.time(), sep=""),"\n","\n")
