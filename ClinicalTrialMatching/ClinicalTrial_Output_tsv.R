# Generate corresponding tsv files
# patient_id.OnCore.tsv and patient_id.NCI.tsv

#---------------------------------------------- 
## GENERATE tsv output for positive matches per patient = ONCORE
#---------------------------------------------- 
if (isTRUE(Internal_match)) {
  
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
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Type","Variant_Detail",
                                  "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                  "OnCore.No","NCT..")
  
  colnames.CNV.original <- c("PatientID","CNV_Gene","Variant_Type","var.anno",
                             "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                             "OnCore.No","NCT..")
  
  colnames.Fusion.original <- c("PatientID","Gene","Variant_Type","Fusion_Detail",
                                "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                "OnCore.No","NCT..")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Type","Variant_Detail",
                       "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                       "OnCore.No","NCT.No")
  
  # Iterate through each patient_id of patient.oncore.matched
  #---------------------------------------------- 
  for (patient_num  in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    if (isTRUE(continue_SNVIndel)) {
      # Extract candidate trials for patient
      Output_SNVIndel_OnCore <- 
        Output_SNVIndel_OnCore[which(Output_SNVIndel_OnCore$PatientID == patient_id),]
      
      # Input missing columns 
      Output_SNVIndel_OnCore$Variant_Type = "SNV/Indel"
      Output_SNVIndel_OnCore$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_OnCore <- Output_SNVIndel_OnCore[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_OnCore) <- colnames.output
    }
    
    if (isTRUE(continue_CNV)) {
      # Extract candidate trials for patient
      Output_CNV_OnCore <- 
        Output_CNV_OnCore[which(Output_CNV_OnCore$PatientID == patient_id),]
      
      # Format back to original terminologies
      Output_CNV_OnCore$var.anno <- gsub("amplification","AMP",Output_CNV_OnCore$var.anno)
      Output_CNV_OnCore$var.anno <- gsub("deletion","DEL",Output_CNV_OnCore$var.anno)
      
      # Input missing columns 
      Output_CNV_OnCore$Variant_Type = "CNV"
      Output_CNV_OnCore$VariantHGVSCoding = NA
      Output_CNV_OnCore$VariantHGVSProtein = NA
      Output_CNV_OnCore$VariantHGVSGenomic = NA
      
      # Extract relevant columns
      Output_CNV_OnCore <- Output_CNV_OnCore[,colnames.CNV.original]
      
      # Rename columns
      colnames(Output_CNV_OnCore) <- colnames.output
    }
    
    if (isTRUE(continue_Fusion)) {
      # Extract candidate trials for patient
      Output_Fusion_OnCore <- 
        Output_Fusion_OnCore[which(Output_Fusion_OnCore$PatientID == patient_id),]
      
      # Input missing columns 
      Output_Fusion_OnCore$Variant_Type = "Fusion"
      Output_Fusion_OnCore$VariantHGVSCoding = NA
      Output_Fusion_OnCore$VariantHGVSProtein = NA
      Output_Fusion_OnCore$VariantHGVSGenomic = NA
      
      # Extract relevant columns
      Output_Fusion_OnCore <- Output_Fusion_OnCore[,colnames.Fusion.original]
      
      # Rename columns
      colnames(Output_Fusion_OnCore) <- colnames.output
    }
    
    # Merge into single dataframe
    #---------------------------------------------- 
    Output_OnCore_FINAL <- data.frame(matrix(ncol = length(colnames.output), nrow = 0))
    colnames(Output_OnCore_FINAL) <- colnames.output
    
    if (isTRUE(continue_SNVIndel)) {
      Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_SNVIndel_OnCore)
    }
    
    if (isTRUE(continue_CNV)) {
      Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_CNV_OnCore)
    }
    
    if (isTRUE(continue_Fusion)) {
      Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_Fusion_OnCore)
    }
    
    if (nrow(Output_OnCore_FINAL) > 0){
      # Write to match results to local computer
      write.table(Output_OnCore_FINAL, file = paste(outdir, patient_id, ".OnCore.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE, na = ".")
    }
    remove(Output_OnCore_FINAL,patient_id)
  }
  remove(patient_num,Output_SNVIndel_OnCore,Output_CNV_OnCore,Output_Fusion_OnCore,
         continue_SNVIndel,continue_CNV,continue_Fusion,colnames.output,
         colnames.SNVIndel.original,colnames.CNV.original,colnames.Fusion.original)
}


#---------------------------------------------- 
## GENERATE tsv output for positive matches per patient = NCI-MATCH
#---------------------------------------------- 
if (isTRUE(NCI_match)) {
  
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
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Type","Variant_Detail",
                                  "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                  "Arm_Name","NCT.No")
  
  colnames.CNV.original <- c("PatientID","CNV_Gene","var.type","var.anno",
                             "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                             "Arm_Name","NCT.No")
  
  colnames.Fusion.original <- c("PatientID","Gene","Variant_Type","Fusion_Detail",
                                "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                "Arm_Name","NCT.No")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Type","Variant_Detail",
                       "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                       "Arm_Name","NCT.No")
  
  # Iterate through each patient_id of patient.oncore.matched
  #---------------------------------------------- 
  for (patient_num  in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    if (isTRUE(continue_SNVIndel_Variant)) {
      # Extract candidate trials for patient
      Output_SNVIndel_Variant <- 
        Output_SNVIndel_Variant[which(Output_SNVIndel_Variant$PatientID == patient_id),]
      
      # Input missing columns 
      Output_SNVIndel_Variant$Variant_Type = "SNV/Indel"
      Output_SNVIndel_Variant$NCT.No = "NCT02465060"
      Output_SNVIndel_Variant$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_Variant <- Output_SNVIndel_Variant[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_Variant) <- colnames.output
    }
    
    if (isTRUE(continue_CNV_Variant)) {
      # Extract candidate trials for patient
      Output_CNV_Variant <- 
        Output_CNV_Variant[which(Output_CNV_Variant$PatientID == patient_id),]
      
      # Format back to original terminologies
      Output_CNV_Variant$var.anno <- gsub("amplification","AMP",Output_CNV_Variant$var.anno)
      Output_CNV_Variant$var.anno <- gsub("deletion","DEL",Output_CNV_Variant$var.anno)
      
      # Input missing columns 
      Output_CNV_Variant$Variant_Type = "CNV"
      Output_CNV_Variant$NCT.No = "NCT02465060"
      Output_CNV_Variant$VariantHGVSCoding = NA
      Output_CNV_Variant$VariantHGVSProtein = NA
      Output_CNV_Variant$VariantHGVSGenomic = NA
      
      # Extract relevant columns
      Output_CNV_Variant <- Output_CNV_Variant[,colnames.CNV.original]
      
      # Rename columns
      colnames(Output_CNV_Variant) <- colnames.output
    }
    
    if (isTRUE(continue_Fusion_Variant)) {
      # Extract candidate trials for patient
      Output_Fusion_Variant <- 
        Output_Fusion_Variant[which(Output_Fusion_Variant$PatientID == patient_id),]
      
      # Input missing columns 
      Output_Fusion_Variant$NCT.No = "NCT02465060"
      Output_Fusion_Variant$VariantHGVSCoding = NA
      Output_Fusion_Variant$VariantHGVSProtein = NA
      Output_Fusion_Variant$VariantHGVSGenomic = NA
      
      # Extract relevant columns
      Output_Fusion_Variant <- Output_Fusion_Variant[,colnames.Fusion.original]
      
      # Rename columns
      colnames(Output_Fusion_Variant) <- colnames.output
    }
    
    if (isTRUE(continue_SNVIndel_NonHotspot)) {
      # Extract candidate trials for patient
      Output_SNVIndel_NonHotspot <- 
        Output_SNVIndel_NonHotspot[which(Output_SNVIndel_NonHotspot$PatientID == patient_id),]
      
      # Input missing columns 
      Output_SNVIndel_NonHotspot$Variant_Type = "SNV/Indel"
      Output_SNVIndel_NonHotspot$NCT.No = "NCT02465060"
      Output_SNVIndel_NonHotspot$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_NonHotspot <- Output_SNVIndel_NonHotspot[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_NonHotspot) <- colnames.output
    }
    
    if (isTRUE(continue_CNV_NonHotspot)) {
      # Extract candidate trials for patient
      Output_CNV_NonHotspot <- 
        Output_CNV_NonHotspot[which(Output_CNV_NonHotspot$PatientID == patient_id),]
      
      # Format back to original terminologies
      Output_CNV_NonHotspot$var.anno <- gsub("amplification","AMP",Output_CNV_NonHotspot$var.anno)
      Output_CNV_NonHotspot$var.anno <- gsub("deletion","DEL",Output_CNV_NonHotspot$var.anno)
      
      # Input missing columns 
      Output_CNV_NonHotspot$NCT.No = "NCT02465060"
      Output_CNV_NonHotspot$VariantHGVSCoding = NA
      Output_CNV_NonHotspot$VariantHGVSProtein = NA
      Output_CNV_NonHotspot$VariantHGVSGenomic = NA
      
      # Extract relevant columns
      Output_CNV_NonHotspot <- Output_CNV_NonHotspot[,colnames.CNV.original]
      
      # Rename columns
      colnames(Output_CNV_NonHotspot) <- colnames.output
    }
    
    # Merge into single dataframe
    #---------------------------------------------- 
    Output_NCI_FINAL <- data.frame(matrix(ncol = length(colnames.output), nrow = 0))
    colnames(Output_NCI_FINAL) <- colnames.output
    
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
    
    if (nrow(Output_NCI_FINAL) > 0){
      # Remove duplicates
      Output_NCI_FINAL <- unique(Output_NCI_FINAL[,])
      
      # Write to match results to local computer
      write.table(Output_NCI_FINAL, file = paste(outdir, patient_id, ".NCI.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE, na = ".")  
    }
    remove(patient_id,Output_NCI_FINAL)
  }
  remove(Output_SNVIndel_Variant,Output_CNV_Variant,Output_Fusion_Variant,Output_SNVIndel_NonHotspot,Output_CNV_NonHotspot,
         continue_SNVIndel_Variant,continue_CNV_Variant,continue_Fusion_Variant,continue_SNVIndel_NonHotspot,continue_CNV_NonHotspot,
         colnames.SNVIndel.original,colnames.CNV.original,colnames.Fusion.original,colnames.output,patient_num)
}
