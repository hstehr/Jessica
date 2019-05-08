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
  
  # Specify empty dataframes
  #---------------------------------------------- 
  if (nrow(Output_SNVIndel_OnCore) > 0) {continue_SNVIndel <- as.logical("TRUE")
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Type","Variant_Detail",
                                  "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                  "OnCore.No","NCT..")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Type","Variant_Detail",
                       "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                       "OnCore.No","NCT.No")
  
  # Iterate through each patient_id of patient.oncore.matched
  #---------------------------------------------- 
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    # Extract patient-specific entries
    Output_SNVIndel_OnCore_Patient <- Output_SNVIndel_OnCore[which(Output_SNVIndel_OnCore$PatientID == patient_id),]
    
    if (isTRUE(continue_SNVIndel & nrow(Output_SNVIndel_OnCore_Patient) > 0)) {
      # Input missing columns 
      Output_SNVIndel_OnCore_Patient$Variant_Type = "SNV/Indel"
      Output_SNVIndel_OnCore_Patient$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_OnCore_Patient <- Output_SNVIndel_OnCore_Patient[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_OnCore_Patient) <- colnames.output
    }
    
    # Merge into single dataframe
    #---------------------------------------------- 
    Output_OnCore_FINAL <- data.frame(matrix(ncol = length(colnames.output), nrow = 0))
    colnames(Output_OnCore_FINAL) <- colnames.output
    
    if (isTRUE(continue_SNVIndel)) {
      Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_SNVIndel_OnCore_Patient)
    }
    
    if (nrow(Output_OnCore_FINAL) > 0){
      # Write to match results to local computer
      write.table(Output_OnCore_FINAL, file = paste(outdir, patient_id, ".OnCore.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE, na = ".")
    }
    remove(Output_OnCore_FINAL,patient_id,Output_SNVIndel_OnCore_Patient)
  }
  remove(patient_num,Output_SNVIndel_OnCore,continue_SNVIndel,colnames.output,colnames.SNVIndel.original)
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
  
  Output_SNVIndel_NonHotspot <- 
    read.csv(paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                   dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Specify empty dataframes
  #---------------------------------------------- 
  if (nrow(Output_SNVIndel_Variant) > 0) {continue_SNVIndel_Variant <- as.logical("TRUE")
  } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
  
  if (nrow(Output_SNVIndel_NonHotspot) > 0) {continue_SNVIndel_NonHotspot <- as.logical("TRUE")
  } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.SNVIndel.original <- c("PatientID","VariantGene","Variant_Type","Variant_Detail",
                                  "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                                  "Arm_Name","NCT.No")
  
  colnames.output <- c("PatientID","Variant_Gene","Variant_Type","Variant_Detail",
                       "VariantHGVSCoding","VariantHGVSProtein","VariantHGVSGenomic",
                       "Arm_Name","NCT.No")
  
  # Iterate through each patient_id of patient.oncore.matched
  #---------------------------------------------- 
  for (patient_num  in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    # Extract patient-specific entries
    Output_SNVIndel_Variant_Patient <- Output_SNVIndel_Variant[which(Output_SNVIndel_Variant$PatientID == patient_id),]
    
    Output_SNVIndel_NonHotspot_Patient <- Output_SNVIndel_NonHotspot[which(Output_SNVIndel_NonHotspot$PatientID == patient_id),]
      
    if (isTRUE(continue_SNVIndel_Variant & nrow(Output_SNVIndel_Variant_Patient) > 0)) {
      # Input missing columns 
      Output_SNVIndel_Variant_Patient$Variant_Type = "SNV/Indel"
      Output_SNVIndel_Variant_Patient$NCT.No = "NCT02465060"
      Output_SNVIndel_Variant_Patient$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_Variant_Patient <- Output_SNVIndel_Variant_Patient[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_Variant_Patient) <- colnames.output
    }
    
    if (isTRUE(continue_SNVIndel_NonHotspot & nrow(Output_SNVIndel_NonHotspot_Patient) > 0)) {
      # Input missing columns 
      Output_SNVIndel_NonHotspot_Patient$Variant_Type = "SNV/Indel"
      Output_SNVIndel_NonHotspot_Patient$NCT.No = "NCT02465060"
      Output_SNVIndel_NonHotspot_Patient$Variant_Detail = NA
      
      # Extract relevant columns
      Output_SNVIndel_NonHotspot_Patient <- Output_SNVIndel_NonHotspot_Patient[,colnames.SNVIndel.original]
      
      # Rename columns
      colnames(Output_SNVIndel_NonHotspot_Patient) <- colnames.output
    }
    
    # Merge into single dataframe
    #---------------------------------------------- 
    Output_NCI_FINAL <- data.frame(matrix(ncol = length(colnames.output), nrow = 0))
    colnames(Output_NCI_FINAL) <- colnames.output
    
    if (isTRUE(continue_SNVIndel_Variant & nrow(Output_SNVIndel_NonHotspot_Patient) > 0)) {
      Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_Variant_Patient)
    }
    
    if (isTRUE(continue_SNVIndel_NonHotspot & nrow(Output_SNVIndel_NonHotspot_Patient) > 0)) {
      Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_NonHotspot_Patient)
    }
    
    if (nrow(Output_NCI_FINAL) > 0) {
      # Remove duplicates
      Output_NCI_FINAL <- unique(Output_NCI_FINAL[,])
      
      # Write to match results to local computer
      write.table(Output_NCI_FINAL, file = paste(outdir, patient_id, ".NCI.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE, na = ".")  
    }
    remove(patient_id,Output_NCI_FINAL,Output_SNVIndel_Variant_Patient,Output_SNVIndel_NonHotspot_Patient)
  }
  remove(Output_SNVIndel_Variant,Output_SNVIndel_NonHotspot,continue_SNVIndel_Variant,continue_SNVIndel_NonHotspot,
         colnames.SNVIndel.original,colnames.output,patient_num)
}
