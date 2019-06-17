## GENERATE detailed match results per patient

# Iterate through each patient_id of patient.list
for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient & specify empty dataframes
  #---------------------------------------------- 
  # SNV Indels
  file_name = paste(tempdir, patient_id, "_SNVIndel.tsv", sep="")
  
  if (file.exists(file_name)) {
    DF_patient_SNVIndel <- read.csv(file = file_name, header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient_SNVIndel) > 0) {continue_SNVIndel <- as.logical("TRUE")
    } else {continue_SNVIndel <- as.logical("FALSE")}
    
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  # CNVs
  file_name = paste(tempdir, patient_id, "_CNV.tsv", sep="")
  
  if (file.exists(file_name)) {
    DF_patient_CNV <- read.csv(file = file_name, header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient_CNV) > 0) {continue_CNV <- as.logical("TRUE")
    } else {continue_CNV <- as.logical("FALSE")}
    
  } else {continue_CNV <- as.logical("FALSE")}
  
  # Fusions
  file_name = paste(tempdir, patient_id, "_Fusion.tsv", sep="")
  
  if (file.exists(file_name)) {
    DF_patient_Fusion <- read.csv(file = file_name, header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient_Fusion) > 0) {continue_Fusion <- as.logical("TRUE")
    } else {continue_Fusion <- as.logical("FALSE")}
    
  } else {continue_Fusion <- as.logical("FALSE")}
  
  remove(file_name)
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.SNVIndel <- c("Output")
  colnames.CNV <- c("Output")
  colnames.Fusion <- c("Output")
  
  email.comment <- paste("This email was generated on ", Sys.time(), ".",sep="")
  trial.comment <- c()
  
  if (isTRUE(Internal_match)) {
    colnames.SNVIndel <- append(colnames.SNVIndel, c("OnCore_SNVIndel_Status"))
    colnames.CNV <- append(colnames.CNV, c("OnCore_CNV_Status"))
    colnames.Fusion <- append(colnames.Fusion, c("OnCore_Fusion_Status"))
    
    email.comment <- paste("OnCore Biomarker Report updated on ", OnCore_Biomarker_Report_timestamp, ". ",
                           email.comment, sep="")
    trial.comment <- append(trial.comment, "Stanford OnCore")
  }
  
  if (isTRUE(NCI_match)) {
    colnames.SNVIndel <- append(colnames.SNVIndel, c("NCI_SNVIndel_Variant_Status","NCI_SNVIndel_NonHotspot_Status"))
    colnames.CNV <- append(colnames.CNV, c("NCI_CNV_Variant_Status","NCI_CNV_NonHotspot_Status"))
    colnames.Fusion <- append(colnames.Fusion, c("NCI_Fusion_Variant_Status"))
    
    email.comment <- paste("Patient Variant Report updated on ", Patient_Variant_Report_timestamp, ". ",
                           email.comment, sep="")
    trial.comment <- append(trial.comment, "NCI-MATCH")
  }
  
  if (isTRUE(length(trial.comment) > 1)) {trial.comment <- paste(trial.comment[[1]],trial.comment[[2]],sep=" and ")
  } else {trial.comment <- trial.comment[[1]]}
  
  # Format dataframe for output
  #---------------------------------------------- 
  if (isTRUE(continue_SNVIndel)) {
    # Specify gene info for output
    DF_patient_SNVIndel$Output <- paste(DF_patient_SNVIndel$VariantGene, " ", DF_patient_SNVIndel$VariantHGVSProtein, sep="")
    
    # Extract relevant columns for output
    DF_patient_SNVIndel <- DF_patient_SNVIndel[,colnames.SNVIndel]
  }
  
  if (isTRUE(continue_CNV)) {
    # Format back to original terminologies
    DF_patient_CNV$var.anno <- gsub("amplification","AMP",DF_patient_CNV$var.anno)
    DF_patient_CNV$var.anno <- gsub("deletion","DEL",DF_patient_CNV$var.anno)
    
    # Specify gene info for output
    DF_patient_CNV$Output <- paste(DF_patient_CNV$CNV_Gene, " (", DF_patient_CNV$var.anno, ")", sep="")
    
    # Extract relevant columns for output
    DF_patient_CNV <- DF_patient_CNV[,colnames.CNV]
  }
  
  if (isTRUE(continue_Fusion)) {
    # Specify gene info for output
    DF_patient_Fusion$Output <- paste(DF_patient_Fusion$Gene, " (", DF_patient_Fusion$Fusion_Detail, " Fusion)", sep="")
    
    # Extract relevant columns for output
    DF_patient_Fusion <- DF_patient_Fusion[,colnames.Fusion]
  }
  
  #----------------------------------------------
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir, patient_id, ".Output.txt", sep=""), 
       append = FALSE, split = FALSE)
  options(max.print=999999)
  
  # Output patient bio
  cat(paste(patient_id, " may qualify for the following ",trial.comment,
            " clinical trial(s) due to mutation(s) identified in the STAMP assay.",sep=""),"\n","\n")
  
  if (isTRUE(Internal_match | NCI_match)) {
    if (isTRUE(continue_SNVIndel)) {
      cat("[Stanford OnCore] Order of matching (SNV Indels): biomarker (gene > condition > detail).","\n")
      cat("[Stanford OnCore] Trial criteria indicated by age group, variant pathogenicity, disease group, disease site and comments need to be manually assessed.","\n","\n")
      
      cat("[NCI-MATCH Inclusion Variants] Order of matching (SNV Indels): Gene > Variant Type > Genomic Region.","\n")
      cat("[NCI-MATCH Nonhotspot Rules] Order of matching (SNV Indels): Gene > function/oncominevariantclass > Exon No.","\n")
      cat("[NCI-MATCH] Trial criteria indicated by age group, variant pathogenicity, nonhotspot rules, IHC results, comments, and disease exclusions need to be manually assessed.","\n","\n")
      
      print(DF_patient_SNVIndel[order(DF_patient_SNVIndel$Output), ], row.names = TRUE)
    } else {
      cat("The STAMP assay did not identify any SNV/Indels in the solid tumor biopsied from the patient.","\n")
    }
    cat("\n","\n")
    
    if (isTRUE(continue_CNV)) {
      cat("[Stanford OnCore] Order of matching (copy number variations): biomarker (gene > condition).","\n")
      cat("[Stanford OnCore] Trial criteria indicated by age group, variant pathogenicity, disease group, disease site and comments need to be manually assessed.","\n","\n")
      
      cat("[NCI-MATCH Inclusion Variants] Order of matching (copy number variations): Gene > Variant Type (ie. amplification, deletion).","\n")
      cat("[NCI-MATCH Nonhotspot Rules] Order of matching (copy number variations): Gene > function/oncominevariantclass.","\n")
      cat("[NCI-MATCH] Trial criteria indicated by age group, variant pathogenicity, nonhotspot rules, IHC results, comments, and disease exclusions need to be manually assessed.","\n","\n")
      
      print(DF_patient_CNV[order(DF_patient_CNV$Output), ], row.names = TRUE)
    } else {
      cat("The STAMP assay did not identify any copy number variations in the solid tumor biopsied from the patient.","\n")
    }
    cat("\n","\n")
    
    if (isTRUE(continue_Fusion)) {
      cat("[Stanford OnCore] Order of matching (fusions): biomarker (gene > condition).","\n")
      cat("[Stanford OnCore] Trial criteria indicated by age group, variant pathogenicity, disease group, disease site and comments need to be manually assessed.","\n","\n")
      
      cat("[NCI-MATCH Inclusion Variants] Order of matching (fusions): Gene > Variant Type (ie. amplification, deletion).","\n")
      cat("[NCI-MATCH] Trial criteria indicated by age group, variant pathogenicity, nonhotspot rules, IHC results, comments, and disease exclusions need to be manually assessed.","\n","\n")
      
      print(DF_patient_Fusion[order(DF_patient_Fusion$Output), ], row.names = TRUE)
    } else {
      cat("The STAMP assay did not identify any fusions in the solid tumor biopsied from the patient.")
    }
    cat("\n","\n","\n")
    
    cat(email.comment,"\n")
  }
  if (exists("DF_patient_SNVIndel")) {remove(DF_patient_SNVIndel)}
  if (exists("DF_patient_CNV")) {remove(DF_patient_CNV)}
  if (exists("DF_patient_Fusion")) {remove(DF_patient_Fusion)}
  remove(patient_id,continue_SNVIndel,continue_CNV,continue_Fusion,
         colnames.SNVIndel,colnames.CNV,colnames.Fusion,email.comment)
  
  sink()
}
remove(patient_num)
