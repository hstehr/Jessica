## GENERATE detailed match results per patient

cat(paste("Timestamp of algorithm matching output START: ", Sys.time(), sep=""),"\n","\n")

# Iterate through each patient_id of patient.list
for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient
  #---------------------------------------------- 
  DF_patient_SNVIndel <- read.csv(file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                                  header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  # Specify empty dataframes
  #---------------------------------------------- 
  if (nrow(DF_patient_SNVIndel) > 0) {continue_SNVIndel <- as.logical("TRUE")
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  # Specify relevant columns depending on trial matches
  #---------------------------------------------- 
  colnames.SNVIndel <- c("Output")
  
  email.comment <- paste("This email was generated on ", Sys.time(), ".",sep="")
  
  if (isTRUE(Internal_match)) {
    colnames.SNVIndel <- append(colnames.SNVIndel, c("OnCore_SNVIndel_Status"))
    
    email.comment <- paste("OnCore Biomarker Report updated on ", OnCore_Biomarker_Report_timestamp, ". ",
                           email.comment, sep="")
  }
  
  if (isTRUE(NCI_match)) {
    colnames.SNVIndel <- append(colnames.SNVIndel, c("NCI_SNVIndel_Variant_Status","NCI_SNVIndel_NonHotspot_Status"))
    
    email.comment <- paste("Patient Variant Report updated on ", Patient_Variant_Report_timestamp, ". ",
                           email.comment, sep="")
  }
  
  # Format dataframe for output
  #---------------------------------------------- 
  if (isTRUE(continue_SNVIndel)) {
    # Specify gene info for output
    DF_patient_SNVIndel$Output <- paste(DF_patient_SNVIndel$VariantGene, " ", DF_patient_SNVIndel$VariantHGVSProtein, sep="")
    
    # Extract relevant columns for output
    DF_patient_SNVIndel <- DF_patient_SNVIndel[,colnames.SNVIndel]
  }
  
  #----------------------------------------------
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir, patient_id, ".Output.txt", sep=""),
       append = FALSE, split = FALSE)
  options(max.print=999999)
  
  # Output patient bio
  DF_patient_INFO <- STAMP_DF[which(STAMP_DF$PatientID == patient_id),]
  
  cat(paste("Ordering Physician: ",DF_patient_INFO$AssayOrderingPhysician[1], sep=""), "\n","\n")
  cat(paste("Patient ", patient_id, 
            " may qualify for the following biomarker-based clinical trial(s) based on NGS of the specimen biopsied. ",
            "Mutations were identified in ", DF_patient_INFO$AssayName[1], 
            " assay (reviewed ", DF_patient_INFO$AssayReportDateReviewed[1], ").", sep=""), 
      "\n", "\t", paste("Demographic: ", DF_patient_INFO$PatientAge[1], "yo ", DF_patient_INFO$PatientGender[1], sep=""), 
      "\n", "\t", paste("Primary tumor site: ", DF_patient_INFO$PrimaryTumorSite[1], sep=""), 
      "\n", "\t", paste("Histological Dx: ", tolower(DF_patient_INFO$HistologicalDx[1]), sep=""), "\n","\n")
  
  remove(DF_patient_INFO)
  
  cat("Feature matching algorithm applied to following clinical trials: Stanford OnCore and NCI-MATCH","\n","\n")
  
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
    
    cat(email.comment,"\n")
  }
  remove(patient_id,DF_patient_SNVIndel,continue_SNVIndel,colnames.SNVIndel,email.comment)
  
  sink()
}
remove(patient_num)
