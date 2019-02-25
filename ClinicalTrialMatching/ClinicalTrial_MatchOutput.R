# Detailed match results (Output.txt) and candidate trial INFO (PositiveMatch.txt)
# Corresponding tsv files (patient_id.OnCore.tsv and patient_id.NCI.tsv)

cat(paste("Timestamp of matching output START: ", Sys.time(), sep=""),"\n","\n")

## GENERATE detailed match results per run
#---------------------------------------------- 
# Iterate through each patient_id of patient.list
for (patient_num  in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(tempdir, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir, patient_id, ".Output.txt", sep=""), 
       append = FALSE, split = FALSE)
  options(max.print=999999)
  
  # Output patient bio
  cat(paste("Ordering Physician: ",DF_patient$AssayOrderingPhysician[1], sep=""), "\n","\n")
  cat(paste("Patient ", patient_id, " may qualify for the following clinical trial(s) based on ", 
            DF_patient$AssayType[1], " of the specimen biopsied on ", DF_patient$SpecimenDateCollected[1], 
            ". Mutations were identified in ", DF_patient$AssayName[1], 
            " assay (reviewed ", DF_patient$AssayReportDateReviewed[1], ").", sep=""), 
      "\n", "\t", paste("Demographic: ", DF_patient$PatientAge[1], "yo ", DF_patient$PatientGender[1], sep=""), 
      "\n", "\t", paste("Primary tumor site: ", DF_patient$PrimaryTumorSite[1], sep=""), 
      "\n", "\t", paste("Histological Dx: ", tolower(DF_patient$HistologicalDx[1]), sep=""), "\n")
  
  DF_patient$Output <- paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSProtein, 
                             " (",DF_patient$VariantPathogenicityStatus, ")", sep="")
  
  colnames.exist <- c("")
  colnames.new <- c("")
  
  if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
    colnames.exist <- append(colnames.exist, "OnCore_Report_Status")
    colnames.new <- append(colnames.new, "OnCore_Report_Match")
  }
  if (isTRUE(exists("DF_Output_Patient_Variant"))) {
    colnames.exist <- append(colnames.exist, "Patient_Variant_Inclusion_Status")
    colnames.new <- append(colnames.new, "Inclusion_Variant_Match")
  }
  if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
    colnames.exist <- append(colnames.exist, "Patient_Variant_NonHotspot_Status")
    colnames.new <- append(colnames.new, "Inclusion_NonHotspot_Match")
  }
  
  colnames.exist <- append(colnames.exist, "Output")
  colnames.new <- append(colnames.new, "Variant_Identified_in_Patient")
  
  colnames.exist <- colnames.exist[colnames.exist != ""]
  colnames.new <- colnames.new[colnames.new != ""]
  
  # Modify colnames for Output
  colnames(DF_patient)[which(colnames(DF_patient) %in% colnames.exist == TRUE)] <- colnames.new
  
  if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
    cat("\n","\n")
    cat("Stanford Internal Clinical Trials:","\n",
        "Matching schema: biomarker (gene > condition > detail) > variant pathogenicity > patient age (adult:  18+ yo) >",
        " disease category (disease group > disease site) based on primary tumor site.", "\n","\n")
    print(DF_patient[order(DF_patient$VariantGene),
                     c("Variant_Identified_in_Patient","OnCore_Report_Match")], row.names = FALSE)
    cat("\n")
    cat("NOTE: Trial criteria indicated by 'Disease group', 'Disease site' and 'Comments' need to be manually assessed.", "\n")
  }
  
  if (isTRUE(exists("DF_Output_Patient_Variant") | exists("DF_Output_Patient_NonHotspot"))) {
    cat("\n","\n")
    cat("NCI-MATCH clinical trials:", "\n",
        "Matching schema: gene > variant type > protein name (inclusion variants only) > ",
        "exon number (inclusion nonhotspots only) > variant pathogenicity > ",
        "patient age (adult: 18+ yo) > exclusion variants identified (gene name > variant type  > protein name) ", 
        "> exclusion nonhotspots identified (gene name > variant type  > exon number) ",
        "> disease exclusion (CTEP categories) based on histological dx and primary tumor site.","\n","\n")
    
    if (isTRUE(exists("DF_Output_Patient_Variant") | exists("DF_Output_Patient_NonHotspot"))) {
      print(DF_patient[order(DF_patient$VariantGene),
                       c("Variant_Identified_in_Patient","Inclusion_Variant_Match",
                         "Inclusion_NonHotspot_Match")], row.names = FALSE)
    } else if (isTRUE(exists("DF_Output_Patient_Variant"))) {
      print(DF_patient[order(DF_patient$VariantGene),
                       c("Variant_Identified_in_Patient","Inclusion_Variant_Match")], row.names = FALSE)
    } else if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
      print(DF_patient[order(DF_patient$VariantGene), 
                       c("Variant_Identified_in_Patient","Inclusion_NonHotspot_Match")], row.names = FALSE)
    }
    
    cat("\n")
    cat("NOTE: Trial criteria i.e. NonHotspot Rules, IHC Results, Comments, and Disease Exclusion Codes need to be manually assessed.","\n","\n","\n")
    cat(paste("OnCore Biomarker Report updated on ", OnCore_Biomarker_Report_timestamp, 
              " and Patient Variant Report updated on ", Patient_Variant_Report_timestamp, 
              ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
  }
  sink()
}


## GENERATE trial notification for positive candidates
#---------------------------------------------- 
if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
  patient.oncore.matched <- sort(unique(DF_Output_OnCore_Biomarker$PatientID))
}

if (isTRUE(exists("DF_Output_Patient_Variant") & exists("DF_Output_Patient_NonHotspot"))) {
  patient.nci.matched <- sort(unique(append(DF_Output_Patient_Variant$PatientID,
                                            DF_Output_Patient_NonHotspot$PatientID)))
} else if (isTRUE(exists("DF_Output_Patient_Variant"))) {
  patient.nci.matched <- sort(unique(DF_Output_Patient_Variant$PatientID))
} else if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
  patient.nci.matched <- sort(unique(DF_Output_Patient_NonHotspot$PatientID))
}

if (isTRUE(exists("DF_Output_OnCore_Biomarker") &
           exists("DF_Output_Patient_Variant") & exists("DF_Output_Patient_NonHotspot"))) {
  patient.list.matched <- sort(unique(append(patient.oncore.matched,patient.nci.matched)))
} else if (isTRUE(exists("DF_Output_Patient_Variant") & exists("DF_Output_Patient_NonHotspot"))) {
  patient.list.matched <- patient.nci.matched
} else if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
  patient.list.matched <- patient.oncore.matched
}

## FUNCTION: print additional ARM criteria for NCI-MATCH trials
Arm_INFO <- function(trial_id) {
  cat("Trial criteria to be manually assessed: ", "\n")  
  
  # Exclusion_NonHotspot_Rules matched by trial_id
  #---------------------------------------------- 
  Exclusion_NonHotspot_Output <- 
    Exclusion_NonHotspot_Rules[Exclusion_NonHotspot_Rules$Arm_Name == trial_id,]
  
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
  }
  
  # Comments matched by trial_id
  #---------------------------------------------- 
  Comments_Output <- Comments[Comments$Arm_Name == trial_id,]
  
  if (nrow(Comments_Output) > 0) {
    
    ## Output applicable Comments criteria
    #---------------------------------------------- 
    cat("Comments:", "\n")
    print(Comments_Output, row.names = FALSE)
    cat("\n")
  }
  
  # Disease_Exclusion_Codes matched by trial_id
  #---------------------------------------------- 
  Disease_Exclusion_Output <- Disease_Exclusion_Codes[Disease_Exclusion_Codes$Arm_Name == trial_id,]
  
  if (nrow(Disease_Exclusion_Output) > 0) {
    
    ## Output applicable Disease_Exclusion_Codes criteria
    #---------------------------------------------- 
    cat("Histological disease exclusion codes:","\n")
    Output.trial <- sort(unique(Disease_Exclusion_Output$SHORT.NAME))
    for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
  }
}

# Iterate through each patient_id of patient.list.matched
for (patient_num  in 1:length(patient.list.matched)) {
  patient_id <- patient.list.matched[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(tempdir, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir, patient_id, ".PositiveMatch.txt", sep=""), 
       append = FALSE, split = FALSE)
  options(max.print=999999)
  
  # Output patient bio
  cat(paste("Ordering Physician: ",DF_patient$AssayOrderingPhysician[1], sep=""), "\n","\n")
  cat(paste("Patient ", patient_id, " may qualify for the following clinical trial(s) based on ", 
            DF_patient$AssayType[1], " of the specimen biopsied on ", DF_patient$SpecimenDateCollected[1], 
            ". Mutations were identified in ", DF_patient$AssayName[1], 
            " assay (reviewed ", DF_patient$AssayReportDateReviewed[1], ").", sep=""), 
      "\n", "\t", paste("Demographic: ", DF_patient$PatientAge[1], "yo ", DF_patient$PatientGender[1], sep=""), 
      "\n", "\t", paste("Primary tumor site: ", DF_patient$PrimaryTumorSite[1], sep=""), 
      "\n", "\t", paste("Histological Dx: ", tolower(DF_patient$HistologicalDx[1]), sep=""), "\n")
  cat("\n","\n")
  
  # Extract candidate trials based on patient_id per cohort 
  #----------------------------------------------
  if (isTRUE(exists("DF_Output_OnCore_Biomarker"))) {
    Matched_Internal <- 
      unique(DF_Output_OnCore_Biomarker[which(DF_Output_OnCore_Biomarker$PatientID == patient_id),])
  }
  
  if (isTRUE(exists("DF_Output_Patient_Variant"))) {
    Matched_NCI.Variants <- 
      unique(DF_Output_Patient_Variant[which(DF_Output_Patient_Variant$PatientID == patient_id),])
  }
  
  if (isTRUE(exists("DF_Output_Patient_NonHotspot"))) {
    Matched_NCI.NonHotspot <- 
      unique(DF_Output_Patient_NonHotspot[which(DF_Output_Patient_NonHotspot$PatientID == patient_id),])
  }
  
  # Output trial INFO from Matched_Internal
  #----------------------------------------------
  if (isTRUE(exists("Matched_Internal"))) {
    if (nrow(Matched_Internal) > 0 ) {
      
      # Extract list of matched OnCore.No
      match.internal.list <- sort(unique(Matched_Internal$OnCore.No))
      
      # Iterate through each matched trial_id per patient 
      #----------------------------------------------
      for (num_internal in 1:length(match.internal.list)) {
        trial_id <- match.internal.list[num_internal]
        
        # Extract trial INFO based on trial_id match in Matched_Internal
        DF_patient_trial <- Matched_Internal[Matched_Internal$OnCore.No == trial_id, ]
        
        # Output notification
        #----------------------------------------------
        cat(paste("No.", num_internal, ": Clinical Trial #", trial_id, ": ", DF_patient_trial$Title, sep=""),"\n",
            "----------------------------------------------------------------------", "\n")
        cat(paste("Principal Investigator: ", DF_patient_trial$PI, " (", DF_patient_trial$PI.Email, ")", sep=""),"\n")
        cat(paste("Primary Clinical Research Coordinator: ", DF_patient_trial$Primary.CRC, 
                  " (", DF_patient_trial$Primary.CRC.Email, ")", sep=""),"\n","\n")
        cat("Biomarker criteria:", "\n", "\t",DF_patient_trial$Biomarker.Description,"\n")
        cat(paste("Disease group criteria: ", DF_patient_trial$Disease.Group, sep=""),"\n")
        cat(paste("Disease site criteria: ", DF_patient_trial$Disease.Sites, sep=""),"\n","\n")
        
        cat("Relevant mutations identified in patient:", "\n")
        Output.patient <-  unique(paste(DF_patient_trial$VariantGene, ": ", DF_patient_trial$VariantHGVSProtein, 
                                        " (", DF_patient_trial$VariantPathogenicityStatus, ")", sep=""))
        for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n")}
        cat("\n","\n")
      }
      
      remove(Matched_Internal)
    } 
  } else {
    num_internal = 0
  }
  
  # Output trial INFO from Matched_NCI.Variants
  #----------------------------------------------
  if (isTRUE(exists("Matched_NCI.Variants"))) {
    if (nrow(Matched_NCI.Variants) > 0) {
      
      # Extract list of matched Arm_No
      match.variant.list <- sort(unique(Matched_NCI.Variants$Arm_Name))
      
      # Iterate through each matched trial_id per patient 
      #----------------------------------------------
      for (num_variant in 1:length(match.variant.list)) {
        trial_id <- match.variant.list[num_variant]
        
        # Extract trial INFO based on trial_id match in Matched_NCI.Variants
        DF_patient_trial <- Matched_NCI.Variants[Matched_NCI.Variants$Arm_Name == trial_id, ]
        
        # Output notification
        #----------------------------------------------
        cat(paste("No.", (num_internal + num_variant), ": NCI-MATCH Trial Treatment ", trial_id, sep=""),"\n", 
            "----------------------------------------------------------------------", "\n")
        cat("Inclusion criteria: ", "\n")
        # Several matches with same HGVS protein nomenclature but different CHR:POS_REF_ALT nomenclature
        Output.trial <- unique(paste(DF_patient_trial$Gene_Name, ": ", DF_patient_trial$Protein, sep=""))
        for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n") }
        
        # Extract trial INFO based on trial_id match in Matched_NCI.NonHotspot (if applicable)
        DF_patient_trial_NonHotspot <- Matched_NCI.NonHotspot[Matched_NCI.NonHotspot$Arm_Name == trial_id, ]
        
        # Output notification
        # For (patient_id/trial_id) overlap between Matched_NCI.Variants and Matched_NCI.NonHotspot
        #----------------------------------------------
        if (nrow(DF_patient_trial_NonHotspot) > 0) {
          cat("Inclusion NonHotspot criteria: ", "\n")
          Output.trial <- unique(paste(DF_patient_trial_NonHotspot$Gene_Name, 
                                       " Exon ", DF_patient_trial_NonHotspot$Exon, 
                                       " (", DF_patient_trial_NonHotspot$Function,")", sep=""))
          for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n") }
          
          # Remove output entry from Matched_NCI.NonHotspot
          Matched_NCI.NonHotspot <- 
            Matched_NCI.NonHotspot[!(Matched_NCI.NonHotspot$PatientID == patient_id & 
                                       Matched_NCI.NonHotspot$Arm_Name == trial_id),]
        }
        
        cat("\n")
        cat("Relevant mutations identified in patient:", "\n")
        Output.patient.IncVar <- unique(paste(DF_patient_trial$VariantGene, " Exon ",
                                              DF_patient_trial$Exon_Number, ": ",
                                              DF_patient_trial$VariantHGVSProtein, " (", 
                                              DF_patient_trial$VariantPathogenicityStatus, ")", sep=""))
        assign("Output.patient.IncVar", Output.patient.IncVar, envir = .GlobalEnv)
        for (entry_num in 1:length(Output.patient.IncVar)) { cat("\t", Output.patient.IncVar[entry_num],"\n") }
        
        cat("\n")
        Arm_INFO(trial_id = trial_id)
        cat("\n","\n")
      }
    }
  } else {
    num_variant = 0
  }
  
  # Output trial INFO from Matched_NCI.NonHotspot
  #----------------------------------------------
  if (isTRUE(exists("Matched_NCI.NonHotspot"))) {
    if (nrow(Matched_NCI.NonHotspot) > 0 ) {
      
      # Extract list of matched Arm_No
      match.nonhotspot.list <- sort(unique(Matched_NCI.NonHotspot$Arm_Name))
      
      # Iterate through each matched trial_id per patient 
      #----------------------------------------------
      for (num_nonhotspot in 1:length(match.nonhotspot.list)) {
        trial_id <- match.nonhotspot.list[num_nonhotspot]
        
        # Extract trial INFO based on trial_id match in Matched_NCI.Variants
        DF_patient_trial <- Matched_NCI.NonHotspot[Matched_NCI.NonHotspot$Arm_Name == trial_id, ]
        
        # Output notification
        #----------------------------------------------
        cat(paste("No.", (num_internal + num_variant + num_nonhotspot), ": NCI-MATCH Trial Treatment ", trial_id, sep=""),"\n", 
            "----------------------------------------------------------------------", "\n")
        cat("Inclusion NonHotspot criteria: ", "\n")
        Output.trial <- unique(paste(DF_patient_trial$Gene_Name, " Exon ", DF_patient_trial$Exon, 
                                     " (", DF_patient_trial$Function, ")", sep=""))
        for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n") }
        
        cat("\n")
        cat("Relevant mutations identified in patient: ","\n")
        Output.patient <- unique(paste(DF_patient_trial$Gene_Name[entry_num], " Exon ",
                                       DF_patient_trial$Exon_Number, ": ",
                                       DF_patient_trial$VariantHGVSProtein[entry_num], 
                                       " (", DF_patient_trial$VariantPathogenicityStatus[entry_num], ")", sep=""))
        for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n") }
        
        cat("\n")
        Arm_INFO(trial_id = trial_id)
        cat("\n","\n")
      } 
      
      remove(match.nonhotspot.list,num_nonhotspot)
    }
  }
  
  cat("NOTES:","\n")
  cat("- Stanford Internal trials: Criteria indicated by Disease group, Disease site and Comments need to be manually assessed.","\n")
  cat("- NCI-MATCH trials: Criteria i.e. IHC Results, Comments and Disease Exclusion Codes need to be manually assessed.","\n")
  cat("\n")
  cat(paste("OnCore Biomarker Report updated on ", OnCore_Biomarker_Report_timestamp, 
            " and Patient Variant Report updated on ", Patient_Variant_Report_timestamp, 
            ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
  sink()
}


## GENERATE tsv output for positive matches
#---------------------------------------------- 
# Iterate through each patient_id of patient.oncore.matched
if (isTRUE(exists("patient.oncore.matched"))) {
  for (patient_num  in 1:length(patient.oncore.matched)) {
    patient_id <- patient.oncore.matched[patient_num]
    
    # Import STAMP entries per patient
    DF_patient_Oncore <- unique(DF_Output_OnCore_Biomarker[DF_Output_OnCore_Biomarker$PatientID == patient_id,
                                                           c("PatientID","VariantGene","VariantHGVSCoding",
                                                             "VariantHGVSProtein","OnCore.No","NCT..")])
    # Rename column name
    colnames(DF_patient_Oncore)[6] <- "NCT"
    
    # Write to match results per patient to local computer
    #----------------------------------------------
    write.table(DF_patient_Oncore, file = paste(outdir, patient_id, ".OnCore.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  remove(patient.oncore.matched,DF_patient_Oncore)
}

# Iterate through each patient_id of patient.nci.matched
if (isTRUE(exists("patient.nci.matched"))) {
  for (patient_num  in 1:length(patient.nci.matched)) {
    patient_id <- patient.nci.matched[patient_num]
    
    # Import STAMP entries per patient
    DF_patient_NCI.Variants <- DF_Output_Patient_Variant[DF_Output_Patient_Variant$PatientID == patient_id,
                                                         c("PatientID","VariantGene","VariantHGVSCoding",
                                                           "VariantHGVSProtein","Arm_Name")]
    DF_patient_NCI.NonHotspot <- DF_Output_Patient_NonHotspot[DF_Output_Patient_NonHotspot$PatientID == patient_id,
                                                              c("PatientID","VariantGene","VariantHGVSCoding",
                                                                "VariantHGVSProtein","Arm_Name")]
    DF_patient_NCI <- unique(rbind(DF_patient_NCI.Variants, DF_patient_NCI.NonHotspot))
    DF_patient_NCI$NCT <- "NCT02465060"
    
    # Write to match results per patient to local computer
    #----------------------------------------------
    write.table(DF_patient_NCI, file = paste(outdir, patient_id, ".NCI.tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  remove(patient.nci.matched,DF_patient_NCI.Variants,DF_patient_NCI.NonHotspot)
}

# ## General CONFIRMATIONS
# #----------------------------------------------
# # Merge all DF_patient tsv
# DF_patient_Full <- data.frame()
# for (patient_num  in 1:length(patient.list)) {
#   patient_id <- patient.list[patient_num]
# 
#   # Import STAMP entries per patient
#   DF_patient <- read.csv(file = paste(tempdir, patient_id, ".tsv", sep=""),
#                          header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
# 
#   DF_patient_Full <- rbind(DF_patient_Full,DF_patient)
# }
# 
# # CONFIRM output generated for all entries
# DF_patient_Full$PatientID[is.na(DF_patient_Full$OnCore_Report_Status)]
# DF_patient_Full$PatientID[is.na(DF_patient_Full$Patient_Variant_Inclusion_Status)]
# DF_patient_Full$PatientID[is.na(DF_patient_Full$Patient_Variant_NonHotspot_Status)]
# 
# # CONFIRM output for "Candidate trial IDENTIFIED"
# OnCore_Confirm <- DF_patient_Full[which(DF_patient_Full$OnCore_Report_Status == "Candidate trial IDENTIFIED"),c(1:(ncol(DF_patient_Full) -2))]
# OnCore_Confirm <- full_join(DF_Output_OnCore_Biomarker,OnCore_Confirm)
# # Not indicated appropriately in DF_Output_OnCore_Biomarker
# OnCore_Confirm[which(OnCore_Confirm$OnCore_Report_Status != "Candidate trial IDENTIFIED"),]
# # Not indicated appropriately in DF_patient
# OnCore_Confirm[which(is.na(OnCore_Confirm$OnCore.No)),]
# 
# NCI.Variants_Confirm <- DF_patient_Full[which(grepl("Candidate trial IDENTIFIED", DF_patient_Full$Patient_Variant_Inclusion_Status)),]
# NCI.Variants_Confirm <- NCI.Variants_Confirm[, !(colnames(NCI.Variants_Confirm) %in% c("OnCore_Report_Status","Patient_Variant_NonHotspot_Status"))]
# NCI.Variants_Confirm <- full_join(DF_Output_Patient_Variant,NCI.Variants_Confirm)
# # Not indicated appropriately in DF_Output_Patient_Variant
# NCI.Variants_Confirm[which(grepl("Candidate trial IDENTIFIED", NCI.Variants_Confirm$Patient_Variant_Inclusion_Status) == FALSE),]
# # Not indicated appropriately in DF_patient
# NCI.Variants_Confirm[which(is.na(NCI.Variants_Confirm$Arm_Name)),]
# 
# NCI.NonHotspot_Confirm <- DF_patient_Full[which(grepl("Candidate trial IDENTIFIED", DF_patient_Full$Patient_Variant_NonHotspot_Status)),]
# NCI.NonHotspot_Confirm <- NCI.NonHotspot_Confirm[, !(colnames(NCI.NonHotspot_Confirm) %in% c("OnCore_Report_Status","Patient_Variant_Inclusion_Status"))]
# NCI.NonHotspot_Confirm <- full_join(DF_Output_Patient_NonHotspot,NCI.NonHotspot_Confirm)
# # Not indicated appropriately in Patient_Variant_NonHotspot_Status
# NCI.NonHotspot_Confirm[which(grepl("Candidate trial IDENTIFIED", NCI.NonHotspot_Confirm$Patient_Variant_NonHotspot_Status) == FALSE),]
# # Not indicated appropriately in DF_patient
# NCI.NonHotspot_Confirm[which(is.na(NCI.NonHotspot_Confirm$Arm_Name)),]
# 
# remove(DF_patient_Full,OnCore_Confirm,NCI.Variants_Confirm,NCI.NonHotspot_Confirm)

if (isTRUE(exists("Output.patient"))) {remove(Output.patient)}
remove(DF_patient,DF_patient_trial,DF_patient_trial_NonHotspot,Matched_NCI.NonHotspot,
       Matched_NCI.Variants,entry_num,match.variant.list,num_internal,num_variant,Output.trial,
       patient_id,patient_num,patient.list.matched,trial_id,Arm_INFO,DF_patient_NCI)

cat(paste("Timestamp of matching output FINISH: ", Sys.time(), sep=""),"\n")
