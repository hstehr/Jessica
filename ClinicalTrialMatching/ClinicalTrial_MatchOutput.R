# Generate graphics illustrating matching rate
# Two outputs generated per patient
# Detailed match results (Output.txt) and candidate trial INFO (PositiveMatch.txt)

cat(paste("Timestamp of matching output START: ", Sys.time(), sep=""),"\n")

## ILLUSTRATION PLOTs: Matching rate
#---------------------------------------------- 
# Determine which patient_id-Arm_Name combo are duplicated
# Compare results from Inclusion_Variants and Inclusion_NonHotspot
DF_NCI_Matched_Dup <- rbind(unique(data.frame(id = DF_Output_Patient_Variant$PatientID,
                                              arm = DF_Output_Patient_Variant$Arm_Name,
                                              trial = "NCI.Inclusion.Variants", stringsAsFactors = FALSE)),
                            unique(data.frame(id = DF_Output_Patient_NonHotspot$PatientID,
                                              arm = DF_Output_Patient_NonHotspot$Arm_Name,
                                              trial = "NCI.Inclusion.NonHotspot", stringsAsFactors = FALSE)))
DF_NCI_Matched_Dup <- DF_NCI_Matched_Dup$id[which(duplicated(DF_NCI_Matched_Dup[,1:2]))]

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
  
  Trial_Count$No.Internal.Trials[row_No] <- 
    as.numeric(length(unique(DF_Output_OnCore_Biomarker$OnCore.No
                             [which(DF_Output_OnCore_Biomarker$PatientID == pt_id)])))
  
  Trial_Count$No.NCI.Trials[row_No] <- 
    as.numeric(length(unique(DF_Output_Patient_Variant$Arm_Name
                             [which(DF_Output_Patient_Variant$PatientID == pt_id)])))
  
  Trial_Count$No.NCI.NonHotspot[row_No] <- 
    as.numeric(length(unique(DF_Output_Patient_NonHotspot$Arm_Name
                             [which(DF_Output_Patient_NonHotspot$PatientID == pt_id)])))
  
  # Remove duplicated patient_id-Arm_Name combo for NCI-MATCH trials (Variant and NonHotspot matches)
  if (Trial_Count$patient_id[row_No] %in% DF_NCI_Matched_Dup) {
    Trial_Count$No.NCI.NonHotspot[row_No] <- Trial_Count$No.NCI.NonHotspot[row_No] -1
  }
}
Trial_Count$No.Total.Trials <- rowSums(Trial_Count[,3:5])

## FUNCTION: save histogram plot in local computer
plot_Trial_Count <- function(DF, trialColumn, fileName_pre, plotTitle, y_max) {
  
  # Plot parameters 
  total_matched <- as.numeric(length(which(DF[[trialColumn]] != 0)))
  x_max = as.numeric(max(DF[[trialColumn]]))
  if (x_max > 2) { image_width = x_max +2
  } else { image_width = x_max +3
  }
  
  tiff(filename = paste(outdir_int, fileName_pre, "_", Sys.Date(), ".tiff", sep=""),
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

plot_Trial_Count(DF = Trial_Count, y_max = 2200,
                 trialColumn = "No.Total.Trials", 
                 fileName_pre = paste("Total_MatchDistribution_Internal_", OnCore_Biomarker_Report_timestamp, 
                                      "_NCI-MATCH_", Patient_Variant_Report_timestamp, "_", filterName, sep=""), 
                 plotTitle = "NCI-MATCH & Stanford Internal Clinical Trial Matching") 

plot_Trial_Count(DF = Trial_Count, y_max = 2200,
                 trialColumn = "No.Internal.Trials", 
                 fileName_pre = paste("OnCore_Biomarker_MatchDistribution_", OnCore_Biomarker_Report_timestamp, "_", 
                                      filterName, sep=""), 
                 plotTitle = "Stanford Internal Clinical Trial Matching") 

plot_Trial_Count(DF = Trial_Count, y_max = 2200,
                 trialColumn = "No.NCI.Trials", 
                 fileName_pre = paste("Patient_Variant_MatchDistribution_",
                                      Patient_Variant_Report_timestamp, "_", 
                                      gsub(".*pa","pa",filterName), sep=""), 
                 plotTitle = "NCI-MATCH Trial Matching (Variants)") 

plot_Trial_Count(DF = Trial_Count, y_max = 2200,
                 trialColumn = "No.NCI.NonHotspot", 
                 fileName_pre = paste("Patient_NonHotspot_MatchDistribution_",
                                      Patient_Variant_Report_timestamp, "_", 
                                      gsub(".*pa","pa",filterName), sep=""), 
                 plotTitle = "NCI-MATCH Trial Matching (NonHotspot)") 

remove(Trial_Count,DF_NCI_Matched_Dup,pt_id,row_No,plot_Trial_Count)


## GENERATE detailed match results per run
#---------------------------------------------- 
# Iterate through each patient_id of patient.list
for (patient_num  in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(outdir_patient, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir_patient, patient_id, "_Output.txt", sep=""), 
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
      "\n", "\t", paste("Specimen site: ", tolower(DF_patient$SpecimenSite[1]), sep=""), 
      "\n", "\t", paste("Dx: ", tolower(DF_patient$DxSummary[1]), sep=""), "\n")
  
  DF_patient$Output <- paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSProtein, 
                             " (",DF_patient$VariantPathogenicityStatus, ")", sep="")
  
  # Modify colnames for Output
  colnames(DF_patient)[which(colnames(DF_patient) == 
                               c("OnCore_Report_Status","Patient_Variant_Inclusion_Status",
                                 "Patient_Variant_NonHotspot_Status","Output"))] <- 
    c("OnCore_Report_Match","Inclusion_Variant_Match",
      "Inclusion_NonHotspot_Match","Variant_Identified_in_Patient")
  
  cat("\n","\n")
  cat("Stanford Internal Clinical Trials:","\n",
      "Matching schema: biomarker (gene > condition > detail) variant pathogenicity > patient age (adult:  18+ yo) >",
      " disease category (disease group > disease site).", "\n","\n")
  print(DF_patient[order(DF_patient$VariantGene),
                   c("Variant_Identified_in_Patient","OnCore_Report_Match")], row.names = FALSE)
  cat("\n")
  cat("NOTE: Trial criteria indicated by 'Disease group', 'Disease site' and 'Comments' need to be manually assessed.", "\n")
  
  cat("\n","\n")
  cat("NCI-MATCH clinical trials:", "\n",
      "Matching schema: gene > variant type > protein name (inclusion variants only) variant pathogenicity > ",
      "patient age (adult: 18+ yo) > exclusion variants identified (gene name > variant type  > protein name).", "\n","\n")
  print(DF_patient[order(DF_patient$VariantGene),
                   c("Variant_Identified_in_Patient","Inclusion_Variant_Match",
                     "Inclusion_NonHotspot_Match")], row.names = FALSE)
  cat("\n")
  cat("NOTE: Trial criteria i.e. NonHotspot Rules, IHC Results, Comments, and Disease Exclusion Codes need to be manually assessed.","\n","\n","\n")
  cat(paste("OnCore Biomarker Report last updated on ", OnCore_Biomarker_Report_timestamp, 
            " ; Patient Variant Report last updated on ", Patient_Variant_Report_timestamp, 
            ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
  sink()
}


## GENERATE trial notification for positive candidates
#---------------------------------------------- 
patient.list.matched <- sort(unique(append(DF_Output_OnCore_Biomarker$PatientID,
                                           append(DF_Output_Patient_Variant$PatientID,
                                                  DF_Output_Patient_NonHotspot$PatientID))))

## FUNCTION: print additional ARM criteria for NCI-MATCH trials
Arm_INFO <- function(trial_id) {
  cat("Trial criteria to be manually assessed: ", "\n")  
  
  # Exclusion_NonHotspot_Rules matched by trial_id
  #---------------------------------------------- 
  Exclusion_NonHotspot_Output <- 
    Exclusion_NonHotspot_Rules[Exclusion_NonHotspot_Rules$Arm_Name == trial_id,]
  
  ## Identify overlap in genes between DF_patient and Exclusion_NonHotspot_Output
  #---------------------------------------------- 
  if (nrow(Exclusion_NonHotspot_Output) > 0) {
    
    genes.patient <- sort(unique(DF_patient$VariantGene))
    exclusion.gene.match <- 
      genes.patient[which(genes.patient %in% unique(Exclusion_NonHotspot_Output$Gene_Name))]
    
    # Remove genes not relevant in Output
    Exclusion_NonHotspot_Output <- Exclusion_NonHotspot_Output[Exclusion_NonHotspot_Output$Gene_Name %in% exclusion.gene.match, ]
    
    if (nrow(Exclusion_NonHotspot_Output) > 0) {
      
      #---------------------------------------------- 
      cat("Exclusion NonHotspot criteria: ", "\n")
      Output.trial <- unique(paste(Exclusion_NonHotspot_Output$Gene_Name, ": Exon ", Exclusion_NonHotspot_Output$Exon, 
                                   " (Function: ", Exclusion_NonHotspot_Output$Function, 
                                   "; Variant Class: ", Exclusion_NonHotspot_Output$oncominevariantclass, ")", sep=""))
      for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\t") }
      cat("\n")
      
      ## Output applicable STAMP entries (DF_patient)
      #---------------------------------------------- 
      ex.gene_num <- which(DF_patient$VariantGene %in% exclusion.gene.match)
      Output.patient <- unique(paste(DF_patient$VariantGene[ex.gene_num], ": ", DF_patient$VariantHGVSProtein[ex.gene_num], 
                                     " (", DF_patient$VariantPathogenicityStatus[ex.gene_num], ")", sep=""))
      
      if (isTRUE(exists("Output.patient.IncVar"))) {
        # Remove patient variants previously printed
        Output.patient <- Output.patient[!which(Output.patient %in% Output.patient.IncVar)]  
      }
      
      if (length(Output.patient) > 0) {
        cat("Additional relevant mutations identified in patient: ","\n")
        for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n") }
      }
      cat("\n")  
    }
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
    cat("Disease_Exclusion_Codes:","\n")
    
    cat("CTEP.CATEGORY:", "\n")
    Output.trial <- sort(unique(Disease_Exclusion_Output$CTEP.CATEGORY))
    for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
    
    cat("CTEP.SUBCATEGORY:", "\n")
    Output.trial <- sort(unique(Disease_Exclusion_Output$CTEP.SUBCATEGORY))
    for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
    
    cat("CTEP.TERM:", "\n")
    Output.trial <- sort(unique(Disease_Exclusion_Output$SHORT.NAME))
    for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num], "\n") }
  }
}

# Iterate through each patient_id of patient.list.matched
patient_num = which(patient.list.matched == "TRF-950")
for (patient_num  in 1:length(patient.list.matched)) {
  patient_id <- patient.list.matched[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(outdir_patient, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  ## Write output to file (detailed match results per variant)
  #----------------------------------------------
  sink(file = paste(outdir_patient, patient_id, "_PositiveMatch.txt", sep=""), 
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
      "\n", "\t", paste("Specimen site: ", tolower(DF_patient$SpecimenSite[1]), sep=""), 
      "\n", "\t", paste("Dx: ", tolower(DF_patient$DxSummary[1]), sep=""), "\n")
  cat("\n","\n")
  
  # Extract candidate trials based on patient_id per cohort 
  #----------------------------------------------
  Matched_Internal <- 
    unique(DF_Output_OnCore_Biomarker[which(DF_Output_OnCore_Biomarker$PatientID == patient_id),])
  Matched_NCI.Variants <- 
    unique(DF_Output_Patient_Variant[which(DF_Output_Patient_Variant$PatientID == patient_id),])
  Matched_NCI.NonHotspot <- 
    unique(DF_Output_Patient_NonHotspot[which(DF_Output_Patient_NonHotspot$PatientID == patient_id),])
  
  # Output trial INFO from Matched_Internal
  #----------------------------------------------
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
  } else {
    num_internal = 0
  }
  
  # Output trial INFO from Matched_NCI.Variants
  #----------------------------------------------
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
                                     ": Exon ", DF_patient_trial_NonHotspot$Exon, 
                                     " (Function: ", DF_patient_trial_NonHotspot$Function, 
                                     "; Variant Class: ", DF_patient_trial_NonHotspot$oncominevariantclass, ")", sep=""))
        for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n") }
        
        # Remove output entry from Matched_NCI.NonHotspot
        Matched_NCI.NonHotspot <- 
          Matched_NCI.NonHotspot[!(Matched_NCI.NonHotspot$PatientID == patient_id & 
                                     Matched_NCI.NonHotspot$Arm_Name == trial_id),]
      }
      
      cat("\n")
      cat("Relevant mutations identified in patient:", "\n")
      Output.patient.IncVar <- unique(paste(DF_patient_trial$VariantGene, ": ",
                                     DF_patient_trial$VariantHGVSProtein, " (", 
                                     DF_patient_trial$VariantPathogenicityStatus, ")", sep=""))
      assign("Output.patient.IncVar", Output.patient.IncVar, envir = .GlobalEnv)
      for (entry_num in 1:length(Output.patient.IncVar)) { cat("\t", Output.patient.IncVar[entry_num],"\n") }
      
      cat("\n")
      Arm_INFO(trial_id = trial_id)
      cat("\n","\n")
    }
  } else {
    num_variant = 0
  }
  
  # Output trial INFO from Matched_NCI.NonHotspot
  #----------------------------------------------
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
      Output.trial <- unique(paste(DF_patient_trial$Gene_Name, ": Exon ", DF_patient_trial$Exon, 
                                   " (Function: ", DF_patient_trial$Function, 
                                   "; Variant Class: ", DF_patient_trial$oncominevariantclass, ")", sep=""))
      for (entry_num in 1:length(Output.trial)) { cat("\t", Output.trial[entry_num],"\n") }
      
      cat("\n")
      cat("Relevant mutations identified in patient: ","\n")
      Output.patient <- unique(paste(DF_patient_trial$Gene_Name[entry_num], ": ", 
                                     DF_patient_trial$VariantHGVSProtein[entry_num], 
                                     " (", DF_patient_trial$VariantPathogenicityStatus[entry_num], ")", sep=""))
      for (entry_num in 1:length(Output.patient)) { cat("\t", Output.patient[entry_num],"\n") }
      
      cat("\n")
      Arm_INFO(trial_id = trial_id)
      cat("\n","\n")
    } 
  }
  
  cat("NOTES:","\n")
  cat("- Stanford Internal trials: Criteria indicated by 'Disease group', 'Disease site' and 'Comments' need to be manually assessed.","\n")
  cat("- NCI-MATCH trials: Criteria i.e. NonHotspot Rules, IHC Results, Comments, and Disease Exclusion Codes need to be manually assessed.","\n")
  cat("\n")
  cat(paste("OnCore Biomarker Report last updated on ", OnCore_Biomarker_Report_timestamp, 
            "; Patient Variant Report last updated on ", Patient_Variant_Report_timestamp, 
            ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
  sink()
}

remove(DF_patient,DF_patient_trial,DF_patient_trial_NonHotspot,Matched_Internal,
       Matched_NCI.NonHotspot,Matched_NCI.Variants,entry_num,match.internal.list,
       match.nonhotspot.list,match.variant.list,num_internal,num_nonhotspot,
       num_variant,Output.patient,Output.patient.IncVar,Output.trial,patient_id,
       patient_num,patient.list.matched,trial_id,Arm_INFO)

cat(paste("Timestamp of matching output FINISH: ", Sys.time(), sep=""),"\n","\n")
