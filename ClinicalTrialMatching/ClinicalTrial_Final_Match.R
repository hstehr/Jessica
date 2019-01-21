setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load relevant files
#---------------------------------------------- 
STAMP_all_variants_QC <- 
  read.csv(file = paste("ClinicalTrialMatching/",Syapse_Export_timestamp, "_syapse_export_DF_STAMP_VariantAnno.csv", sep=""),
           header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")

int_file_02 = paste(getwd(), "/ClinicalTrialMatching/OnCore_Biomarker_Matched_", OnCore_Biomarker_Report_timestamp, 
                    "_",filterName_initial, filterName_int, "_FINAL.csv", sep="")
OnCore_Biomarker_Matched <- read.csv(file = int_file_02, header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

int_file_03 = paste(getwd(), "/ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_Matched_Variants_FINAL.csv", sep="")
Patient_Variant_Matched <- read.csv(file = int_file_03, header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

int_file_04 = paste(getwd(), "/ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_Matched_NonHotspot_FINAL.csv", sep="")
Patient_NonHotspot_Matched <- read.csv(file = int_file_04, header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## Additional criteria for NCI-MATCH trial ARMs
#---------------------------------------------- 
# Import DF_Exclusion_Variants
int_file_05 = paste(getwd(), "/ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_.xlsx", sep="")

PATIENT_VARIANT_REPORT_QC <- import_list(int_file_05, setclass = "tbl")
invisible(capture.output(lapply(names(PATIENT_VARIANT_REPORT_QC),
                                function(x) assign(x, PATIENT_VARIANT_REPORT_QC[[x]],
                                                   envir = .GlobalEnv))))
remove(PATIENT_VARIANT_REPORT_QC,DF_Exclusion_Variants,DF_Inclusion_NonHotspot_Rules,DF_Inclusion_Variants)

## Extract only variants with pathogenicity status indicated in "pathogenic_accepted"
#---------------------------------------------- 
if (isTRUE(pathogenic_FILTER)) {
  print("pathogenic Filter: ON")
  print("Patient_Variant_Matched")
  print(paste("Pre-FILTER: n= ", nrow(Patient_Variant_Matched), " entries (n= ",
              length(unique(Patient_Variant_Matched$sys.uniqueId)), " patients)", sep=""))
  TempDF <- Patient_Variant_Matched[!Patient_Variant_Matched$smpl.pathogenicityStatus %in% pathogenic_accepted,]
  print(paste("REMOVED: n= ", nrow(TempDF), " entries (n= ", length(unique(TempDF$sys.uniqueId)), " patients)", sep=""))
  cat("\n")

  Patient_Variant_Matched <- 
    Patient_Variant_Matched[Patient_Variant_Matched$smpl.pathogenicityStatus %in% pathogenic_accepted,]
  print(paste("POST-FILTER: n= ", nrow(Patient_Variant_Matched), " entries (n= ", 
              length(unique(Patient_Variant_Matched$sys.uniqueId)), " patients)", sep=""))
  
  Patient_NonHotspot_Matched <- 
    Patient_NonHotspot_Matched[Patient_NonHotspot_Matched$smpl.pathogenicityStatus %in% pathogenic_accepted,]
  
  print("OnCore_Biomarker_Matched")
  print(paste("Pre-FILTER: n= ", nrow(OnCore_Biomarker_Matched), " entries (n= ",
              length(unique(OnCore_Biomarker_Matched$sys.uniqueId)), " patients)", sep=""))
  TempDF <- OnCore_Biomarker_Matched[!OnCore_Biomarker_Matched$smpl.pathogenicityStatus %in% pathogenic_accepted,]
  print(paste("REMOVED: n= ", nrow(TempDF), " entries (n= ", length(unique(TempDF$sys.uniqueId)), " patients)", sep=""))
  cat("\n")
  
  OnCore_Biomarker_Matched <- 
    OnCore_Biomarker_Matched[OnCore_Biomarker_Matched$smpl.pathogenicityStatus %in% pathogenic_accepted,]
  print(paste("POST-FILTER: n= ", nrow(OnCore_Biomarker_Matched), " entries (n= ", 
              length(unique(OnCore_Biomarker_Matched$sys.uniqueId)), " patients)", sep=""))
  
  remove(TempDF)
} else {
  print("pathogenic Filter: OFF")
  print("Patient_Variant_Matched")
  print(paste("Current: n= ", nrow(Patient_Variant_Matched), " entries (n= ",
              length(unique(Patient_Variant_Matched$sys.uniqueId)), " patients)", sep=""))
  cat("\n")
  print("OnCore_Biomarker_Matched")
  print(paste("Current: n= ", nrow(OnCore_Biomarker_Matched), " entries (n= ",
              length(unique(OnCore_Biomarker_Matched$sys.uniqueId)), " patients)", sep=""))
  
}

## ILLUSTRATION PLOTs: Matching rate
#---------------------------------------------- 
## FUNCTION
# Save histogram plot in local computer
plot_Trial_Count <- function(DF, trialColumn, fileName_pre, plotTitle, y_max) {
  
  # Plot parameters 
  total_matched <- as.numeric(length(which(DF[[trialColumn]] != 0)))
  x_max = as.numeric(max(DF[[trialColumn]]))
  if (x_max > 2) {
    image_width = x_max +2
  } else {
    image_width = x_max +3
  }
  
  tiff(filename = paste(getwd(), "/ClinicalTrialMatching/Retrospective_Analysis/", 
                        fileName_pre, "_", Sys.Date(), ".tiff", sep=""),
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

# Determine which patient_id-Arm_Name combo are duplicated
# Compare results from Inclusion_Variants and Inclusion_NonHotspot
DF_NCI_Matched_Dup <- rbind(unique(data.frame(id = Patient_Variant_Matched$sys.uniqueId,
                                          arm = Patient_Variant_Matched$Arm_Name,
                                          trial = "NCI.Inclusion.Variants", stringsAsFactors = FALSE)),
                        unique(data.frame(id = Patient_NonHotspot_Matched$sys.uniqueId,
                                          arm = Patient_NonHotspot_Matched$Arm_Name,
                                          trial = "NCI.Inclusion.NonHotspot", stringsAsFactors = FALSE)))
DF_NCI_Matched_Dup <- DF_NCI_Matched_Dup$id[which(duplicated(DF_NCI_Matched_Dup[,1:2]))]

# Confirmation check of correct adjustment to "No.NCI.NonHotspot"
print(paste("Patient_id is duplicated at most ONCE for NCI-MATCH trials: ", 
      (length(DF_NCI_Matched_Dup) == length(unique(DF_NCI_Matched_Dup))),
      sep=""))
cat("\n","\n")

# Tally No. trial matches from each source
Trial_Count <- data.frame(patient_id = STAMP_all_variants_QC$sys.uniqueId,
                          patient_age = STAMP_all_variants_QC$current.age,
                          stringsAsFactors = FALSE)
Trial_Count$patient_age <- as.numeric(Trial_Count$patient_age)
Trial_Count <- unique(Trial_Count)

# Subset only ADULT patients (18+yo)
if (isTRUE(pathogenic_FILTER)) {
  Trial_Count <- Trial_Count[Trial_Count$patient_age >= 18,]
}

Trial_Count$No.Internal.Trials <- NA
Trial_Count$No.NCI.Trials <- NA
Trial_Count$No.NCI.NonHotspot <- NA

for (row_No in 1:nrow(Trial_Count)) {
  pt_id <- Trial_Count$patient_id[row_No]
  
  Trial_Count$No.Internal.Trials[row_No] <- 
    as.numeric(length(unique(OnCore_Biomarker_Matched$OnCore.No
                             [which(OnCore_Biomarker_Matched$sys.uniqueId == pt_id)])))
  
  Trial_Count$No.NCI.Trials[row_No] <- 
    as.numeric(length(unique(Patient_Variant_Matched$Arm_Name
                             [which(Patient_Variant_Matched$sys.uniqueId == pt_id)])))
  
  Trial_Count$No.NCI.NonHotspot[row_No] <- 
    as.numeric(length(unique(Patient_NonHotspot_Matched$Arm_Name
                             [which(Patient_NonHotspot_Matched$sys.uniqueId == pt_id)])))
  
  # Remove duplicated patient_id-Arm_Name combo for NCI-MATCH trials (Variant and NonHotspot matches)
  if (Trial_Count$patient_id[row_No] %in% DF_NCI_Matched_Dup) {
    Trial_Count$No.NCI.NonHotspot[row_No] <- Trial_Count$No.NCI.NonHotspot[row_No] -1
  }
  
  remove(pt_id)
}
Trial_Count$No.Total.Trials <- rowSums(Trial_Count[,3:5])

plot_Trial_Count(DF = Trial_Count, y_max = 2180,
                 trialColumn = "No.Total.Trials", 
                 fileName_pre = paste("Total_MatchDistribution_Internal", OnCore_Biomarker_Report_timestamp, 
                                      "_NCI-MATCH_", Patient_Variant_Report_timestamp, "_", filterName, sep=""), 
                 plotTitle = "NCI-MATCH & Stanford Internal Clinical Trial Matching") 

plot_Trial_Count(DF = Trial_Count, y_max = 2180,
                 trialColumn = "No.Internal.Trials", 
                 fileName_pre = paste("OnCore_Biomarker_MatchDistribution_", OnCore_Biomarker_Report_timestamp, "_", 
                                      filterName_initial, filterName_int, sep=""), 
                 plotTitle = "Stanford Internal Clinical Trial Matching") 

plot_Trial_Count(DF = Trial_Count, y_max = 2180,
                 trialColumn = "No.NCI.Trials", 
                 fileName_pre = paste("Patient_Variant_MatchDistribution_",
                                      Patient_Variant_Report_timestamp, "_", pathogenic_pre, sep=""), 
                 plotTitle = "NCI-MATCH Trial Matching (Variants)") 

plot_Trial_Count(DF = Trial_Count, y_max = 2180,
                 trialColumn = "No.NCI.NonHotspot", 
                 fileName_pre = paste("Patient_NonHotspot_MatchDistribution_",
                                      Patient_Variant_Report_timestamp, "_", pathogenic_pre, sep=""), 
                 plotTitle = "NCI-MATCH Trial Matching (NonHotspot)") 

remove(Trial_Count,row_No,DF_NCI_Matched_Dup,plot_Trial_Count)

## Generate email notification
#---------------------------------------------- 
## FUNCTION
# Print STAMP entry to screen for Exclusion_NonHotspot_Rules
Arm_Exclusion <- function(patient_DF, matched_DF, exclusion_DF) {
  patient_id <- sort(unique(matched_DF$sys.uniqueId))
  
  ## Extract all STAMP entries for each patient
  patient_DF <- patient_DF[patient_DF$sys.uniqueId == patient_id,]
  patient_gene <- unique(patient_DF$base.gene)
  
  exclusion_gene <- unique(exclusion_DF$Gene_Name)
  DF_exclusion_output <- patient_DF[patient_DF$base.gene %in% exclusion_gene,]
  DF_exclusion_output <- DF_exclusion_output[!(DF_exclusion_output$smpl.hgvsProtein %in% matched_DF$smpl.hgvsProtein),]
  
  assign("DF_exclusion_output", DF_exclusion_output, envir = .GlobalEnv)
}

# Print additional ARM criteria for NCI-MATCH trials
Arm_INFO <- function(trial_id) {
  
  # Extract rows with ARM_Name of interest
  DF_Exclusion_NonHotspot_Patient <- 
    DF_Exclusion_NonHotspot_Rules[DF_Exclusion_NonHotspot_Rules$Arm_Name == trial_id,]
  DF_IHC_Patient <- DF_IHC_Results[DF_IHC_Results$Arm_Name == trial_id,]
  DF_Comments_Patient <- DF_Comments[DF_Comments$Arm_Name == trial_id,]
  DF_Disease_Exclusion_Codes_Patient <- DF_Disease_Exclusion_Codes[DF_Disease_Exclusion_Codes$Arm_Name == trial_id,]
  
  # Determine if any rows of interested extracted 
  if (nrow(DF_Exclusion_NonHotspot_Patient) + nrow(DF_IHC_Patient) +
      nrow(DF_Comments_Patient) + nrow(DF_Disease_Exclusion_Codes_Patient) > 0) {
    cat("\n")
    cat("Trial criteria to be manually assessed: ", "\n")  
    
    # Print rows from DF_Exclusion_NonHotspot_Rules if match by Gene_Name
    if (nrow(DF_Exclusion_NonHotspot_Patient) > 0) {
      if (trial_num == 0) {
        Arm_Exclusion(patient_DF = STAMP_all_variants_QC, 
                      matched_DF = DF_patient_NCI, 
                      exclusion_DF = DF_Exclusion_NonHotspot_Patient)
      } else {
        Arm_Exclusion(patient_DF = STAMP_all_variants_QC, 
                      matched_DF = DF_patient_NCI.NonHotspot, 
                      exclusion_DF = DF_Exclusion_NonHotspot_Patient)
      }
      
      if (nrow(DF_exclusion_output) > 0) {
        cat("Exclusion NonHotspot criteria: ", "\n")
        for (entry_num in 1:nrow(DF_Exclusion_NonHotspot_Patient)) {
          cat("\t", paste(DF_Exclusion_NonHotspot_Patient$Gene_Name[entry_num], 
                          ":Exon ", DF_Exclusion_NonHotspot_Patient$Exon[entry_num], 
                          " (Function: ", DF_Exclusion_NonHotspot_Patient$Function[entry_num], 
                          "; Variant Class: ", DF_Exclusion_NonHotspot_Patient$oncominevariantclass[entry_num], 
                          ")", sep=""),"\t")
        }
        cat("\n", paste("Additional relevant mutated transcript(s) identified in patient ", patient_id, ":", sep=""),"\n")
        for (entry_num in 1:nrow(DF_exclusion_output)) {
          cat("\t", paste(DF_exclusion_output$base.gene[entry_num], ":", DF_exclusion_output$smpl.hgvsProtein[entry_num], 
                          " (", DF_exclusion_output$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
        }
        cat("\n")
      }
    }
    
    # Print rows from DF_IHC_Results
    if (nrow(DF_IHC_Patient) > 0) {
      cat("IHC results: ", "\n")
      for (entry_num in 1:nrow(DF_IHC_Patient)) {
        cat("\t", paste(DF_IHC_Patient$Gene[entry_num], 
                        "(Status: ", DF_IHC_Patient$Status_POSITIVE_NEGATIVE_INDETERMINATE[entry_num], 
                        "; Variant: ", DF_IHC_Patient$Variant_PRESENT_NEGATIVE_EMPTY[entry_num], 
                        ")", sep=""),"\t")
      }
      cat("\n")
    }
    
    # Print rows from DF_Comments
    if (nrow(DF_Comments_Patient) > 0) {
      cat("Comments:", "\n")
      print(DF_Comments_Patient, row.names = FALSE)
      cat("\n")
    }
    
    # Print rows from DF_Disease_Exclusion_Codes
    if (nrow(DF_Disease_Exclusion_Codes_Patient) > 0) {
      cat("\n")
      cat("Disease_Exclusion_Codes:", "\n")
      cat("CTEP.CATEGORY:", "\n")
      print(unique(DF_Disease_Exclusion_Codes_Patient$CTEP.CATEGORY))
      
      cat("CTEP.SUBCATEGORY:", "\n")
      print(unique(DF_Disease_Exclusion_Codes_Patient$CTEP.SUBCATEGORY))
      
      cat("CTEP.TERM:", "\n")
      print(unique(DF_Disease_Exclusion_Codes_Patient$SHORT.NAME))
      
      cat("MedDRA.CODE:", "\n")
      print(unique(DF_Disease_Exclusion_Codes_Patient$MedDRA.CODE))
      cat("\n")
    }
  }
}

# Generate output per unique patient
patient.list <- sort(unique(append(OnCore_Biomarker_Matched$sys.uniqueId, 
                                   append(Patient_Variant_Matched$sys.uniqueId,
                                          Patient_NonHotspot_Matched$sys.uniqueId))))

for (id_num  in 1:length(patient.list)) {
  patient_id <- patient.list[id_num]
  
  # Extract rows of matched trials per patient
  DF_patient_Internal <- unique(OnCore_Biomarker_Matched[which(OnCore_Biomarker_Matched$sys.uniqueId == patient_id),])
  DF_patient_NCI <- unique(Patient_Variant_Matched[which(Patient_Variant_Matched$sys.uniqueId == patient_id),])
  DF_patient_NCI.NonHotspot <- unique(Patient_NonHotspot_Matched[which(Patient_NonHotspot_Matched$sys.uniqueId == patient_id),])
  
  # Output trials from OnCore_Biomarker_Matched
  if (nrow(DF_patient_Internal) > 0 ) {
    
    ## Write output to file (Internal Stanford Trial)
    outdir_sub_internal = paste(outdir, "Internal_", OnCore_Biomarker_Report_timestamp, "_", 
                                gsub("_pathogenicFILTER[_][[:alpha:]]{,3}$", "", filterName), "/", sep="")
    if (!dir.exists(outdir_sub_internal)){dir.create(outdir_sub_internal)} 
    
    txt_filename <- paste(outdir_sub_internal, patient_id, ".txt", sep="")
    sink(file = txt_filename, 
         append = FALSE, 
         split = FALSE)
    
    options(max.print=999999)
    
    # Extract patient INFO
    DF_patient_INFO <- DF_patient_Internal[,c(1:33,39,42,69:71)]
    DF_patient_INFO <- unique(DF_patient_INFO)
    
    # Output patient bio
    cat(paste("Ordering Physician: ",DF_patient_INFO$smpl.hasOrderingPhysician, sep=""), "\n","\n")
    cat(paste("Patient ", patient_id, " may qualify for the following clinical trial(s) based on ", 
              DF_patient_INFO$smpl.csmAssay, " of the specimen biopsied on ", DF_patient_INFO$smpl.dateCollected, 
              ". Mutations were identified in ", DF_patient_INFO$smpl.assayName, 
              " assay (reviewed ", DF_patient_INFO$smpl.reportDateReviewed, ").", sep=""), 
        "\n", "\t", paste("Demographic: ", DF_patient_INFO$current.age, "yo ", DF_patient_INFO$smpl.gender, sep=""), 
        "\n", "\t", paste("Primary tumor site: ", DF_patient_INFO$smpl.primaryTumorSite, sep=""), 
        "\n", "\t", paste("Specimen site: ", DF_patient_INFO$smpl.specimenSite, sep=""), 
        "\n", "\t", paste("Dx: ", DF_patient_INFO$smpl.ppDiagnosticSummary, sep=""), "\n","\n","\n")
    remove(DF_patient_INFO)
    
    # Extract list of matched OnCore.No
    trial.list.internal <- sort(unique(DF_patient_Internal$OnCore.No))
    
    # Output matched INFO per trial 
    for (trial_num_internal in 1:length(trial.list.internal)) {
      
      trial_id <- trial.list.internal[trial_num_internal]
      
      # Extract trial INFO
      trial_INFO <- unique(DF_patient_Internal[DF_patient_Internal$OnCore.No == trial_id, 
                                               which(colnames(DF_patient_Internal) == "OnCore.No"):
                                                 which(colnames(DF_patient_Internal) == "Title")])
      
      # Output trial INFO
      cat(paste("No.", trial_num_internal, ": Clinical Trial #", trial_id, ": ", trial_INFO$Title, sep=""),"\n",
          "----------------------------------------------------------------------", "\n")
      cat(paste("Principal Investigator: ", trial_INFO$PI, " (", trial_INFO$PI.Email, ")", sep=""),"\n")
      cat(paste("Primary Clinical Research Coordinator: ", trial_INFO$Primary.CRC, 
                " (", trial_INFO$Primary.CRC.Email, ")", sep=""),
          "\n","\n")
      
      cat("Biomarker criteria:", "\n", "\t", trial_INFO$Biomarker.Description,"\n")
      cat(paste("Disease group criteria: ", trial_INFO$Disease.Group, sep=""),"\n")
      cat(paste("Disease site criteria: ", trial_INFO$Disease.Sites, sep=""),
          "\n","\n")
      
      DF_patient_trial <- DF_patient_Internal[DF_patient_Internal$OnCore.No == trial_id,]
      DF_patient_trial_ouput <- unique(DF_patient_trial[, c("base.gene","smpl.hgvsProtein","smpl.pathogenicityStatus",
                                                            "Biomarker_GeneName","Biomarker_Condition","Biomarker_Detail","Biomarker_Comment")])
      
      cat(paste("Relevant mutation(s) identified in patient ", patient_id, ":", sep=""),"\n")
      for (entry_num in 1:nrow(DF_patient_trial_ouput)) {
        cat("\t", paste(DF_patient_trial_ouput$base.gene[entry_num], ":", DF_patient_trial_ouput$smpl.hgvsProtein[entry_num], 
                        " (", DF_patient_trial_ouput$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
      }
      cat("\n","\n")
      remove(trial_INFO,trial_id,entry_num,DF_patient_trial,DF_patient_trial_ouput)
    }
  
    remove(trial_num_internal)
    
    cat(paste("NOTE: Trial criteria indicated by 'Disease group', 'Disease site' and 'Comments' need to be manually assessed. OnCore Biomarker Report last updated on ", 
              OnCore_Biomarker_Report_timestamp, ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
    sink()     
  }
  
  if (isTRUE(nrow(DF_patient_NCI) > 0  | nrow(DF_patient_NCI.NonHotspot) > 0)) {
    
    trial_num = 0
    assign("trial_num", trial_num, envir = .GlobalEnv)
    
    ## Write output to file (NCI-MATCH Trial)
    outdir_sub_NCI = paste(outdir, "NCI-MATCH_", Patient_Variant_Report_timestamp,"/", sep="")
    if (!dir.exists(outdir_sub_NCI)){dir.create(outdir_sub_NCI)} 
    
    txt_filename <- paste(outdir_sub_NCI, patient_id, ".txt", sep="")
    sink(file = txt_filename, 
         append = FALSE, 
         split = FALSE)
    
    options(max.print=999999)
    
    # Extract patient INFO
    DF_patient_INFO <- rbind(DF_patient_NCI[,c(1:33,39,42,69:71)],
                             DF_patient_NCI.NonHotspot[,c(1:33,39,42,69:71)])
    DF_patient_INFO <- unique(DF_patient_INFO)
    
    # Output patient bio
    cat(paste("Ordering Physician: ",DF_patient_INFO$smpl.hasOrderingPhysician, sep=""), "\n","\n")
    cat(paste("Patient ", patient_id, " may qualify for the following clinical trial(s) based on ", 
              DF_patient_INFO$smpl.csmAssay, " of the specimen biopsied on ", DF_patient_INFO$smpl.dateCollected, 
              ". Mutations were identified in ", DF_patient_INFO$smpl.assayName, 
              " assay (reviewed ", DF_patient_INFO$smpl.reportDateReviewed, ").", sep=""), 
        "\n", "\t", paste("Demographic: ", DF_patient_INFO$current.age, "yo ", DF_patient_INFO$smpl.gender, sep=""), 
        "\n", "\t", paste("Primary tumor site: ", DF_patient_INFO$smpl.primaryTumorSite, sep=""), 
        "\n", "\t", paste("Specimen site: ", DF_patient_INFO$smpl.specimenSite, sep=""), 
        "\n", "\t", paste("Dx: ", DF_patient_INFO$smpl.ppDiagnosticSummary, sep=""), "\n","\n","\n")
    remove(DF_patient_INFO)
     
    # Output trials from Patient_Variant_Matched
    if (nrow(DF_patient_NCI) > 0) {
      
      # Extract list of matched Arm_Name (Inclusion Variants)
      trial.list.NCI <- sort(unique(DF_patient_NCI$Arm_Name))
      
      # Output matched INFO per trial 
      for (trial_num_NCI in 1:length(trial.list.NCI)) {
        
        # Extract and output trial INFO (Inclusion Variants)
        trial_id <- trial.list.NCI[trial_num_NCI]
        DF_patient_trial <- DF_patient_NCI[DF_patient_NCI$Arm_Name == trial_id,]
        DF_patient_trial_ouput <- unique(DF_patient_trial[, c("base.gene","smpl.hgvsProtein","smpl.pathogenicityStatus",
                                                              "Gene_Name","Protein")])
        
        cat(paste("No.", trial_num_NCI, ": NCI-MATCH Trial Treatment ", trial_id, sep=""),"\n", 
            "----------------------------------------------------------------------", "\n")
        cat("Inclusion criteria: ", "\n")
        for (entry_num in 1:nrow(DF_patient_trial_ouput)) {
          cat("\t", paste(DF_patient_trial_ouput$Gene_Name[entry_num], ":", 
                          DF_patient_trial_ouput$Protein[entry_num], sep=""),"\n")
        }
        
        # Extract trial INFO (Inclusion NonHotspot)
        DF_Inclusion_NonHotspot_Patient <- 
          DF_patient_NCI.NonHotspot[DF_patient_NCI.NonHotspot$Arm_Name == trial_id,]
        
        # Output trial INFO (Inclusion NonHotspot)
        if (nrow(DF_Inclusion_NonHotspot_Patient) > 0) {
          cat("Inclusion NonHotspot criteria: ", "\n")
          for (entry_num in 1:nrow(DF_Inclusion_NonHotspot_Patient)) {
            cat("\t", paste(DF_Inclusion_NonHotspot_Patient$Gene_Name[entry_num], ":Exon ", 
                            DF_Inclusion_NonHotspot_Patient$Exon[entry_num], " (Function: ",
                            DF_Inclusion_NonHotspot_Patient$Function[entry_num], "; Variant Class: ",
                            DF_Inclusion_NonHotspot_Patient$oncominevariantclass[entry_num], ")", sep=""),"\n")
          }
          
          # Adjust patient-specific Inclusion NonHotspot dataframe
          DF_patient_NCI.NonHotspot <- 
            DF_patient_NCI.NonHotspot[!(DF_patient_NCI.NonHotspot$sys.uniqueId == patient_id & 
                                          DF_patient_NCI.NonHotspot$Arm_Name == trial_id),]
        }
        
        cat("\n")
        cat(paste("Relevant mutation(s) identified in patient ", patient_id, ":", sep=""),"\n")
        for (entry_num in 1:nrow(DF_patient_trial_ouput)) {
          cat("\t", paste(DF_patient_trial_ouput$base.gene[entry_num], ":", 
                          DF_patient_trial_ouput$smpl.hgvsProtein[entry_num], 
                          " (", DF_patient_trial_ouput$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
        }
        
        Arm_INFO(trial_id = trial_id)
        cat("\n","\n")
        remove(trial_id,entry_num,DF_patient_trial,DF_Inclusion_NonHotspot_Patient,DF_patient_trial_ouput)
      }
    } else {
      trial_num_NCI = 0
    }
    
    # Output trials from Patient_NonHotspot_Matched
    if (nrow(DF_patient_NCI.NonHotspot) > 0 ) {
      
      # Extract list of matched Arm_Name (Inclusion NonHotspot)
      trial.list.NonHotspot <- sort(unique(DF_patient_NCI.NonHotspot$Arm_Name))
      
      # Output matched INFO per trial 
      for (trial_num in 1:length(trial.list.NonHotspot)) {
        assign("trial_num", trial_num, envir = .GlobalEnv)
        
        # Extract and output trial INFO (Inclusion NonHotspot)
        trial_id <- trial.list.NonHotspot[trial_num]
        DF_patient_trial <- DF_patient_NCI.NonHotspot[DF_patient_NCI.NonHotspot$Arm_Name == trial_id,]
        
        cat(paste("No.", (trial_num_NCI + trial_num), ": NCI-MATCH Trial Treatment ", 
                  unique(DF_patient_trial$Arm_Name), sep=""),"\n",
            "----------------------------------------------------------------------", "\n")
        cat("Inclusion NonHotspot criteria: ", "\n")
        for (entry_num in 1:nrow(DF_patient_NCI.NonHotspot)) {
          cat("\t", paste(DF_patient_NCI.NonHotspot$Gene_Name[entry_num], 
                          ":Exon ", DF_patient_NCI.NonHotspot$Exon[entry_num], 
                          " (Function: ", DF_patient_NCI.NonHotspot$Function[entry_num], 
                          "; Variant Class: ", DF_patient_NCI.NonHotspot$oncominevariantclass[entry_num], 
                          ")", sep=""),"\n")
        }
        
        cat("\n", paste("Relevant mutated transcript(s) identified in patient ", patient_id, ":", sep=""),"\n")
        for (entry_num in 1:nrow(DF_patient_NCI.NonHotspot)) {
          cat("\t", paste(DF_patient_NCI.NonHotspot$base.gene[entry_num], ":", DF_patient_NCI.NonHotspot$smpl.hgvsProtein[entry_num], 
                          " (", DF_patient_NCI.NonHotspot$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
        }
        
        Arm_INFO(trial_id = trial_id)
        cat("\n","\n")
        remove(trial_id,entry_num,DF_patient_trial)
      } 
    }
    
    cat(paste("NOTE: NCI-MATCH trial criteria i.e. NonHotspot Rules, IHC Results, Comments, and Disease Exclusion Codes need to be manually assessed.",
      " Patient Variant Report last updated on ", Patient_Variant_Report_timestamp, ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
    sink() 
  }
}

remove(patient_id,patient.list,trial_num,trial.list.internal,trial.list.NCI,txt_filename,
       DF_patient_Internal,DF_patient_NCI,DF_patient_NCI.NonHotspot,Arm_INFO,STAMP_all_variants_QC,
       id_num,trial_num_NCI,trial.list.NonHotspot,Patient_NonHotspot_Matched,Patient_Variant_Matched,
       DF_Comments,DF_Exclusion_NonHotspot_Rules,DF_IHC_Results,OnCore_Biomarker_Matched,
       Arm_Exclusion,DF_Disease_Exclusion_Codes,DF_exclusion_output,outdir_sub_internal,outdir_sub_NCI)

# Delete intermediate file
#----------------------------------------------
int_file_01 = paste(getwd(), "/ClinicalTrialMatching/Biomarker_Report_LongFormat_", OnCore_Biomarker_Report_timestamp, ".csv", sep="")

if (isTRUE(deleteIntermediateFile)) {
  if (file.exists(int_file_01)){file.remove(int_file_01)}
  if (file.exists(int_file_02)){file.remove(int_file_02)}
  if (file.exists(int_file_03)){file.remove(int_file_03)}
  if (file.exists(int_file_04)){file.remove(int_file_04)}
  if (file.exists(int_file_05)){file.remove(int_file_05)}
}

remove(int_file_01,int_file_02,int_file_03,int_file_04,int_file_05,outdir)
