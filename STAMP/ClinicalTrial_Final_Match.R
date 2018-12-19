## Internal INPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched_merged.csv
## NCI INPUT: ClinicalTrialMatching/DF_Output_Patient_Variant_Matched.csv
## Generate notification email for ordering physician
## Output: ClinicalTrialMatching/Patient_Email_Retrospective/Match_ClinicalTrial_patient_id.txt

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load Library
#----------------------------------------------
library("ggplot2")
library("dplyr")

## Load relevant files
#---------------------------------------------- 
OnCore_Biomarker_Matched <-
  read.csv(file = "ClinicalTrialMatching/OnCore_Biomarker_Matched_FINAL.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

Patient_Variant_Matched <-
  read.csv(file = "ClinicalTrialMatching/Patient_Variant_Matched_FINAL.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

OnCore_Biomarker_Report_timestamp <- as.Date("2018-10-01")
Patient_Variant_Report_timestamp <- as.Date("2018-09-08")

## Matching rate: Plot number of candidate trials per patient
#---------------------------------------------- 
DF_OnCore_Biomarker <- data.frame(table(OnCore_Biomarker_Matched$sys.uniqueId))
for (row_No in 1:nrow(DF_OnCore_Biomarker)) {
  DF_OnCore_Biomarker$No.Internal.Trials[row_No] <- 
    length(unique(OnCore_Biomarker_Matched$OnCore.No[OnCore_Biomarker_Matched$sys.uniqueId == 
                                                       DF_OnCore_Biomarker$Var1[row_No]]))
}

DF_Patient_Variant <- data.frame(table(Patient_Variant_Matched$sys.uniqueId))
for (row_No in 1:nrow(DF_Patient_Variant)) {
  DF_Patient_Variant$No.NCI.Trials[row_No] <- 
    length(unique(Patient_Variant_Matched$Arm_Name[Patient_Variant_Matched$sys.uniqueId == 
                                                     DF_Patient_Variant$Var1[row_No]]))
}

DF_OnCore_Biomarker$Var1 <- as.character(DF_OnCore_Biomarker$Var1)
DF_Patient_Variant$Var1 <- as.character(DF_Patient_Variant$Var1)

DF <- full_join(DF_OnCore_Biomarker, DF_Patient_Variant, by = "Var1")
DF$No.Internal.Trials[is.na(DF$No.Internal.Trials)] <- 0
DF$No.NCI.Trials[is.na(DF$No.NCI.Trials)] <- 0
DF$No.Candidate.Trial <- DF$No.Internal.Trials + DF$No.NCI.Trials

DF$No.Internal.Trials[DF$No.Internal.Trials == 0] <- NA
DF$No.NCI.Trials[DF$No.NCI.Trials == 0] <- NA

tiff(filename = "STAMP/Total_MatchDistribution.tiff",
     width = 7, height = 7, units = "in", res = 150)

total_count <- as.numeric(length(which(!is.na(DF$No.Candidate.Trial))))
ggplot(DF, aes(x = DF$No.Candidate.Trial)) +
  geom_histogram(bins = 5, col="gray") +
  labs(title="Patient Match for NCI-MATCH & Stanford Internal Clinical Trials",
       subtitle = paste("N = ", total_count, sep="")) +
  
  scale_x_continuous(name="No. Clinical Trials per Patient", breaks = seq(0, 5, 1)) +
  scale_y_continuous(name="Number of Individuals", breaks = seq(0,500,25)) + 
  
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold",size=12),
        plot.subtitle = element_text(hjust=1, face="bold",size=12),
        legend.text=element_text(size=10),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))

dev.off()

tiff(filename = "STAMP/OnCore_Biomarker_MatchDistribution.tiff",
     width = 7, height = 7, units = "in", res = 150)

total_count <- as.numeric(length(which(!is.na(DF$No.Internal.Trials))))
ggplot(DF, aes(x = DF$No.Internal.Trials)) +
  geom_histogram(bins = 3, col="gray") +
  labs(title="Patient Match for Stanford Internal Clinical Trials",
       subtitle = paste("N = ", total_count, sep="")) +
  
  scale_x_continuous(name="No. Clinical Trials per Patient", breaks = seq(0, 3, 1)) +
  scale_y_continuous(name="Number of Individuals", breaks = seq(0,500,25)) + 
  
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold",size=12),
        plot.subtitle = element_text(hjust=1, face="bold",size=12),
        legend.text=element_text(size=10),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))

dev.off()

tiff(filename = "STAMP/Patient_Variant_MatchDistribution.tiff",
     width = 7, height = 7, units = "in", res = 150)

total_count <- as.numeric(length(which(!is.na(DF$No.NCI.Trials))))
ggplot(DF, aes(x = DF$No.NCI.Trials)) +
  geom_histogram(bins = 3, col="gray") +
  labs(title="Patient Match for NCI-MATCH Trials",
       subtitle = paste("N = ", total_count, sep="")) +
  
  scale_x_continuous(name="No. Clinical Trials per Patient", breaks = seq(0, 3, 1)) +
  scale_y_continuous(name="Number of Individuals", breaks = seq(0,500,25)) + 
  
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold",size=12),
        plot.subtitle = element_text(hjust=1, face="bold",size=12),
        legend.text=element_text(size=10),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))

dev.off()

print(paste("No patients with n=1 internal trials matched: ", length(which(DF$No.Internal.Trials == 1)), sep="")) # 341
print(paste("No patients with n=2 internal trials matched: ", length(which(DF$No.Internal.Trials == 2)), sep="")) # 50
print(paste("No patients with n=3 internal trials matched: ", length(which(DF$No.Internal.Trials == 3)), sep="")) # 2

print(paste("No patients with n=1 NCI trials matched: ", length(which(DF$No.NCI.Trials == 1)), sep="")) # 185
print(paste("No patients with n=2 NCI trials matched: ", length(which(DF$No.NCI.Trials == 2)), sep="")) # 33
print(paste("No patients with n=3 NCI trials matched: ", length(which(DF$No.NCI.Trials == 3)), sep="")) # 2

print(paste("No patients with n=1 TOTAL trials matched: ", length(which(DF$No.Candidate.Trial == 1)), sep="")) # 327
print(paste("No patients with n=2 TOTAL trials matched: ", length(which(DF$No.Candidate.Trial == 2)), sep="")) # 119
print(paste("No patients with n=3 TOTAL trials matched: ", length(which(DF$No.Candidate.Trial == 3)), sep="")) # 39
print(paste("No patients with n=4 TOTAL trials matched: ", length(which(DF$No.Candidate.Trial == 4)), sep="")) # 3
print(paste("No patients with n=5 TOTAL trials matched: ", length(which(DF$No.Candidate.Trial == 5)), sep="")) # 2

remove(DF,DF_OnCore_Biomarker,DF_Patient_Variant,total_count,row_No)

# Generate output per unique patient = 490
#---------------------------------------------- 
patient.list <- sort(unique(append(OnCore_Biomarker_Matched$sys.uniqueId, Patient_Variant_Matched$sys.uniqueId)))

for (id_num  in 1:length(patient.list)) {
  patient_id <- patient.list[id_num]
  
  DF_patient_Internal <- OnCore_Biomarker_Matched[which(OnCore_Biomarker_Matched$sys.uniqueId == patient_id),]
  DF_patient_NCI <- Patient_Variant_Matched[which(Patient_Variant_Matched$sys.uniqueId == patient_id),]
  
  trial.list.internal <- sort(unique(DF_patient_Internal$OnCore.No))
  trial.list.NCI <- sort(unique(DF_patient_NCI$Arm_Name))
  
  ## Write output to file
  txt_filename <- paste("ClinicalTrialMatching/Patient_Email_Retrospective/Match_ClinicalTrial_", patient_id, ".txt", sep="")
  sink(file = txt_filename, 
       append = TRUE, 
       split = FALSE)
  
  options(max.print=999999)
  
  cat(paste("Ordering Physician: ", unique(DF_patient_Internal$smpl.hasOrderingPhysician), sep=""))
  cat("\n","\n")
  cat(paste("Patient ", patient_id, " (", unique(DF_patient_Internal$current.age), "yo ", unique(DF_patient_Internal$smpl.gender),
            " dx with ", unique(DF_patient_Internal$smpl.ppDiagnosticSummary), ") may qualify for following clinical trials based on ", 
            unique(DF_patient_Internal$smpl.csmAssay), " of the ", unique(DF_patient_Internal$smpl.specimenSite), " (biopsied on ",
            unique(DF_patient_Internal$smpl.dateCollected), "). Mutations identified in ", unique(DF_patient_Internal$smpl.assayName), 
            " assay (reviewed ", unique(DF_patient_Internal$smpl.reportDateReviewed), ").", sep=""),"\n")
  cat("\n","\n")

  if (nrow(DF_patient_Internal) > 0 ) {
    for (trial_num in 1:length(trial.list.internal)) {
      
      trial_id <- trial.list.internal[trial_num]
      DF_patient_trial <- DF_patient_Internal[DF_patient_Internal$OnCore.No == trial_id,]
      
      cat(paste("No.", trial_num, ": Clinical Trial #", unique(DF_patient_trial$OnCore.No), ": ", 
                unique(DF_patient_trial$Title), sep=""),"\n")
      cat("----------------------------------------------------------------------", "\n")
      cat(paste("Principal Investigator: ", unique(DF_patient_trial$PI), " (", unique(DF_patient_trial$PI.Email), ")", sep=""),"\n")
      cat(paste("Primary Clinical Research Coordinator: ", unique(DF_patient_trial$Primary.CRC), " (",
                unique(DF_patient_trial$Primary.CRC.Email), ")", sep=""),"\n")
      cat("\n")
      
      cat(paste("Biomarker criteria: ", unique(DF_patient_trial$Biomarker.Description), sep=""),"\n")
      cat(paste("Disease site criteria: ", unique(DF_patient_trial$Disease.Sites), sep=""),"\n")
      cat(paste("Disease group criteria: ", unique(DF_patient_trial$Disease.Group), sep=""),"\n")
      cat("\n")
      
      cat(paste("Relevant mutation(s) identified in patient ", patient_id, ":", sep=""),"\n")
      for (entry_num in 1:nrow(DF_patient_trial)) {
        cat(paste(DF_patient_trial$base.gene[entry_num], ":", DF_patient_trial$smpl.hgvsProtein[entry_num], 
                  " (", DF_patient_trial$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
      }
      cat("\n","\n")
    }
  } else {
    trial_num = 0
  }

  if (nrow(DF_patient_NCI) > 0 ) {
    for (trial_num_NCI in 1:length(trial.list.NCI)) {
      
      trial_id <- trial.list.NCI[trial_num_NCI]
      DF_patient_trial <- DF_patient_NCI[DF_patient_NCI$Arm_Name == trial_id,]
      
      cat(paste("No.", (trial_num_NCI + trial_num), ": NCI-MATCH Trial Treatment ", 
                unique(DF_patient_trial$Arm_Name), sep=""),"\n")
      cat("----------------------------------------------------------------------", "\n")
      ## Need to incorporate contact personnel information????

      cat("Inclusion criteria: ", "\n")
      for (entry_num in 1:nrow(DF_patient_NCI)) {
        cat(paste(DF_patient_NCI$Gene_Name[entry_num], ":", DF_patient_NCI$Protein[entry_num], sep=""),"\t")
      }
      ## How incorporate other criteria???? 
      cat("\n","\n")
      
      cat(paste("Relevant mutation(s) identified in patient ", patient_id, ":", sep=""),"\n")
      for (entry_num in 1:nrow(DF_patient_NCI)) {
        cat(paste(DF_patient_NCI$base.gene[entry_num], ":", DF_patient_NCI$smpl.hgvsProtein[entry_num], 
                  " (", DF_patient_NCI$smpl.pathogenicityStatus[entry_num], ")", sep=""),"\n")
      }
      cat("\n","\n")
    }
  }
  
  cat(paste("NOTE: Any comments for clinical trial(s) above need to be manually assessed. OnCore Biomarker Report last updated on ", 
            OnCore_Biomarker_Report_timestamp, "; Patient Variant Report last updated on ", Patient_Variant_Report_timestamp,
            ". This email was generated on ", Sys.time(), ".", sep=""),"\n")
  sink() 
}

remove(entry_num,id_num,patient_id,patient.list,trial_id, trial_num,
       trial.list.internal,trial.list.NCI,txt_filename,DF_patient_Internal,DF_patient_trial,
       DF_patient_NCI)
