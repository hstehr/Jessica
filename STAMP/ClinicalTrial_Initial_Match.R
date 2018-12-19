## Clinical trial INPUT: OnCore_Biomarker_Report, Patient_Variant_Report
## Patient INPUT: ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv
## Function: Match STAMP entries to active internal clinical trials based on 
## age.group (ADULT), biomarker.gene, biomarker.condition, disease_site
## Function: Match STAMP entries to active NCI clinical trials based on 
## age.group (ADULT), Gene_Name, Variant_Type
## OnCore_Biomarker_Report OUTPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched.csv
## Patient_Variant_Report OUTPUT: ClinicalTrialMatching/Patient_Variant_Matched.csv

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load Library
#----------------------------------------------
library("plyr")
library("eeptools")

## Load relevant active clinical trial files
#----------------------------------------------
OnCore_Biomarker_Report <- 
  read.csv(file = "ClinicalTrialMatching/Biomarker_Report_LongFormat.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

DF_Inclusion_Variants <- 
  read.csv(file = "ClinicalTrialMatching/Patient_Variant_Report_Inclusion_Variants.csv",
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

## Load relevant patient mutation files
#----------------------------------------------
STAMP_all_variants_QC <- 
  read.csv(file = "ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv",
           header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")
ncol_STAMP <- as.numeric(ncol(STAMP_all_variants_QC))

## Output file of matched mutations with clinical trials
#----------------------------------------------
DF_Output_OnCore_Biomarker <- data.frame(matrix(NA, ncol = ncol_STAMP))
colnames(DF_Output_OnCore_Biomarker) <- colnames(STAMP_all_variants_QC)

DF_Output_Patient_Variant <- data.frame(matrix(NA, ncol = (ncol_STAMP +ncol(DF_Inclusion_Variants))))
colnames(DF_Output_Patient_Variant) <- c(colnames(STAMP_all_variants_QC), colnames(DF_Inclusion_Variants))

## Extract gene from DFs
#----------------------------------------------
gene.OnCore_Biomarker <- data.frame(Gene=sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName)),
                                    Source="OnCore_Biomarker_Report")
genes.Inclusion_Variants <- data.frame(Gene=sort(unique(DF_Inclusion_Variants$Gene_Name)),
                                       Source="Inclusion_Variants")

DF_Gene_clinicaltrials <- rbind(genes.Inclusion_Variants,gene.OnCore_Biomarker)
DF_Gene_clinicaltrials$Gene <- as.character(DF_Gene_clinicaltrials$Gene)
gene.clinicaltrials <- sort(unique(DF_Gene_clinicaltrials$Gene))
remove(genes.Inclusion_Variants,gene.OnCore_Biomarker)

#################################################################################################

patient.list <- unique(STAMP_all_variants_QC$sys.uniqueId)

for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  ## Extract STAMP entries for individual patient = DF_patient
  DF_patient <- STAMP_all_variants_QC[which(STAMP_all_variants_QC$sys.uniqueId == patient_id),]
  gene.patient <- sort(unique(DF_patient$base.gene))
  patient_age <- unique(DF_patient$current.age)
  
  # Temporary patient output
  Patient_Output <- STAMP_all_variants_QC[which(STAMP_all_variants_QC$sys.uniqueId == patient_id),
                                          c("smpl.assayName","smpl.reportDateReviewed","base.gene",
                                            "smpl.specimenSite")]
  Patient_Output$OnCore_Report_Status <- NA
  Patient_Output$Patient_Variant_Report_Status <- NA
  
  for (gene_num in 1:length(gene.patient)) {
    gene_id <- gene.patient[gene_num]
    
    # Clinical trial status for mutation from OnCore_Biomarker_Report and Patient_Variant_Report = TEMPORARY
    # print(DF_Gene_clinicaltrials[which(DF_Gene_clinicaltrials$Gene == gene_id),], row.names = FALSE)
    
    ## Match base.gene from patient with Gene_Name from clinical trial reports
    #----------------------------------------------
    if (gene_id %in% gene.clinicaltrials) {
      
      ## Extract matches from  OnCore_Biomarker_Report = DF_Gene_OnCore_Biomarker
      #----------------------------------------------
      if (gene_id %in% OnCore_Biomarker_Report$Biomarker_GeneName) {

        DF_Gene_OnCore_Biomarker <-
          OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$Biomarker_GeneName == gene_id),]

        ## Exclude non-adult patients -- assume Age.Group for all clinical trials = ADULT
        ## Age group definition: https://clinicaltrials.gov/ct2/about-studies/glossary
        # Child (birth-17) | Adult (18-64) | Older Adult (65+)
        #----------------------------------------------
        if (patient_age >= 18 & patient_age <=64) {

          ## Match smpl.specimenSite from patient with Disease.Site
          ## Most general site is given priority
          #----------------------------------------------
          diseasesite.list <- unique(DF_Gene_OnCore_Biomarker$Disease.Site)
          diseasesite.patient <- unique(DF_patient[DF_patient$base.gene == gene_id,"smpl.specimenSite"])
          ################################################################
          # How deal with missing data about specimen location???????
          ################################################################
          if (is.na(diseasesite.patient)) {
            diseasesite.NA <- TRUE
            diseasesite.patient <- "Any Site"
          } else {
            diseasesite.NA <- FALSE
          }

          if ("Any Site" %in% diseasesite.list) {
            disease_site_name <- "Any Site"
            disease_site_match <- TRUE
          } else {
            for (site_num in 1:length(diseasesite.list)) {
              if (isTRUE(grepl(diseasesite.list[site_num], diseasesite.patient))) {
                disease_site_name <- diseasesite.list[site_num]
                disease_site_match <- TRUE
              } else {
                disease_site_match <- FALSE
              }
            }
          }

          if (isTRUE(disease_site_match)) {
            ############################################################
            # Match by disease group?????
            ############################################################

            # Corresponding row in patient file
            if (isTRUE(diseasesite.NA)) {
              pt_rowNo <- which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                  STAMP_all_variants_QC$base.gene == gene_id &
                                  is.na(STAMP_all_variants_QC$smpl.specimenSite))
            } else {
              pt_rowNo <- which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                  STAMP_all_variants_QC$base.gene == gene_id &
                                  STAMP_all_variants_QC$smpl.specimenSite == diseasesite.patient)
            }

            # Corresponding OnCore.No in OnCore_Biomarker_Report
            Candidate_ClinicalTrials_No <- DF_Gene_OnCore_Biomarker[
              which(DF_Gene_OnCore_Biomarker$Biomarker_GeneName == gene_id &
                      DF_Gene_OnCore_Biomarker$Disease.Site == disease_site_name), "OnCore.No"]

            # Extract OnCore.No to corresponding match of patient mutation in Output file
            DF_Output_pre <- STAMP_all_variants_QC[pt_rowNo, c(1:ncol_STAMP)]

            for (Core_No in 1:length(Candidate_ClinicalTrials_No)) {
              DF_Output_pre[[ncol_STAMP +1]] <- disease_site_name
              DF_Output_pre[[ncol_STAMP +1 +Core_No]] <- Candidate_ClinicalTrials_No[Core_No]
            }

            # Append to Output file
            DF_Output_OnCore_Biomarker <- rbind.fill(DF_Output_OnCore_Biomarker, DF_Output_pre)

          } else {
            Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
              "Disease site DOES NOT satisfy criteria of clinical trials"
          }

        } else {
          Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
            "Patient age DOES NOT satisfy criteria of clinical trials"
        }

      } else {
        Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
          "Mutation is NOT found in OnCore_Biomarker_Report"
      }
      
      ## Extract matches from DF_Inclusion_Variants = DF_Gene_Patient_Variant
      #----------------------------------------------
      if (gene_id %in% DF_Inclusion_Variants$Gene_Name) {
        
        DF_Gene_Patient_Variant <-
          DF_Inclusion_Variants[which(DF_Inclusion_Variants$Gene_Name == gene_id),]
        
        ## Assume Age.Group for all clinical trials = ADULT?????????
        #----------------------------------------------
        if (patient_age >= 18 & patient_age <=64) {
          
          ## Match var.type from patient with Variant_Type
          #----------------------------------------------
          gene_var_type <- DF_patient$var.type[DF_patient$base.gene == gene_id]
          
          # var_num = 1
          for (var_num in 1:length(gene_var_type)) {
            gene_var_type_id <- gene_var_type[var_num]
            
            DF_Gene_Patient_Variant_TYPE <- 
              DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Variant_Type == gene_var_type_id,]
            
            if (nrow(DF_Gene_Patient_Variant_TYPE) > 0) {
              
              # Corresponding row in patient file
              pt_rowNo <- as.numeric(which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                             STAMP_all_variants_QC$base.gene == gene_id &
                                             STAMP_all_variants_QC$var.type == gene_var_type_id))
              
              # Extract DF_Inclusion_Variants row to corresponding match of patient mutation in Output file
              pt_rowNo_extract <- STAMP_all_variants_QC[pt_rowNo, c(1:ncol_STAMP)]
              
              # pt_rowNo_num=1
              for (pt_rowNo_num in 1:nrow(pt_rowNo_extract)) {
                
                DF_Output_pre <- data.frame(matrix(NA, ncol = (ncol_STAMP +ncol(DF_Gene_Patient_Variant_TYPE))))
                colnames(DF_Output_pre) <- c(colnames(pt_rowNo_extract), colnames(DF_Gene_Patient_Variant_TYPE))
                
                for (var_row_No in 1:nrow(DF_Gene_Patient_Variant_TYPE)) {
                  DF_Output_pre[var_row_No,] <- cbind(pt_rowNo_extract[pt_rowNo_num,],
                                                      DF_Gene_Patient_Variant_TYPE[var_row_No,])
                }
                
                # Append to Output file
                DF_Output_Patient_Variant <- rbind.fill(DF_Output_Patient_Variant, DF_Output_pre)
              }
              
            } else {
              Patient_Output$Patient_Variant_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
                "Patient variant type DOES NOT satisfy criteria of clinical trials"
            }
          }
          
        } else {
          Patient_Output$Patient_Variant_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
            "Patient age DOES NOT satisfy criteria of clinical trials"
        }
          
      } else {
          Patient_Output$Patient_Variant_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
            "Mutation is NOT found in Patient_Variant_Report"
        }
        
      } else {
      Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Mutation is NOT found in OnCore_Biomarker_Report"
      Patient_Output$Patient_Variant_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Mutation is NOT found in OnCore_Biomarker_Report"
    }
  }
  
  print(paste("Patient ", patient_id, ":"))
  # print(Patient_Output[,c(3:6)], row.names = FALSE)
  # cat("\n")
}

remove(site_num,DF_Gene_OnCore_Biomarker, DF_Output_pre,DF_patient,
       Candidate_ClinicalTrials_No,Core_No,gene_id,gene_num,
       disease_site_match,disease_site_name,diseasesite.list,
       diseasesite.NA,diseasesite.patient,gene.patient,Patient_Output,
       patient_age,patient_id,patient_num,pt_rowNo,ncol_STAMP,
       gene_var_type,gene_var_type_id,pt_rowNo_num,var_num,var_row_No,
       pt_rowNo_extract,DF_Gene_Patient_Variant_TYPE,DF_Gene_Patient_Variant)

## Remove rows that are all empty
DF_Output_OnCore_Biomarker <- 
  DF_Output_OnCore_Biomarker[rowSums(is.na(DF_Output_OnCore_Biomarker)) != ncol(DF_Output_OnCore_Biomarker),]
DF_Output_Patient_Variant <- 
  DF_Output_Patient_Variant[rowSums(is.na(DF_Output_Patient_Variant)) != ncol(DF_Output_Patient_Variant),]

## Write to local computer
#----------------------------------------------
write.csv(DF_Output_OnCore_Biomarker,
          file = "ClinicalTrialMatching/OnCore_Biomarker_Matched.csv",
          na = "NA",
          row.names = FALSE)

write.csv(DF_Output_Patient_Variant,
          file = "ClinicalTrialMatching/Patient_Variant_Matched.csv",
          na = "NA",
          row.names = FALSE)
