setwd("~/Documents/ClinicalDataScience_Fellowship/")

## Load relevant active clinical trial files
#----------------------------------------------
OnCore_Biomarker_Report <- 
  read.csv(file = paste("ClinicalTrialMatching/Biomarker_Report_LongFormat_", OnCore_Biomarker_Report_timestamp,".csv", sep=""),
           header = TRUE, na.strings = "NA",stringsAsFactors = FALSE,sep = ",")

PATIENT_VARIANT_REPORT_QC <- 
  import_list(paste("ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_.xlsx", sep=""), setclass = "tbl")
invisible(capture.output(lapply(names(PATIENT_VARIANT_REPORT_QC),
                                function(x) assign(x, PATIENT_VARIANT_REPORT_QC[[x]],
                                                   envir = .GlobalEnv))))
remove(PATIENT_VARIANT_REPORT_QC,DF_Comments,DF_Exclusion_NonHotspot_Rules,
       DF_Exclusion_Variants,DF_IHC_Results,DF_Disease_Exclusion_Codes)

## Load relevant patient mutation files
#----------------------------------------------
STAMP_all_variants_QC <- 
  read.csv(file = paste("ClinicalTrialMatching/", Syapse_Export_timestamp, "_syapse_export_DF_STAMP_VariantAnno.csv", sep=""),
           header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")
ncol_STAMP <- as.numeric(ncol(STAMP_all_variants_QC))

## Output file of matched mutations with clinical trials
#----------------------------------------------
DF_Output_OnCore_Biomarker <- data.frame(matrix(NA, ncol = ncol_STAMP))
colnames(DF_Output_OnCore_Biomarker) <- colnames(STAMP_all_variants_QC)

DF_Output_Patient_Variant <- data.frame(matrix(NA, ncol = (ncol_STAMP +ncol(DF_Inclusion_Variants))))
colnames(DF_Output_Patient_Variant) <- c(colnames(STAMP_all_variants_QC), colnames(DF_Inclusion_Variants))

DF_Output_Patient_NonHotspot <- data.frame(matrix(NA, ncol = (ncol_STAMP +ncol(DF_Inclusion_NonHotspot_Rules))))
colnames(DF_Output_Patient_NonHotspot) <- c(colnames(STAMP_all_variants_QC), colnames(DF_Inclusion_NonHotspot_Rules))

## Extract gene from DFs
#----------------------------------------------
genes.OnCore_Biomarker <- data.frame(Gene=sort(unique(OnCore_Biomarker_Report$Biomarker_GeneName)),
                                     Source="OnCore_Biomarker_Report")
genes.Inclusion_Variants <- data.frame(Gene=sort(unique(DF_Inclusion_Variants$Gene_Name)),
                                       Source="Inclusion_Variants")
genes.Inclusion_NonHotspot <- data.frame(Gene=sort(unique(DF_Inclusion_NonHotspot_Rules$Gene_Name)),
                                         Source="Inclusion_NonHotspot_Rules")

DF_Gene_clinicaltrials <- rbind(genes.OnCore_Biomarker,
                                genes.Inclusion_Variants,genes.Inclusion_NonHotspot)
DF_Gene_clinicaltrials$Gene <- as.character(DF_Gene_clinicaltrials$Gene)
gene.clinicaltrials <- sort(unique(DF_Gene_clinicaltrials$Gene))

#################################################################################################

# # Identify test parameters
# for (patient_num in 1:length(patient.list)) {
#   patient_id <- patient.list[patient_num]
#   Patient_Output <- STAMP_all_variants_QC[which(STAMP_all_variants_QC$sys.uniqueId == patient_id),
#                                           c("smpl.assayName","smpl.reportDateReviewed","base.gene",
#                                             "smpl.specimenSite")]
#   if (nrow(Patient_Output) > 25) {
#     print(patient_id)
#     print(length(unique(Patient_Output$base.gene)))
#   }
# }

## Assume Age.Group for all clinical trials = ADULT
# Child (birth-17) | Adult (18-64) | Older Adult (65+)
# https://clinicaltrials.gov/ct2/about-studies/glossary
#----------------------------------------------
if (isTRUE(adult.group_FILTER)) {
  STAMP_all_variants_QC <- STAMP_all_variants_QC[STAMP_all_variants_QC$current.age >= 18,]
}
patient.list <- unique(STAMP_all_variants_QC$sys.uniqueId)

# Generate output per unique patient
for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  ## Extract STAMP entries for individual patient = DF_patient
  DF_patient <- STAMP_all_variants_QC[which(STAMP_all_variants_QC$sys.uniqueId == patient_id),]
  gene.patient <- sort(unique(DF_patient$base.gene))
  patient_age <- as.numeric(unique(DF_patient$current.age))
  
  # Patient output to screen
  Patient_Output <- STAMP_all_variants_QC[which(STAMP_all_variants_QC$sys.uniqueId == patient_id),
                                          c("base.gene","smpl.hgvsCoding","smpl.hgvsProtein",
                                            "var.type","var.anno","smpl.pathogenicityStatus",
                                            "smpl.primaryTumorSite","primaryTumorSite.category")]
  Patient_Output$OnCore_Report_Status <- NA
  Patient_Output$Patient_Variant_Inclusion_Status <- NA
  Patient_Output$Patient_Variant_NonHotspot_Status <- NA
  
  # Generate output per unique gene
  for (gene_num in 1:length(gene.patient)) {
    gene_id <- gene.patient[gene_num]
    
    ## Extract matches from  OnCore_Biomarker_Report = DF_Gene_OnCore_Biomarker
    #----------------------------------------------
    if (gene_id %in% genes.OnCore_Biomarker$Gene) {
      
      DF_Gene_OnCore_Biomarker <-
        OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$Biomarker_GeneName == gene_id),]
      
      ## Match Biomarker_Condition
      #----------------------------------------------
      biomarker.condition.trial <- unique(DF_Gene_OnCore_Biomarker$Biomarker_Condition)
      biomarker.condition.patient <- unique(DF_patient[DF_patient$base.gene == gene_id,"var.anno"])
      
      # Remove clinical trials with Biomarker_Condition without patient match 
      for (anno_num in 1:length(biomarker.condition.trial)) {
        if (!(biomarker.condition.trial[anno_num] %in% biomarker.condition.patient)) {
          DF_Gene_OnCore_Biomarker <- DF_Gene_OnCore_Biomarker[DF_Gene_OnCore_Biomarker$Biomarker_Condition !=
                                                                 biomarker.condition.trial[anno_num],]
        }
      }
      
      for (anno_num_pt in 1:length(biomarker.condition.patient)) {
        if (biomarker.condition.patient[anno_num_pt] %in% biomarker.condition.trial) {
          
          ## Match primaryTumorSite.category from patient with Disease.Group.category from clinical trial
          #----------------------------------------------
          Disease.Group.category.trial <- unique(DF_Gene_OnCore_Biomarker$Disease.Group.category)
          primaryTumorSite.category.patient <- unique(DF_patient[DF_patient$base.gene == gene_id,"primaryTumorSite.category"])
          disease_group_name <- NA
          disease_group_match <- NA
          
          ## Modify based on whether disease.group_FILTER == TRUE
          if (isTRUE(disease.group_FILTER)) {
            if ("any site" %in% Disease.Group.category.trial) {
              disease_group_name <- "any site"
              disease_group_match <- TRUE
            } else {
              for (site_num in 1:length(Disease.Group.category.trial)) {
                if (isTRUE(grepl(Disease.Group.category.trial[site_num], primaryTumorSite.category.patient))) {
                  disease_group_name <- append(disease_group_name, Disease.Group.category.trial[site_num])
                  disease_group_match <- append(disease_group_match, "TRUE")
                } else {
                  disease_group_match <- append(disease_group_match, "FALSE")
                }
              }
            }
            
            if ("TRUE" %in% disease_group_match) {
              disease_group_match <- TRUE
            } else {
              disease_group_match <- FALSE
            }
            disease_group_name <- disease_group_name[which(!is.na(disease_group_name))]
            
          } else {
            disease_group_name <- Disease.Group.category.trial
          }
          
          if (isTRUE(disease_group_match == TRUE | is.na(disease_group_match))) {
            
            # Corresponding row in patient file
            pt_rowNo <- which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                STAMP_all_variants_QC$base.gene == gene_id &
                                STAMP_all_variants_QC$var.anno == biomarker.condition.patient[anno_num_pt])
            
            # Corresponding OnCore.No in OnCore_Biomarker_Report
            Candidate_ClinicalTrials_No <- unique(DF_Gene_OnCore_Biomarker[
              which(DF_Gene_OnCore_Biomarker$Biomarker_GeneName == gene_id &
                      DF_Gene_OnCore_Biomarker$Disease.Group.category %in% disease_group_name), 
              c("Disease.Group.category","OnCore.No")])
            
            # Extract OnCore.No to corresponding match of patient mutation in Output file
            DF_Output_pre <- STAMP_all_variants_QC[pt_rowNo, c(1:ncol_STAMP)]
            
            for (Core_No in 1:nrow(Candidate_ClinicalTrials_No)) {
              col_start = as.numeric(ncol_STAMP +1 +2*(Core_No -1))
              col_end = as.numeric(ncol_STAMP +2 +2*(Core_No -1))
              
              DF_Output_pre[, c(col_start:col_end)] <- Candidate_ClinicalTrials_No[Core_No,c(1:2)]
            }
            
            # Append to Output file
            DF_Output_OnCore_Biomarker <- rbind.fill(DF_Output_OnCore_Biomarker, DF_Output_pre)
            
            Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
              "Candidate match (phase 1)"
            
          } else {
            Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
              "Disease Group criteria NOT satisfied"
          }
          
        } else {
          Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id &
                                                      Patient_Output$var.anno == biomarker.condition.patient[anno_num_pt])] <-
            "Biomarker Condition criteria NOT satisfied"
        }
      }
      
    } else {
      Patient_Output$OnCore_Report_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Mutation is NOT found"
    }
    
    ## Extract matches from DF_Inclusion_Variants = DF_Gene_Patient_Variant
    #----------------------------------------------
    if (gene_id %in% genes.Inclusion_Variants$Gene) {

      DF_Gene_Patient_Variant <-
        DF_Inclusion_Variants[which(DF_Inclusion_Variants$Gene_Name == gene_id),]

      ## Match var.type from patient with Variant_Type
      #----------------------------------------------
      var.type.trial <- unique(DF_Gene_Patient_Variant$Variant_Type)
      var.type.patient <- unique(DF_patient$var.type[DF_patient$base.gene == gene_id])

      # Remove clinical trials with Variant_Type without patient match
      for (anno_num in 1:length(var.type.trial)) {
        if (!(var.type.trial[anno_num] %in% var.type.patient)) {
          DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Variant_Type !=
                                                               var.type.trial[anno_num],]
        }
      }

      for (var_num in 1:length(var.type.patient)) {
        gene_var_type_id <- var.type.patient[var_num]

        DF_Gene_Patient_Variant_TYPE <-
          DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Variant_Type == gene_var_type_id,]

        if (nrow(DF_Gene_Patient_Variant_TYPE) > 0) {

          # Corresponding row in patient file
          pt_rowNo <- as.numeric(which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                         STAMP_all_variants_QC$base.gene == gene_id &
                                         STAMP_all_variants_QC$var.type == gene_var_type_id))

          # Extract DF_Inclusion_Variants row to corresponding match of patient mutation in Output file
          pt_rowNo_extract <- STAMP_all_variants_QC[pt_rowNo, c(1:ncol_STAMP)]

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

          Patient_Output$Patient_Variant_Inclusion_Status[which(Patient_Output$base.gene == gene_id)] <-
            "Candidate match (phase 1)"

        } else {
          Patient_Output$Patient_Variant_Inclusion_Status[which(Patient_Output$base.gene == gene_id)] <-
            "Variant type criteria NOT satisfied"
        }
      }

    } else {
      Patient_Output$Patient_Variant_Inclusion_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Mutation is NOT found"
    }

    ## Extract matches from DF_Inclusion_Variants = DF_Gene_Patient_NonHotspot
    #----------------------------------------------
    if (gene_id %in% genes.Inclusion_NonHotspot$Gene) {

      DF_Gene_Patient_NonHotspot <-
        DF_Inclusion_NonHotspot_Rules[which(DF_Inclusion_NonHotspot_Rules$Gene_Name == gene_id),]

      # Corresponding row in patient file
      pt_rowNo <- as.numeric(which(STAMP_all_variants_QC$sys.uniqueId == patient_id &
                                     STAMP_all_variants_QC$base.gene == gene_id))

      # Extract DF_Inclusion_NonHotspot_Rules row to corresponding match of patient mutation in Output file
      pt_rowNo_extract <- STAMP_all_variants_QC[pt_rowNo, c(1:ncol_STAMP)]

      for (pt_rowNo_num in 1:nrow(pt_rowNo_extract)) {

        DF_Output_pre <- data.frame(matrix(NA, ncol = (ncol_STAMP +ncol(DF_Gene_Patient_NonHotspot))))
        colnames(DF_Output_pre) <- c(colnames(pt_rowNo_extract), colnames(DF_Gene_Patient_NonHotspot))

        for (var_row_No in 1:nrow(DF_Gene_Patient_NonHotspot)) {
          DF_Output_pre[var_row_No,] <- cbind(pt_rowNo_extract[pt_rowNo_num,],
                                              DF_Gene_Patient_NonHotspot[var_row_No,])
        }

        # Append to Output file
        DF_Output_Patient_NonHotspot <- rbind.fill(DF_Output_Patient_NonHotspot, DF_Output_pre)
      }

      Patient_Output$Patient_Variant_NonHotspot_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Candidate match (phase 1)"

    } else {
      Patient_Output$Patient_Variant_NonHotspot_Status[which(Patient_Output$base.gene == gene_id)] <-
        "Mutation is NOT found"
    }
  }
  
  ## Write output to file
  txt_filename <- paste(outdir, patient_id, "_Match_Output.txt", sep="")
  sink(file = txt_filename, 
       append = FALSE, 
       split = FALSE)
  
  options(max.print=999999)
  
  cat(paste("Patient: ", patient_id, ":", sep=""),"\n")
  cat("\n")
  print(Patient_Output, row.names = FALSE)
  cat("\n")
  
  sink() 
}

remove(DF_Gene_OnCore_Biomarker, DF_Output_pre,DF_patient,Candidate_ClinicalTrials_No,Core_No,gene_id,
       gene_num,gene.patient,Patient_Output,txt_filename,var.type.patient,patient_age,patient_id,patient_num,
       pt_rowNo,ncol_STAMP,gene_var_type_id,pt_rowNo_num,var_num,var_row_No,pt_rowNo_extract,
       DF_Gene_Patient_Variant_TYPE,DF_Gene_Patient_Variant,biomarker.condition.patient,biomarker.condition.trial,
       disease_group_match,disease_group_name,Disease.Group.category.trial,primaryTumorSite.category.patient,
       anno_num,var.type.trial,genes.Inclusion_NonHotspot,genes.Inclusion_Variants,genes.OnCore_Biomarker,
       DF_Gene_Patient_NonHotspot,DF_Gene_clinicaltrials,gene.clinicaltrials,patient.list)

## Remove rows that are all empty
DF_Output_OnCore_Biomarker <- 
  DF_Output_OnCore_Biomarker[rowSums(is.na(DF_Output_OnCore_Biomarker)) != ncol(DF_Output_OnCore_Biomarker),]
DF_Output_Patient_Variant <- 
  DF_Output_Patient_Variant[rowSums(is.na(DF_Output_Patient_Variant)) != ncol(DF_Output_Patient_Variant),]
DF_Output_Patient_NonHotspot <- 
  DF_Output_Patient_NonHotspot[rowSums(is.na(DF_Output_Patient_NonHotspot)) != ncol(DF_Output_Patient_NonHotspot),]

remove(DF_Inclusion_NonHotspot_Rules,DF_Inclusion_Variants,OnCore_Biomarker_Report,STAMP_all_variants_QC,
       col_end,col_start)
if (isTRUE(exists("site_num"))){remove(site_num)}

## Write to local computer
#----------------------------------------------
write.csv(DF_Output_OnCore_Biomarker,
          file = paste("ClinicalTrialMatching/OnCore_Biomarker_Matched_", OnCore_Biomarker_Report_timestamp, "_", filterName_initial, ".csv", sep=""),
          na = "NA", row.names = FALSE)

write.csv(DF_Output_Patient_Variant,
          file = paste("ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_Matched_Variants.csv", sep=""), 
          na = "NA", row.names = FALSE)

write.csv(DF_Output_Patient_NonHotspot,
          file = paste("ClinicalTrialMatching/Patient_Variant_Report_", Patient_Variant_Report_timestamp, "_QC_Matched_NonHotspot.csv", sep=""), 
          na = "NA", row.names = FALSE)
