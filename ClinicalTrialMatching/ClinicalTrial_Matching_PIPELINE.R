rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

#################################
## Customize variables 
#################################
Age_Calculation_timestamp <- format(as.Date("2018-12-31"), format= "%Y-%m-%d") # Default = Sys.Date() 
Syapse_Export_timestamp <- format(as.Date("2018-10-18"), format= "%Y-%m-%d")
OnCore_Biomarker_Report_timestamp <- format(as.Date("2018-10-01"), format= "%Y-%m")
Patient_Variant_Report_timestamp <- format(as.Date("2018-09-06"), format= "%Y-%m-%d")
# Patient_Variant_Report_timestamp <- format(as.Date("2018-12-11"), format= "%Y-%m-%d")

deleteIntermediateFile = TRUE
adult.group_FILTER = TRUE # APPLY TO ALL
pathogenic_accepted <- c("Pathogenic", "Likely Pathogenic")

# Specify individual filters
disease.group_FILTER = TRUE
disease.site_FILTER = TRUE # Dependent on disease.group_FILTER == TRUE
pathogenic_FILTER = TRUE

# # Run every combination of filters
# # Export script from "## PIPELINE" onward as file = "~/Desktop/Untitled.R"
# disease.group_FILTER_all = c("TRUE","FALSE")
# disease.site_FILTER_all = c("TRUE","FALSE")
# pathogenic_FILTER_all = c("TRUE","FALSE")
# 
# for (group_type in 1:length(disease.group_FILTER_all)) {
#   disease.group_FILTER <- as.logical(disease.group_FILTER_all[group_type])
# 
#   if (isTRUE(disease.group_FILTER)) {
#     for (site_type in 1:length(disease.site_FILTER_all)) {
#       disease.site_FILTER <- as.logical(disease.site_FILTER_all[site_type])
#       for (patho_type in 1:length(pathogenic_FILTER_all)) {
#         pathogenic_FILTER <- as.logical(pathogenic_FILTER_all[patho_type])
#         source("~/Desktop/Untitled.R")
#       }
#     }
#   } else {
#     for (patho_type in 1:length(pathogenic_FILTER_all)) {
#       pathogenic_FILTER <- as.logical(pathogenic_FILTER_all[patho_type])
#       source("~/Desktop/Untitled.R")
#     }
#   }
# }
# 
# remove(disease.group_FILTER_all,disease.site_FILTER_all,pathogenic_FILTER_all,group_type,site_type,patho_type)

#################################
## PIPELINE
#################################
# Load Libraries 
#----------------------------------------------
suppressMessages(library("easypackages"))
suppressMessages(libraries("plyr", "tidyverse", "dplyr", "ggplot2", "eeptools", "splitstackshape", 
                           "reshape", "tidyr", "Biobase", "stringr", "rio", "openxlsx"))
# tidyverse_conflicts()     # Conflicts with dplyr

## Filter indication
#----------------------------------------------
filterName_initial <- "diseaseFILTER_"
if (isTRUE(disease.group_FILTER)) {
  filterName_initial <- paste(filterName_initial, "groupON", sep="")
} else {
  filterName_initial <- paste(filterName_initial, "groupOFF", sep="")
}
if (isTRUE(disease.group_FILTER & disease.site_FILTER)) {
  filterName_int <- "siteON"
} else {
  filterName_int <- "siteOFF"
}
if (isTRUE(pathogenic_FILTER)) {
  pathogenic_pre = "pathogenicFILTER_ON"
} else {
  pathogenic_pre = "pathogenicFILTER_OFF"
}
filterName <- paste(filterName_initial,filterName_int, "_", pathogenic_pre, sep="")

# Clean up Patient data from Syapse
#----------------------------------------------
## QC-parameters: smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding
source(paste("STAMP/", Syapse_Export_timestamp, "_syapse_export_all_variants_QC.R", sep=""))
remove(DF_Full)     ## Output: "syapse_export_all_variants_QC.csv"

## Subset STAMP entries > classify mutations > assign current age > classify primaryTumorSite
## Classification #1: Synonymous, Upstream, Intronic, SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Classification #2: "MUTATION","OTHER"
source(paste(Syapse_Export_timestamp, "_syapse_export_all_variants_QC_STAMP_VariantAnno.R", sep=""))
## Output: "Mutation_Hotspot/syapse_export_DF_[STAMP_4Map | NAprotein].csv"
remove(DF_STAMP_VariantAnno)     ## Output: "syapse_export_DF_STAMP_VariantAnno.csv""

# Clean up Clinical trial data
#----------------------------------------------
## QC-parameters: Biomarker.Description, Disease.Sites (most general), Disease.Group
source("ClinicalTrialMatching/Biomarker_Report_Convert2Long.R")
remove(OnCore_Biomarker_Report)     ## Output: "Biomarker_Report_LongFormat.csv"

## Partition criteria into separate files and general QC
source("Patient_Variant_Report_Convert2Long.R")
remove(list_of_datasets)     ## Output: "Patient_Variant_Report_QC.xlsx"

# PRE_Matching
#----------------------------------------------
## Function (Internal): age.group (ADULT), biomarker.gene, biomarker.condition ("MUTATION"), disease_group (diseaseFILTER)
## Function (NCI-MATCH): age.group (ADULT), Gene_Name, Variant_Type (Inclusion ONLY)
outdir = paste(getwd(), "/Retrospective_Analysis/Internal_", OnCore_Biomarker_Report_timestamp, 
               "_NCI-MATCH_", Patient_Variant_Report_timestamp, "_", filterName_initial, "/", sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
source("ClinicalTrial_Initial_Match.R")
remove(DF_Output_OnCore_Biomarker,DF_Output_Patient_Variant,DF_Output_Patient_NonHotspot)
## OUTPUT: OnCore_Biomarker_Matched.csv, Patient_Variant_Matched.csv, Patient_NonHotspot_Matched.csv

# Merge patient entries with candidate clinical trials 
#----------------------------------------------
## Extract candidate clinical trials based on disease_site (FILTER) and Biomarker_Detail
source("ClinicalTrialMatching/Biomarker_Report_Extract4Match.R")
remove(DF_Output_Biomarker_Matched_FINAL)     ## Output: OnCore_Biomarker_Matched_FINAL.csv

## Match NCI-MATCH "Protein" with STAMP "hgvs.Protein"
source("ClinicalTrialMatching/Patient_Variant_Report_Extract4Match.R")
remove(DF_Output_Patient_Variant_Matched,DF_Output_Patient_NonHotspot_Matched)
## Output: Patient_Variant_Report_QC_Matched_FINAL.xlsx

# FINAL_Matching
#----------------------------------------------
## Generate summary text file of candidate clinical trial(s) per patient
outdir = paste(getwd(), "/ClinicalTrialMatching/Retrospective_Analysis/Internal_", OnCore_Biomarker_Report_timestamp, 
               "_NCI-MATCH_", Patient_Variant_Report_timestamp, "_", filterName, "/", sep="")
if (!dir.exists(outdir)){dir.create(outdir)} 
source("ClinicalTrialMatching/ClinicalTrial_Final_Match.R")     ## Output: Match_ClinicalTrial_patient_id.txt

# Remove Syapse-related files
#----------------------------------------------
int_file_01 = paste(getwd(), "/ClinicalTrialMatching/",Syapse_Export_timestamp, "_syapse_export_DF_STAMP_VariantAnno.csv", sep="")
int_file_02 = paste(getwd(), "/STAMP/",Syapse_Export_timestamp, "_syapse_export_all_variants_QC.csv", sep="")

if (file.exists(int_file_01)){file.remove(int_file_01)}
if (file.exists(int_file_02)){file.remove(int_file_02)}

remove(int_file_01,int_file_02)
