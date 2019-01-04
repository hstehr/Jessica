rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

#################################
## Customize variables 
#################################
Syapse_Export_timestamp <- format(as.Date("2018-10-18"), format= "%Y-%m-%d")
STAMP_annotation_timestamp <- format(as.Date("2016-08-23"), format= "%Y-%m-%d")
deleteIntermediateFile = TRUE
saveStaticPlots = TRUE

# Classify pathogenicity statuses
pathogenic = c("Pathogenic", "Likely Pathogenic")
vus = c("Unknown significance", "Unknown")
benign = c("Likely Benign")

#################################
## PIPELINE
#################################
# Load Libraries 
#----------------------------------------------
suppressMessages(library("easypackages"))
suppressMessages(libraries("dplyr","eeptools","ggplot2", "ggpubr","rio","plotly","devtools","jsonlite"))

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
remove(DF_STAMP_VariantAnno, Age_Calculation_timestamp)

# Generate lollipop plots for SNV and Frameshift/In-Frame mutations
#----------------------------------------------
# Initialization for online plotting
Sys.setenv("plotly_username"="jwrchen")
Sys.setenv("plotly_api_key"="R71FseH6JGAfY8qq6zsa")
options(browser = 'false')
source("Mutation_Hotspot/mutation_hotspot.R")

## Remove Syapse-related files
#----------------------------------------------
setwd("~/Documents/ClinicalDataScience_Fellowship/")
int_file_01 = paste(getwd(), "/ClinicalTrialMatching/", Syapse_Export_timestamp,
                    "_syapse_export_DF_STAMP_VariantAnno.csv", sep="")
int_file_02 = paste(getwd(), "/STAMP/", Syapse_Export_timestamp,
                    "_syapse_export_all_variants_QC.csv", sep="")
int_file_03 = paste(getwd(), "/Mutation_Hotspot/", Syapse_Export_timestamp,
                    "_syapse_export_DF_NAprotein.csv", sep="")
int_file_04 = paste(getwd(), "/Mutation_Hotspot/", Syapse_Export_timestamp,
                    "_syapse_export_DF_STAMP_4Map.csv", sep="")

if (isTRUE(deleteIntermediateFile)) {
  if (file.exists(int_file_01)){file.remove(int_file_01)}
  if (file.exists(int_file_02)){file.remove(int_file_02)}
  if (file.exists(int_file_03)){file.remove(int_file_03)}
  if (file.exists(int_file_04)){file.remove(int_file_04)}
}

remove(int_file_01,int_file_02,int_file_03,int_file_04)
