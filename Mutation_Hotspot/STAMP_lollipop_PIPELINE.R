rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/Mutation_Hotspot/")

#################################
## Customize variables 
#################################
Syapse_Export_timestamp <- format(as.Date("2018-10-18"), format= "%Y-%m-%d")
STAMP_annotation_timestamp <- format(as.Date("2016-08-23"), format= "%Y-%m-%d")
deleteIntermediateFile = TRUE
saveStaticPlots = TRUE
saveDynamicPlots = FALSE

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
source(paste(Syapse_Export_timestamp, "_Syapse_Export_QC.R", sep=""))
  ## QC-parameters: smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding
  ## Calculate patient age > merge with primary tumor site data > filter for STAMP entries 
  ## Output: "Syapse_Export_QC.tsv"

source(paste(Syapse_Export_timestamp, "_Syapse_VariantAnnotate.R", sep=""))
  ## Classification (var.type): Synonymous, Upstream, Intronic, SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
  ## Output: "Syapse_Export_DF_STAMP_4Map.tsv"

# Generate lollipop plots for SNV and Frameshift/In-Frame mutations
#----------------------------------------------
# Initialization for online plotting
Sys.setenv("plotly_username"="jwrchen")
Sys.setenv("plotly_api_key"="R71FseH6JGAfY8qq6zsa")
options(browser = 'false')
source("mutation_hotspot.R")

## Remove Syapse-related files
#----------------------------------------------
setwd("/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/Mutation_Hotspot/")
int_file_01 = paste(getwd(), Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep="")
int_file_02 = paste(getwd(), Syapse_Export_timestamp, "_Syapse_Export_DF_STAMP_4Map.tsv", sep="")

if (isTRUE(deleteIntermediateFile)) {
  if (file.exists(int_file_01)){file.remove(int_file_01)}
  if (file.exists(int_file_02)){file.remove(int_file_02)}
}

remove(int_file_01,int_file_02)
