#!/bin/bash

## Directories
data_root="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/"
script_root="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/"
outdir="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/trials"

STAMP="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv"

# OnCore trials: to turn off matching, set args == "FALSE"
OnCore="~/Documents/ClinicalDataScience_Fellowship/Biomarker_Report/Biomarker_Report_2019-03.csv"
# Remove OnCore.No: to indicate none, set args == "NULL"
OnCore_ArmRemove="ECOG-ACRIN-EAY131"

# NCI-MATCH trials: to turn off matching, set args == "FALSE"
NCI="~/Documents/ClinicalDataScience_Fellowship/PATIENT_VARIANT_REPORT/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-22.xlsx"
# Remove NCI-MATCH ARMS: to indicate none, set args == "NULL"
NCI_ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

# Transcript exon files
stamp_reference_transcripts="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/STAMP/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt"
exons_ensembl="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/STAMP/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt"
histoDx_key="~/Desktop/Version4b_RetrospectiveAnalysis_20190502/HistologicalDx_CTEP.csv"

Rscript ClinicalTrial_Matching_PIPELINE.R $data_root $STAMP $OnCore $OnCore_ArmRemove $NCI $NCI_ArmRemove $script_root $outdir $stamp_reference_transcripts $exons_ensembl $histoDx_key