#!/bin/bash

## Directories
data_root="/Users/ClinicalTrialMatching/"
script_root="/Users/ClinicalTrialMatching/"
outdir="/Users/ClinicalTrialMatching/trials"

STAMP="~/Desktop/STAMP/2019-04-30_syapse_export.csv"
OnCore="~/Desktop/Biomarker_Report_2019-03.csv"
NCI="~/Desktop/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-22.xlsx"
NCI_ArmRemove="ARM-F,ARM-G,ARM-S2,ARM-Z1C,ARM-Z1F"

# Transcript exon files
stamp_reference_transcripts="/Users/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt"
exons_ensembl="/Users/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt"
histoDx_key="~/Desktop/HistologicalDx_CTEP.csv"

Rscript ClinicalTrial_Matching_PIPELINE.R $data_root $STAMP $OnCore $NCI $NCI_ArmRemove $script_root $outdir $stamp_reference_transcripts $exons_ensembl $histoDx_key
