#!/bin/bash

data_root="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching"

outdir=${data_root}
outdir_anno="Version3"
STAMP=${data_root}/../STAMP/2018-10-18_syapse_export_all_variants_patientNameAndMrnRemoved.csv
OnCore=${data_root}/Biomarker_Report_2018-10.csv
NCI=${data_root}/PATIENT_VARIANT_REPORT_TEMPLATE_2018-12-11.xlsx

Rscript ClinicalTrial_Matching_PIPELINE.R $outdir $outdir_anno $STAMP $OnCore $NCI  
