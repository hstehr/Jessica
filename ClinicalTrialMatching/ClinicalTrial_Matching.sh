#!/bin/bash

# Directories
data_root="/Users/jessicachen/STAMP_v2.4_reports"
script_root="/Users/jessicachen/ClinicalTrialMatching"

# Reference files
stamp_reference=${script_root}/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt
exons_ensembl=${script_root}/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt
histoDx_key=${script_root}/HistologicalDx_CTEP.csv

# Clinical trials 
OnCore=${script_root}/Biomarker_Report_2018-10.csv
NCI=${script_root}/PATIENT_VARIANT_REPORT_TEMPLATE_2019-02-25.xlsx

# Extract STAMP sequence files and patient ID
ls ${data_root}/reports/*.variant_report.txt | awk -F"/" '{print $NF}' > ${data_root}/trial_names.tsv
lookup=${data_root}/trial_names.tsv

for i in $(seq 1 $(wc -l ${lookup}| awk '{print $1}'))
do
	STAMP=$(cat ${lookup} | awk -F',' -v row=$i '(NR==row){print $0}')
	patient=$(cat ${lookup} | awk -F',' -v row=$i '(NR==row){print $0}' | cut -f1 -d ".")

 	Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $data_root/reports/$STAMP $patient $OnCore $NCI $script_root $stamp_reference $exons_ensembl $histoDx_key
done

rm ${data_root}/trial_names.tsv
