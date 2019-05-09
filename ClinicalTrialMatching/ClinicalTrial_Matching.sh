#!/bin/bash

# Directories
data_root="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/STAMP_v2.4_reports"
script_root="/Users/jessicachen/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts"
outdir=${data_root}/trials

# OnCore trials: to turn off matching, set args == "FALSE"
OnCore=${script_root}/Biomarker_Report_2019-03.csv
# Remove OnCore.No: to indicate none, set args == "NULL"
OnCore_ArmRemove="ECOG-ACRIN-EAY131"

# NCI-MATCH trials: to turn off matching, set args == "FALSE"
NCI=${script_root}/PATIENT_VARIANT_REPORT_TEMPLATE_2019-04-08.xlsx
# Remove NCI-MATCH ARMS: to indicate none, set args == "NULL"
NCI_ArmRemove="ARM-F,ARM-G,ARM-K1,ARM-M,ARM-S2,ARM-Z1C,ARM-Z1F"

# Reference files
stamp_reference=${script_root}/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt
exons_ensembl=${script_root}/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt
histoDx_key=${script_root}/HistologicalDx_CTEP.csv

# Extract STAMP sequence files and patient ID
ls ${data_root}/reports/*.variant_report.txt | awk -F"/" '{print $NF}' > ${data_root}/trial_names.tsv
while read STAMP
do
	patient=$(basename $STAMP .variant_report.txt)
	CNV=${patient}.cnvs
	Fusion=${patient}.fusions.filtered.txt 

	echo $patient

 	Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root \
 	$patient $data_root/reports/$STAMP $data_root/reports/$CNV $data_root/reports/$Fusion \
 	$OnCore $OnCore_ArmRemove $NCI $NCI_ArmRemove \
 	$script_root $outdir \
 	$stamp_reference $exons_ensembl $histoDx_key

done < ${data_root}/trial_names.tsv  

rm ${data_root}/trial_names.tsv