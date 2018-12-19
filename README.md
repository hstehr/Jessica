# Patient database
#----------------------------------------------
20181113_syapse_export_all_variants_QC.R
## Clean up 20181018_syapse_export_all_variants_patientNameAndMrnRemoved.csv
## Specific parameters that are quality checked
## smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding
## Output: "20181113_syapse_export_all_variants_QC.csv"

20181113_syapse_export_all_variants_QC_STAMP_VariantAnno.R
## Subset 20181113_syapse_export_all_variants_QC.csv based on mutation type = 9918 total STAMP entries
## Specific categories of mutations
## Synonymous, Upstream, Intronic, SNV, Frameshift, Indels, Insertions, Deletions, Duplications
## Synonymous, upstream & intronic mutations are not mapped to lollipop plots = 305 entries
## Output: "Mutation_Hotspot/20181114_syapse_export_DF_STAMP_4Map.csv" = 9592 entries
## Output: "Mutation_Hotspot/20181114_syapse_export_DF_NAprotein.csv" >> missing smpl.hgvsProtein = 21 entries

## Classification of mutations into "AMPLIFICATION", "DELETION","FUSION","MUTATION"
## Output: "ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv"


# Active clinical trials WITHIN stanford 
#----------------------------------------------
Biomarker_Report_Convert2Long.R
## Clean up Biomarker Report_2018-10_OnCore_Biomarker_Report.csv
## Specific parameters converted to long format 
## Biomarker.Description, Disease.Sites
## Output: "Biomarker_Report_LongFormat.csv"

Patient_Variant_Report_Convert2Long.R
## Clean up PATIENT_VARIANT_REPORT_TEMPLATE.xlsx
## Output: "Patient_Variant_Report_Inclusion_NonHotspot_Rules.csv"
## Output: "Patient_Variant_Report_Exclusion_NonHotspot_Rules.csv"
## Output: "Patient_Variant_Report_Exclusion_Variants.csv"
## Output: "Patient_Variant_Report_Inclusion_Variants.csv"
## Output: "Patient_Variant_Report_IHC_Results.csv"
## Output: "Patient_Variant_Report_Comments.csv"

ClinicalTrial_Initial_Match.R
## Clinical trial INPUT: OnCore_Biomarker_Report, Patient_Variant_Report
## Patient INPUT: ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv
## Function: Match STAMP entries to active internal clinical trials based on 
## age.group (ADULT), biomarker.gene, biomarker.condition, disease_site
## Function: Match STAMP entries to active NCI clinical trials based on 
## age.group (ADULT), Gene_Name, Variant_Type
## OnCore_Biomarker_Report OUTPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched.csv
## Patient_Variant_Report OUTPUT: ClinicalTrialMatching/Patient_Variant_Matched.csv

Biomarker_Report_Extract4Match.R
## Clinical trial INPUT: ClinicalTrialMatching/Biomarker_Report_LongFormat.csv
## Matched patient INPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched.csv
## Extract candidate clinical trials based on gene_id and disease_site of STAMP entries
## Match STAMP entries based on Biomarker.Description of clinical trials
## Output: ClinicalTrialMatching/OnCore_Biomarker_Matched_merged.csv

Patient_Variant_Report_Extract4Match.R
## Clinical trial INPUT: Patient_Variant_Report
## Matched patient INPUT: ClinicalTrialMatching/Patient_Variant_Matched.csv
## Extract candidate NCI clinical trials based on "hgvs.Protein" components of STAMP entries
## Match STAMP entries based on "Protein" componenets of clinical trials
## Output: ClinicalTrialMatching/DF_Output_Patient_Variant_Matched.csv

ClinicalTrial_Final_Match.R
## Internal INPUT: ClinicalTrialMatching/OnCore_Biomarker_Matched_merged.csv
## NCI INPUT: ClinicalTrialMatching/DF_Output_Patient_Variant_Matched.csv
## Generate notification email for ordering physician
## Output: ClinicalTrialMatching/Patient_Email_Retrospective/Match_ClinicalTrial_patient_id.txt
