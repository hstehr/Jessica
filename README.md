# Installation Instructions
#### Run R from command line
```
$ R
```

#### Installation of relevant R packages from CRAN
```
> install.packages("plyr", dependencies = TRUE)
> install.packages("dplyr", dependencies = TRUE)
> install.packages("Biobase", dependencies = TRUE)
> install.packages("eeptools", dependencies = TRUE)
> install.packages("splitstackshape", dependencies = TRUE)
> install.packages("reshape", dependencies = TRUE)
> install.packages("rio", dependencies = TRUE)
> install.packages("stringr", dependencies = TRUE)
> install.packages("openxlsx", dependencies = TRUE)
```

#### Confirm installation of R packages 
```
> library(plyr)
> library(dplyr)
> library(Biobase)
> library(eeptools)
> library(splitstackshape)
> library(reshape)
> library(rio)
> library(stringr)
> library(openxlsx)
```

#### Exit R 
```
> quit()
```

# Clinical Trial Matching PIPELINE 
#### Arguments
###### args[1]: Directory to save output to.
###### args[2]: Location of STAMP entries.
###### args[3]: Patient ID.
###### args[4]: Location of OnCore Report (Stanford Internal Clinical Trials). To turn off matching, set args == FALSE.
###### args[5]: Location of Patient Variant Report (NCI-MATCH Clinical Trials). To turn off matching, set args == FALSE.
###### args[6]: Directory of pipeline scripts.
###### args[7]: Location of stamp_reference_transcripts.
###### args[8]: Location of exons_ensembl.
###### args[9]: Location of disease exclusion key.


#### EXAMPLE - match to both set of clinical trials.
```
#!/bin/bash

# Directories
data_root="/Users/jessicachen/STAMP_v2.4_reports"
script_root="/Users/jessicachen/ClinicalTrialMatching"

# Reference files
stamp_reference=${script_root}/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt
exons_ensembl=${script_root}/nsembl-Gene-Exon-Annotations/exons_ensembl75.txt
histoDx_key=${script_root}/HistologicalDx_CTEP.csv

# Clinical trials 
OnCore=${script_root}/Biomarker_Report_YYYY-MM.csv
NCI=${script_root}/PATIENT_VARIANT_REPORT_TEMPLATE_YYYY-MM-DD.xlsx

# Extract STAMP sequence files and patient ID
ls ${data_root}/reports/*.variant_report.txt | awk -F"/" '{print $NF}' > ${data_root}/trial_names.tsv
lookup=${data_root}/trial_names.tsv

STAMP_file=paste(data.root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
patient_id="LastNameFirstName_PatientID"

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $STAMP_file $patient_id $OnCore $NCI $script_root $stamp_reference $exons_ensembl $histoDx_key
```

#### EXAMPLE - match to Internal clinical trials ONLY.
```
#!/bin/bash

# Directories
data_root="/Users/jessicachen/STAMP_v2.4_reports"
script_root="/Users/jessicachen/ClinicalTrialMatching"

# Reference files
stamp_reference=${script_root}/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt
exons_ensembl=${script_root}/nsembl-Gene-Exon-Annotations/exons_ensembl75.txt
histoDx_key=${script_root}/HistologicalDx_CTEP.csv

# Clinical trials 
OnCore=${script_root}/Biomarker_Report_YYYY-MM.csv
NCI="FALSE"

# Extract STAMP sequence files and patient ID
ls ${data_root}/reports/*.variant_report.txt | awk -F"/" '{print $NF}' > ${data_root}/trial_names.tsv
lookup=${data_root}/trial_names.tsv

STAMP_file=paste(data.root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
patient_id="LastNameFirstName_PatientID"

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $STAMP_file $patient_id $OnCore $NCI $script_root $stamp_reference $exons_ensembl $histoDx_key
```

#### EXAMPLE - match to NCI-MATCH clinical trials ONLY.
```
#!/bin/bash

# Directories
data_root="/Users/jessicachen/STAMP_v2.4_reports"
script_root="/Users/jessicachen/ClinicalTrialMatching"

# Reference files
stamp_reference=${script_root}/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt
exons_ensembl=${script_root}/nsembl-Gene-Exon-Annotations/exons_ensembl75.txt
histoDx_key=${script_root}/HistologicalDx_CTEP.csv

# Clinical trials 
OnCore="FALSE"
NCI=${script_root}/PATIENT_VARIANT_REPORT_TEMPLATE_YYYY-MM-DD.xlsx

# Extract STAMP sequence files and patient ID
ls ${data_root}/reports/*.variant_report.txt | awk -F"/" '{print $NF}' > ${data_root}/trial_names.tsv
lookup=${data_root}/trial_names.tsv

STAMP_file=paste(data.root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
patient_id="LastNameFirstName_PatientID"

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $STAMP_file $patient_id $OnCore $NCI $script_root $stamp_reference $exons_ensembl $histoDx_key
```

# Mutation plots (STAMP entries)
```
source("Mutation_Hotspot/STAMP_lollipop_PIPELINE.R")
```
