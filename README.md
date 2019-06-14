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
> install.packages("openxlsx", dependencies = TRUE
> install.packages("tidyr", dependencies = TRUE)
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
> library(tidyr)
```

#### Exit R 
```
> quit()
```

# Clinical Trial Matching PIPELINE 
#### Render file executable for user 
```
$ chmod u+x ClinicalTrial_Matching.sh
$ ClinicalTrial_Matching.sh
```

#### Arguments
###### args[1]: Directory location of "REPORTS" folder.
###### args[2]: Patient ID - string format.
###### args[3]: File location of SNV/Indel entries.
###### args[4]: File location of CNV entries.
###### args[5]: File location of Fusion entries.
###### args[6]: File location of OnCore Report (Stanford Internal Clinical Trials). To turn off matching, set args == FALSE.
###### args[7]: Names of OnCore Arms to remove
###### args[8]: File location of Patient Variant Report (NCI-MATCH Clinical Trials). To turn off matching, set args == FALSE.
###### args[9]: Names of NCI-MATCH Arms to remove.
###### args[10]: Directory location of pipeline scripts.
###### args[11]: File location of OUTPUT directory. 
###### args[12]: File location of stamp_reference_transcripts file.
###### args[13]: File location of exons_ensembl file.
###### args[14]: File location of disease exclusion key file.


#### EXAMPLE - match to both set of clinical trials.
```
#!/bin/bash

script_root="/Users/jessicachen/Documents"
data_root="/Users/jessicachen/STAMP_v2.4_reports"
outdir=${data_root}/trials

OnCore_file=paste(script_root,"/Biomarker_Report_YYYY-MM.csv",sep="")
OnCore_ArmRemove="ECOG-ACRIN-EAY131"

NCI_file=paste(script_root,"/PATIENT_VARIANT_REPORT_TEMPLATE_YYYY-MM-DD.xlsx",sep="")
NCI_ArmRemove="ARM-Z1C,ARM-Z1F"

stamp_reference_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx_key=paste(script_root,"/HistologicalDx_CTEP.csv",sep="")

patient_id="LastNameFirstName_PatientID"
STAMP_file=paste(data_root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
CNV_file=paste(data_root,"/reports/LastNameFirstName_PatientID.cnvs",sep="")
Fusion_file=paste(data_root,"/reports/LastNameFirstName_PatientID.fusions.filtered.txt",sep="")

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $patient_id $STAMP_file $CNV_file $Fusion_file $OnCore_file $OnCore_ArmRemove $NCI_file $NCI_ArmRemove $script_root $outdir $stamp_reference_file $exons_ensembl_file $histoDx_key
```

#### EXAMPLE - match to Internal clinical trials ONLY.
```
#!/bin/bash

script_root="/Users/jessicachen/Documents"
data_root="/Users/jessicachen/STAMP_v2.4_reports"
outdir=${data_root}/trials

OnCore_file=paste(script_root,"/Biomarker_Report_YYYY-MM.csv",sep="")
OnCore_ArmRemove="ECOG-ACRIN-EAY131"

NCI_file="FALSE"
NCI_ArmRemove="NULL"

stamp_reference_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx_key=paste(script_root,"/HistologicalDx_CTEP.csv",sep="")

patient_id="LastNameFirstName_PatientID"
STAMP_file=paste(data_root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
CNV_file=paste(data_root,"/reports/LastNameFirstName_PatientID.cnvs",sep="")
Fusion_file=paste(data_root,"/reports/LastNameFirstName_PatientID.fusions.filtered.txt",sep="")

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $patient_id $STAMP_file $CNV_file $Fusion_file $OnCore_file $OnCore_ArmRemove $NCI_file $NCI_ArmRemove $script_root $outdir $stamp_reference_file $exons_ensembl_file $histoDx_key
```

#### EXAMPLE - match to NCI-MATCH clinical trials ONLY.
```
#!/bin/bash

script_root="/Users/jessicachen/Documents"
data_root="/Users/jessicachen/STAMP_v2.4_reports"
outdir=${data_root}/trials

OnCore_file="FALSE"
OnCore_ArmRemove="NULL"

NCI_file=paste(script_root,"/PATIENT_VARIANT_REPORT_TEMPLATE_YYYY-MM-DD.xlsx",sep="")
NCI_ArmRemove="ARM-Z1C,ARM-Z1F"

stamp_reference_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt",sep="")
exons_ensembl_file=paste(script_root,"/Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt",sep="")
histoDx_key=paste(script_root,"/HistologicalDx_CTEP.csv",sep="")

patient_id="LastNameFirstName_PatientID"
STAMP_file=paste(data_root,"/reports/LastNameFirstName_PatientID.variant_report.txt",sep="")
CNV_file=paste(data_root,"/reports/LastNameFirstName_PatientID.cnvs",sep="")
Fusion_file=paste(data_root,"/reports/LastNameFirstName_PatientID.fusions.filtered.txt",sep="")

Rscript ${script_root}/ClinicalTrial_Matching_PIPELINE.R $data_root $patient_id $STAMP_file $CNV_file $Fusion_file $OnCore_file $OnCore_ArmRemove $NCI_file $NCI_ArmRemove $script_root $outdir $stamp_reference_file $exons_ensembl_file $histoDx_key
```

## Visualization Graphics for STAMPEDE 
```
source("STAMPEDE_Visualizations/STAMPEDE_Visualization_PIPELINE.R")
```

## Retrospective analysis of STAMP entries (SNV Indels, CNVs, Fusions)
```
source("ClinicalTrialMatching/SyapseExport_RetrospectiveAnalysis/ClinicalTrial_Matching_ARGS.R")

```
