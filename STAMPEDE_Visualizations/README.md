## STAMPEDE Documentation: Data Visualizations

### User Access
#### Initialization for online plotting https://plot.ly/r/getting-started/
1.	Create a plotly account - it is required to publish charts online. https://plot.ly/api_signup
2.	Set your authentication credentials in the R session https://plot.ly/settings/api 
```
Sys.setenv("plotly_username"="your_plotly_username")
Sys.setenv("plotly_api_key"="your_api_key")
```
3.	Publish your graphs to Plotly with api_create and use filename to title the file in your Plotly account. Set online plot privacy (i.e.public, private or secret). https://plot.ly/r/getting-started/#online-plot-privacy 
library(plotly)
```
p <- plot_ly(midwest, x = ~percollege, color = ~state, type = "box")
api_create(p, filename = "r-docs-midwest-boxplots",
fileopt= “overwrite”, sharing= “public”)
```
4.	Suppress auto open of the created URL in the browser when publishing graphs via api_create().
```
options(browser = 'false')
```

### Dependencies
#### Local dependencies
Installation of following R packages
```
install.packages("plotly", dependencies = TRUE)
install.packages("easypackages")
install.packages("dplyr")
install.packages("eeptools")
install.packages("gridExtra")
install.packages("reshape")
install.packages("gtable")
install.packages("grid")
install.packages("ggplot2") 
install.packages("ggpubr")
install.packages("rio") 
install.packages("devtools")
install.packages("jsonlite")
```

#### Remote dependencies 
Wifi connection is required.

### Input Files
#### Copy Number Variant (CNV)
*Manually curated by Henning Stehr

| FIELD NAME  | EXAMPLE VALUE |
| ------------- | ------------- |
| sys.uniqueId  | TRF-XXXX  |
| Gene  | NKX2  |
| Variant  | AMP  |

#### Single Number Variation (SNV) / Insertions and Deletions (Indels)
*Exported from Syapse

| FIELD NAME	| EXAMPLE VALUE |
| ------------- | ------------- |
| sys.uniqueId | TRF-XXXX |
| base.dob | MM/DD/YY |
| smpl.gender | Female |
| smpl.histologicalDiagnosis | Squamous cell carcinoma |
| smpl.primaryTumorSite | Lung |
| smpl.specimenType | formalin-fixed paraffin embedded tissue (FFPE) |
| smpl.percentTumor | 40 |
| smpl.assayName | STAMP - Solid Tumor Actionable Mutation Panel (130 genes) |
| smpl.hasOrderingPhysician | LAST_NAME, FIRST_NAME (PHYS-XXX) |
| sys.date_changed | YYYY-MM-DD 00:36:15.863000+00:00 |
| sys.date_created | YYYY-MM-DD 16:40:47.334000+00:00 |
| smpl.reportDateReviewed | YYYY-MM-DD |
| smpl.dateReceived | YYYY-MM-DD |
| smpl.amendedString | NA |
| smpl.amendmentReason | NA |
| sys.label | TP53 c.743G>A (p.Arg248Gln) |
| base.gene | TP53 |
| smpl.pathogenicityStatus | Pathogenic |
| smpl.hgvsCoding | c.743G>A |
| smpl.hgvsProtein | p.Arg248Gln |
| smpl.transcript | NM_000546.5 |
| smpl.genomicDescription | chr17:g.7577538C>T |
| base.chromosome | chr17 |
| smpl.hgvsGenomic | g.7577538C>T |
| smpl.chromosomePositionStart | 7577538 |

#### Fusion
*Exported from Syapse

| FIELD NAME	| EXAMPLE VALUE |
| ------------- | ------------- |
| sys.uniqueId | TRF-XXXX |
| smpl.gender | Male |
| base.dob | YYYY-MM-DD |
| smpl.histologicalDiagnosis | Adenocarcinoma |
| smpl.primaryTumorSite | lung |
| smpl.specimenType | formalin-fixed paraffin embedded tissue (FFPE) |
| smpl.percentTumor | 70 |
| smpl.assayName | STAMP - Solid Tumor Actionable Mutation Panel (130 genes) |
| smpl.hasOrderingPhysician | LAST_NAME, FIRST_NAME (PHYS-XX) |
| sys.date_changed | YYYY-MM-DD 16:08:06.236000+00:00 |
| sys.date_created | YYYY-MM-DD 17:43:04.335000+00:00 |
| smpl.reportDateReviewed | YYYY-MM-DD |
| smpl.dateReceived | YYYY-MM-DD |
| smpl.amendedString | None |
| smpl.amendmentReason | None |
| sys.label | SND1-BRAF |
| smpl.fusionGene1 | SND1 |
| smpl.fusionGene2 | BRAF |

### Source Code
#### Source code location
> https://github.com/SHCMolPathLab/Jessica
> https://stanfordmedicine.box.com/s/m15thsbokgqjvxlciso5wwi393uomaw9 

#### Source code modifications 
Customize parameters for deployment i.e. plotly login information, location of files directories and reference files.
```
STAMPEDE_Visualization_PIPELINE.R
```

#### Modification to datasets
Each R script performs quality control and applies appropriate filters on the corresponding dataset i.e. SNV/Indels, Fusions, CNVs.
```
CNV_QC.R
Fusion_QC.R
SNVIndel_QC.R
TestOrder_QC.R
```

R script generates functions.
```
STAMPEDE_Functions.R
```

Each R script specifies functions to be iterated through per dataset i.e. entire dataset, each gene, each primary tumor site. 
```
STAMPEDE_EntireDataset.R
STAMPEDE_PerGene.R
STAMPEDE_PerSite.R
```

R script extracts total gene and primary tumor site counts for website.
```
Count_SiteGene.R
```

### Extraction of iframe codes
iframe code per plot is found on “Share” component of Plotly file. 

Example
iframe code to embed for plot (https://plot.ly/~jwrchen/3328/#/)
```
<iframe width="900" height="800" frameborder="0" scrolling="no" src="//plot.ly/~jwrchen/3322.embed"></iframe>
```

R script generates iframe codes for STAMPEDE website – each iteration requires manual curation.
```
iframe_codes.R
```

Comma-delimited file containing directory and iframe code of generated plots - this is done via manual curation.
```
mappings.csv
```
