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

R script generates iframe codes for website - this is done via manual curation.
```
iframe_codes.R
```
