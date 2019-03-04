## Install relevant R packages from CRAN
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

## Confirm installation of R packages 
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

## Clinical Trial Matching (NCI-MATCH and Stanford Internal)
```
bash ClinicalTrialMatching/ClinicalTrial_Matching.sh
```

## Mutation plots (STAMP entries)
```
source("Mutation_Hotspot/STAMP_lollipop_PIPELINE.R")
```
