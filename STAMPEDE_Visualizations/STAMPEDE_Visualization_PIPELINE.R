## Generate data visualizations of STAMP database
## Filter for entries from "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"

# Load Libraries 
#----------------------------------------------
suppressMessages(library("easypackages"))
suppressMessages(libraries("dplyr","eeptools","gridExtra","reshape","gtable","grid","plotly",
                           "ggplot2", "ggpubr","rio","devtools","jsonlite"))

## Filters
#----------------------------------------------
saveStaticPlots = TRUE
saveDynamicPlots = TRUE
deleteIntermediateFile = TRUE

# Initialization for online plotting (plotly)
#----------------------------------------------
Sys.setenv("plotly_username"="jwrchen")
Sys.setenv("plotly_api_key"="R71FseH6JGAfY8qq6zsa")
options(browser = 'false')

## Set location of directories
#----------------------------------------------
data.root = "~/Documents/ClinicalDataScience_Fellowship/"
pipeline.root ="~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/"
script.dir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/"

## Set location of reference files 
#----------------------------------------------
stamp_reference_transcripts = paste(pipeline.root, 
                                    "Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt", sep="")
exons_ensembl = paste(pipeline.root, "Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt", sep="")
AA_key=paste(pipeline.root,"AminoAcid_Conversion.csv",sep="")

histoDx.key = paste(pipeline.root,"HistologicalDx_CTEP.tsv",sep="")

# Annotated version has gene domain length for lollipop diagrams
STAMPv2_Annotation <- suppressMessages(import_list(
  paste(data.root,"STAMP/2016-08-23_STAMPv2_Annotation.xlsx",sep=""),setclass = "tbl"))

## Specify temp folder and output files
#----------------------------------------------
tempdir = paste(data.root,"STAMPEDE_Visualizations_temp/",sep="")
if (!dir.exists(tempdir)){dir.create(tempdir)} 

## Write output to file
#----------------------------------------------
err.output = paste(script.dir, "Syapse.err.txt",sep="")
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

out.output = paste(script.dir, "Syapse.out.txt",sep="")
sink(file = out.output, append = FALSE, split = FALSE)
options(max.print=999999)

# Classify pathogenicity statuses
#----------------------------------------------
pathogenic = c("Pathogenic", "Likely Pathogenic")
vus = c("Unknown significance", "Unknown")
benign = c("Likely Benign")

## Generate reference tables
#----------------------------------------------
AA_key_table <- read.csv(file = AA_key, sep = ",")

stamp_reference_transcripts <- 
  read.csv(file = stamp_reference_transcripts, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <-
  read.csv(file = exons_ensembl, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl,
                                  by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)

## Generate annotation tables
#----------------------------------------------
# Generate individual DF per worksheet in excel file
invisible(capture.output(lapply(names(STAMPv2_Annotation), 
                                function(x) assign(x, STAMPv2_Annotation[[x]], envir = .GlobalEnv))))

# Re-format protein domain INFO 
colnames(Domains) <- unlist(Domains[2,])
Domains <- data.frame(Domains[c(3:nrow(Domains)),c(2:ncol(Domains))])

# Re-format protein length INFO 
colnames(Genes) <- unlist(Genes[2,])
Genes <- data.frame(Genes[c(3:nrow(Genes)),c(2:ncol(Genes))])

fusion.gene.list.full <- sort(unique(Genes$Name[which(Genes$Fusions == "yes")]))
cnv.gene.list.full <- sort(unique(Genes$Name[which(Genes$CNV == "yes")]))

# Load Disease & Histological dx ontologies
#----------------------------------------------
setwd(pipeline.root)
source("DiseaseGroupCategories.R")
source("HistologicalDx_CTEP_Match.R")

#################################
## QC STAMP datasets = SNV/Indels, Fusions, CNVs
#################################
## STAMP database non-QC PIPELINE
source("SyapseExport_RetrospectiveAnalysis/STAMP_TestOrder_QC.R")

setwd(script.dir)
## STAMP database QC PIPELINE
source("SNVIndel_QC.R")

## STAMP Fusion database QC PIPELINE
source("Fusion_QC.R")

## STAMP CNV database QC PIPELINE
source("CNV_QC.R")

#################################
## Data Visualizations
#################################
# Specify distribution of timeline
#----------------------------------------------
month.number <- c("2014-01","2014-02","2014-03","2014-04","2014-05","2014-06","2014-07","2014-08","2014-09","2014-10","2014-11","2014-12",
                  "2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08","2015-09","2015-10","2015-11","2015-12",
                  "2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08","2016-09","2016-10","2016-11","2016-12",
                  "2017-01","2017-02","2017-03","2017-04","2017-05","2017-06","2017-07","2017-08","2017-09","2017-10","2017-11","2017-12",
                  "2018-01","2018-02","2018-03","2018-04","2018-05","2018-06","2018-07","2018-08","2018-09","2018-10","2018-11","2018-12",
                  "2019-01","2019-02","2019-03","2019-04")

month.alpha <- c("Jan 2014","Feb 2014","Mar 2014","Apr 2014","May 2014","Jun 2014","Jul 2014","Aug 2014","Sep 2014","Oct 2014","Nov 2014","Dec 2014",
                 "Jan 2015","Feb 2015","Mar 2015","Apr 2015","May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
                 "Jan 2016","Feb 2016","Mar 2016","Apr 2016","May 2016","Jun 2016","Jul 2016","Aug 2016","Sep 2016","Oct 2016","Nov 2016","Dec 2016",
                 "Jan 2017","Feb 2017","Mar 2017","Apr 2017","May 2017","Jun 2017","Jul 2017","Aug 2017","Sep 2017","Oct 2017","Nov 2017","Dec 2017",
                 "Jan 2018","Feb 2018","Mar 2018","Apr 2018","May 2018","Jun 2018","Jul 2018","Aug 2018","Sep 2018","Oct 2018","Nov 2018","Dec 2018",
                 "Jan 2019","Feb 2019","Mar 2019","Apr 2019")

# Load functions
#----------------------------------------------
source("STAMPEDE_Functions.R")

## Entire Dataset
#----------------------------------------------
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_ALL/"
cat(paste("Outdirectory: ", outdir, sep=""),"\n","\n")
if (isTRUE(saveStaticPlots)) {if (!dir.exists(outdir)){dir.create(outdir)}}

cohort_id="all"
source("STAMPEDE_EntireDataset.R")

## PrimaryTumorSite: Iterate through list
#----------------------------------------------
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_PerSite/"
cat(paste("Outdirectory: ", outdir, sep=""),"\n","\n")
if (isTRUE(saveStaticPlots)) {if (!dir.exists(outdir)){dir.create(outdir)}}

source("STAMPEDE_PerSite.R")

## Gene: Iterate through list
#----------------------------------------------
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_PerGene/"
cat("\n",paste("Outdirectory: ", outdir, sep=""))

outdir.lollipop = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_Lollipop/"
cat("\n",paste("Outdirectory (lollipop plots): ", outdir.lollipop, sep=""),"\n","\n")

if (isTRUE(saveStaticPlots)) {
  if (!dir.exists(outdir)){dir.create(outdir)} 
  if (!dir.exists(outdir.lollipop)){dir.create(outdir.lollipop)} 
}

source("STAMPEDE_PerGene.R")

remove(pathogenic,vus,benign)
closeAllConnections()

# Delete temporary directory
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)} 

#################################
## Reference TABLES
#################################
outdir = script.dir

# Generate count tables for genes and primary tumor sites
source("Count_SiteGene.R")

# Generate table of iframe codes **MANUAL curation**
# Fusion genes are limited to those in fusion.gene.list.full
source("iframe_codes.R")
