## Generate data visualizations of STAMP database
## Filter for entries from "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"
## Incorporate lollipop plots per gene into pipeline
rm(list=ls())

# Load Libraries 
#----------------------------------------------
suppressMessages(library("easypackages"))
suppressMessages(libraries("dplyr","eeptools","gridExtra","reshape","gtable","grid","plotly"))
suppressMessages(libraries("ggplot2", "ggpubr","rio","devtools","jsonlite"))

## Filters
#----------------------------------------------
saveStaticPlots = TRUE
saveDynamicPlots = TRUE
deleteIntermediateFile = TRUE

## Load file locations
#----------------------------------------------
data.root = "~/Documents/ClinicalDataScience_Fellowship/STAMP/"
script.root = "~/Documents/ClinicalDataScience_Fellowship/ClinicalTrialMatching/PIPELINE_scripts/"

STAMP.file = paste(data.root, "2018-10-18_syapse_export_all_variants_patientNameAndMrnRemoved.csv", sep="")
stamp_reference_transcripts = paste(script.root, "Ensembl-Gene-Exon-Annotations/stamp_reference_transcripts.txt", sep="")
exons_ensembl = paste(script.root, "Ensembl-Gene-Exon-Annotations/exons_ensembl75.txt", sep="")
aminoAcid_conversion = paste(data.root, "AminoAcid_Conversion.csv", sep="")

## Specify temp folder
#----------------------------------------------
tempdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations_temp/"
if (!dir.exists(tempdir)){dir.create(tempdir)} 

## Load files 
#----------------------------------------------
STAMP_DF <- 
  read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = ",")

histoDx.key = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-0223_HistologicalDx_CTEP.csv"

AminoAcid_Conversion <- 
  read.csv(file = aminoAcid_conversion, header = TRUE, na.strings = c(""," "), stringsAsFactors = FALSE, sep = ",")

stamp_reference_transcripts <- 
  read.csv(file = stamp_reference_transcripts, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

exons_ensembl <-
  read.csv(file = exons_ensembl, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Merge Gene-Exon Key Table
stamp_reference_full <- left_join(stamp_reference_transcripts,exons_ensembl,
                                  by = c("Transcript" = "enst"))
remove(stamp_reference_transcripts,exons_ensembl)

## Specify script location
setwd(script.root)

## Disease Ontology 
source("DiseaseGroupCategories.R")

## Histological Dx Ontology
source("HistologicalDx_CTEP_Match.R")

## Timestamp
Syapse_Export_timestamp <- 
  format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",
                      sub(".*/", "", STAMP.file))), format= "%Y-%m-%d")

# Specify output file
out.ouput = paste("~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/",
                  Syapse_Export_timestamp,"_Syapse.out.txt",sep="")
err.output = paste("~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/",
                   Syapse_Export_timestamp,"_Syapse.err.txt",sep="")

## Write output to file
#----------------------------------------------
sink(file = err.output, append = FALSE, split = FALSE)
options(max.print=999999)
sink()

sink(file = out.ouput, append = FALSE, split = FALSE)
options(max.print=999999)

## Print parameters to output
#----------------------------------------------
cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n",
    "\t", "Static plots GENERATED: ", saveStaticPlots, "\n",
    "\t", "Dynamic plots GENERATED: ", saveDynamicPlots, "\n")

#################################
## STAMP database QC
#################################

# Clean up patient data from Syapse
#----------------------------------------------
source("../Syapse_Export_QC.R")
source("../Syapse_VariantAnnotate.R")

# Filter entries for complete cases
#----------------------------------------------
STAMP_DF$AssayReportDateReviewed <- as.Date(STAMP_DF$AssayReportDateReviewed, format = "%m/%d/%y")

# Percent tumor
STAMP_DF_rep <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-01-31_syapse_export_all.csv", 
                         header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = ",")
STAMP_DF_rep <- unique(STAMP_DF_rep[,c("UNIQUE_ID","TUMOR_PURITY","DATE_REPORTED")])
STAMP_DF_rep <- STAMP_DF_rep[which(!is.na(STAMP_DF_rep$TUMOR_PURITY)),]
STAMP_DF_rep <- STAMP_DF_rep[which(grepl("comment", STAMP_DF_rep$TUMOR_PURITY, ignore.case = TRUE) == FALSE),]
STAMP_DF_rep$TUMOR_PURITY <- gsub("^less than ", "", STAMP_DF_rep$TUMOR_PURITY)
STAMP_DF_rep <- STAMP_DF_rep[which(!STAMP_DF_rep$TUMOR_PURITY %in% c("None", "N/A", "CMML", "MDS")),]
STAMP_DF_rep <- STAMP_DF_rep[which(!STAMP_DF_rep$TUMOR_PURITY %in% c("10-Jun","10-May","15-Oct","5-Mar","20-Oct","2-Jan")),]
STAMP_DF_rep$TUMOR_PURITY <- gsub("%$", "", STAMP_DF_rep$TUMOR_PURITY)
STAMP_DF_rep$TUMOR_PURITY <- gsub("^>", "", STAMP_DF_rep$TUMOR_PURITY)
STAMP_DF_rep$TUMOR_PURITY <- gsub("^<", "", STAMP_DF_rep$TUMOR_PURITY)
STAMP_DF_rep$TUMOR_PURITY <- gsub("(^[[:digit:]]+)([-][[:digit:]]+$)", "\\1", STAMP_DF_rep$TUMOR_PURITY)
STAMP_DF_rep$TUMOR_PURITY <- ceiling(as.numeric(STAMP_DF_rep$TUMOR_PURITY))

# Manual edit by "DATE_REPORTED" = "AssayReportDateReviewed"
STAMP_DF_rep$TUMOR_PURITY[which(STAMP_DF_rep$UNIQUE_ID == "TRF-1238")] <- as.numeric("15")
STAMP_DF_rep$TUMOR_PURITY[which(STAMP_DF_rep$UNIQUE_ID == "TRF-2378")] <- as.numeric("60") 
STAMP_DF_rep$TUMOR_PURITY[which(STAMP_DF_rep$UNIQUE_ID == "TRF-2775")] <- as.numeric("0") 
# Unclear which is accurate >> remove entry
STAMP_DF_rep <- STAMP_DF_rep[which(STAMP_DF_rep$UNIQUE_ID != "TRF-2506"),]

colnames(STAMP_DF_rep) <- c("PatientID","TUMOR_PURITY","DATE_REPORTED")
STAMP_DF_rep <- unique(STAMP_DF_rep[,c(1:2)])
STAMP_DF <- left_join(STAMP_DF, STAMP_DF_rep, by = "PatientID", all.x=TRUE)

# Specimen type
STAMP_DF_rep <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-10-18_syapse_export_all_variants_patientNameAndMrnRemoved.csv", 
                         header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = ",")
colnames(STAMP_DF_rep) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(STAMP_DF_rep))
STAMP_DF_rep <- unique(STAMP_DF_rep[,c("sys.uniqueId","smpl.specimenType","smpl.reportDateReviewed")])
STAMP_DF_rep <- STAMP_DF_rep[which(!STAMP_DF_rep$smpl.specimenType %in% c("None","other")),]
colnames(STAMP_DF_rep) <- c("PatientID","Specimen_Type","AssayReportDateReviewed")
STAMP_DF_rep <- unique(STAMP_DF_rep[,c(1:2)])
STAMP_DF <- left_join(STAMP_DF, STAMP_DF_rep, by = "PatientID", all.x=TRUE)

# Collapse similar primary tumor site
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver and hepatocellular (liver)"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"

# Complete cases 
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$AssayReportDateReviewed),]
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$PrimaryTumorSite),]
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$Specimen_Type),]
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$TUMOR_PURITY),]

STAMP_DF <- STAMP_DF[which(STAMP_DF$var.type != "Synonymous"),]
STAMP_DF <- STAMP_DF[which(STAMP_DF$PrimaryTumorSite != "unknown"),]
STAMP_DF <- STAMP_DF[which(STAMP_DF$PrimaryTumorSite != "other primary site"),]
STAMP_DF <- STAMP_DF[which(grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*", STAMP_DF$VariantHGVSProtein)),]
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"),]

remove(STAMP_DF_rep)
cat(paste("No. STAMP entries = ", nrow(STAMP_DF), " and No. patients = ", length(unique(STAMP_DF$PatientID)),sep=""),"\n")

write.table(STAMP_DF, file = paste(tempdir, Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#################################
## Data Visualizations
#################################

# Initialization for online plotting
#----------------------------------------------
Sys.setenv("plotly_username"="jwrchen")
Sys.setenv("plotly_api_key"="R71FseH6JGAfY8qq6zsa")
options(browser = 'false')

## FUNCTIONS
#----------------------------------------------
site_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","PrimaryTumorSite")])
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(PrimaryTumorSite) %>% tally())
  colnames(DF_tabulate) <- c("PrimaryTumorSite","No.Patients")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Patients, decreasing = TRUE),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$PrimaryTumorSite <- factor(DF_tabulate$PrimaryTumorSite, 
                                         levels = DF_tabulate$PrimaryTumorSite[order(-DF_tabulate$No.Patients)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Patients)/10) * 10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else {y_increment = 10
    }
  } else {y_increment = 50
  }
  
  # Plot parameters
  height.table = 14
  width.table = 5
  
  height.plot = 15
  width.plot = 25
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Patients) == 1)) {comment = "patient"
  } else {comment = "patients"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=No.Patients)) +
    geom_bar(stat="identity") +
    
    labs(title = "Quantification of Primary Tumor Sites",
         subtitle = paste("N = ", sum(DF_tabulate$No.Patients), " ", comment, sep="")) +
    
    xlab("Primary Tumor Site") + 
    scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("site_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id, "_table.tiff", sep=""),
         width = width.table, height = height.table, units = "in", res = 350)
    grid.arrange(Output.table)
    dev.off()
    
    tiff(filename = paste(outdir, file_id,"_graph.tiff", sep=""),
         width = width.plot, height = height.plot, units = "in", res = 350)
    grid.arrange(plot)
    dev.off()    
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

gene_count_fxn <- function (DF, cohort, outdir) {
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF %>% group_by(VariantGene) %>% tally())
  colnames(DF_tabulate) <- c("Gene","No.Occurrence")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Occurrence, decreasing = TRUE),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$Gene <- factor(DF_tabulate$Gene, 
                             levels = DF_tabulate$Gene[order(-DF_tabulate$No.Occurrence)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Occurrence)/10) * 10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else {y_increment = 10
    }
  } else {y_increment = 50
  }
  
  # Plot parameters
  height.table = 30
  width.table = 5
  
  height.plot = 15
  width.plot = 40
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Occurrence) == 1)) {comment = "occurrence"
  } else {comment = "occurrences"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Gene, y=No.Occurrence)) +
    geom_bar(stat="identity") +
    
    labs(title = "Quantification of Genes",
         subtitle = paste("N = ", sum(DF_tabulate$No.Occurrence), " ", comment, sep="")) +
    
    xlab("Gene Name") + 
    scale_y_continuous(name="Number of Occurrences", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("gene_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    # List of 125 genes does not aesthetically fit into image
    # tiff(filename = paste(outdir, file_id, "_table.tiff", sep=""),
    #      width = width.table, height = height.table, units = "in", res = 350)
    # grid.arrange(Output.table)
    # dev.off()
    
    tiff(filename = paste(outdir, file_id,"_graph.tiff", sep=""),
         width = width.plot, height = height.plot, units = "in", res = 350)
    grid.arrange(plot)
    dev.off()    
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If less than 20 unique genes, output all genes
## If more than 20 unique genes 
# > apply cutoff = round down value of 20th top gene to nearest 5
# > if cutoff is less than n=2, remove genes with n=1
top_gene_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- DF[,c("PatientID","VariantGene")]
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(VariantGene) %>% tally())
  colnames(DF_tabulate) <- c("Gene","Frequency")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Frequency, decreasing = TRUE),]
  
  if (isTRUE(nrow(DF_tabulate) > 20)) {
    cutoff = 5*floor(DF_tabulate$Frequency[[20]]/5) 
    
    if (isTRUE(cutoff < 2)) {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$Frequency >= 2),]
    } else {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$Frequency >= cutoff),]
    }
  }
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 18))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$Gene <- factor(DF_tabulate$Gene, levels = DF_tabulate$Gene[order(-DF_tabulate$Frequency)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$Frequency)/10) * 10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else {y_increment = 10
    }
  } else {y_increment = 50
  }
  
  # Plot parameters
  height = 15
  if (nrow(DF_tabulate) <= 2) {width = 12
  } else if (nrow(DF_tabulate) <= 10) {width = 15
  } else {width = 20
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$Frequency) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Gene, y=Frequency)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Mutated Genes",
         subtitle = paste("N = ", sum(DF_tabulate$Frequency), " / ", nrow(DF), " total STAMP ", comment, sep="")) +
    
    xlab("Gene") +
    scale_y_continuous(name="Frequency of Mutations", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("top_gene_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Top_Gene_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If all variants have n=1, exit function >> no plot generated
## If at least 1 variant have (n > 1) 
## If less than 20 unique variants, output all variants
## If more than 20 unique variants
# > apply cutoff = round down value of 20th top variant to nearest 5
# > if cutoff is less than n=2, remove variants with n=1
top_variant_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- DF[,c("PatientID","VariantGene","VariantHGVSProtein")]
  DF_Fxn$VariantDetail <- paste(DF_Fxn$VariantGene, DF_Fxn$VariantHGVSProtein, sep=" ")
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(VariantDetail) %>% tally())
  colnames(DF_tabulate) <- c("Variant","Frequency")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Frequency, decreasing = TRUE),]
  
  continue_checkpoint <- NA
  if (isTRUE(max(DF_tabulate$Frequency) == 1)) {
    print(paste(cohort, ": all variants have a frequency of n=1", sep=""))
    continue_checkpoint <- as.logical("FALSE")
    
  } else {
    if (isTRUE(nrow(DF_tabulate) > 20)) {
      cutoff = 5*floor(DF_tabulate$Frequency[[20]]/5) 
      
      if (isTRUE(cutoff < 2)) {
        DF_tabulate <- DF_tabulate[which(DF_tabulate$Frequency >= 2),]
      } else {
        DF_tabulate <- DF_tabulate[which(DF_tabulate$Frequency >= cutoff),]
      }
    }
  }
  
  if (isTRUE(is.na(continue_checkpoint))) {
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                                     base_size = 18))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 1, l = 1, r = ncol(Output.table))
    
    # HISTOGRAM
    #----------------------------------------------
    DF_tabulate$Variant <- factor(DF_tabulate$Variant, 
                                  levels = DF_tabulate$Variant[order(-DF_tabulate$Frequency)])
    
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$Frequency)/10) * 10
    if (isTRUE(ymax <= 30)) {
      if (isTRUE(ymax <= 20)) {y_increment = 1
      } else {y_increment = 2
      }
    } else {y_increment = 5
    }
    
    # Plot parameters
    height = 15
    if (nrow(DF_tabulate) <= 10) {
      if (nrow(DF_tabulate) <= 2) {width = 20
      } else {width = 25
      }
    } else {width = 30
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$Frequency) == 1)) {comment = "entry"
    } else {comment = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Variant, y=Frequency)) +
      geom_bar(stat="identity") +
      
      labs(title = "Top Mutated Variants",
           subtitle = paste("N = ", sum(DF_tabulate$Frequency), " / ", nrow(DF), " total STAMP ", comment, sep="")) +
      
      xlab("Variant") +
      scale_y_continuous(name="Frequency of Mutations", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text.y=element_text(size=14),
            axis.text.x=element_text(size=14,angle = 45, hjust = 1),
            axis.title=element_text(size=14,face="bold"))
    
    # Save to local computer
    file_id = paste("top_variant_count_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Top_Variant_Count/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
}

## If less than 20 unique sites, output all sites
## If more than 20 unique sites 
# > apply cutoff = round down value of 20th top site to nearest 5
# > if cutoff is less than n=2, remove sites with n=1
top_site_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- DF[,c("PatientID","PrimaryTumorSite")]
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(PrimaryTumorSite) %>% tally())
  colnames(DF_tabulate) <- c("PrimaryTumorSite","No.Individuals")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Individuals, decreasing = TRUE),]
  
  if (isTRUE(nrow(DF_tabulate) > 20)) {
    cutoff = 5*floor(DF_tabulate$No.Individuals[[20]]/5) 
    
    if (isTRUE(cutoff < 2)) {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Individuals >= 2),]
    } else {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Individuals >= cutoff),]
    }
  }
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 18))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$PrimaryTumorSite <- factor(DF_tabulate$PrimaryTumorSite, 
                                         levels = DF_tabulate$PrimaryTumorSite[order(-DF_tabulate$No.Individuals)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Individuals)/10) * 10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else {y_increment = 50
    }
  } else {y_increment = 100
  }
  
  # Plot parameters
  height = 15
  if (nrow(DF_tabulate) <= 2) {width = 12
  } else if (nrow(DF_tabulate) <= 10) {width = 20
  } else {width = 35
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Individuals) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=No.Individuals)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Primary Tumor Sites Biopsied",
         subtitle = paste("N = ", sum(DF_tabulate$No.Individuals), " / ", nrow(DF), " total STAMP ", comment, sep="")) +
    
    xlab("Primary Tumor Site") +
    scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("top_site_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Top_Site_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

gender_age_distribution_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","PatientGender","PatientAge")])
  
  # Y-axis parameters
  ymax <- ceiling(max(table(DF_Fxn$PatientAge))/10)*10
  if (isTRUE(ymax <= 30)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else {y_increment = 2
    }
  } else {y_increment = 5
  }
  
  # Plot parameters
  height = 7.5
  width = 15
  
  # Subtitle parameters
  if (isTRUE(nrow(DF_Fxn) == 1)) {comment = "patient"
  } else {comment = "patients"
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency = age of patient
  DF_tabulate_full <- data.frame(DF_Fxn %>% group_by(PatientAge,PatientGender) %>% tally())
  colnames(DF_tabulate_full) <- c("PatientAge","Gender","No.Patients")
  row_append <- data.frame(PatientAge=setdiff(seq(1,100), unique(DF_tabulate_full$PatientAge)),
                           Gender=unique(DF_tabulate_full$Gender[[1]]),
                           No.Patients=0, stringsAsFactors = FALSE)
  DF_tabulate_full <- rbind(DF_tabulate_full,row_append)
  DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$PatientAge, decreasing = FALSE),]
  
  gender.missing <- setdiff(c("Female","Male"), unique(DF_tabulate_full$Gender))
  if (length(gender.missing) > 0) {
    DF_tabulate_full <- rbind(DF_tabulate_full,
                              data.frame(PatientAge=50, Gender=gender.missing,No.Patients=0, stringsAsFactors = FALSE))
  }
  
  # Specify Age.Cohort
  DF_Fxn$Age.Cohort <- NA
  for (row_No in 1:nrow(DF_Fxn)) {
    if (isTRUE(DF_Fxn$PatientAge[row_No] < 18)) {DF_Fxn$Age.Cohort[row_No] <- "Child (< 18yo)"
    } else if (isTRUE(DF_Fxn$PatientAge[row_No] < 65)) {DF_Fxn$Age.Cohort[row_No] <- "Adult (18-64yo)"
    } else {DF_Fxn$Age.Cohort[row_No] <- "Older Adult (>= 65yo)"
    }
  }
  
  # Tabulate frequency = age cohorts
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Age.Cohort,PatientGender) %>% tally())
  colnames(DF_tabulate) <- c("Age.Cohort","Gender","No.Patients")
  # Convert table to wide format
  DF_tabulate <- data.frame(cast(DF_tabulate, Age.Cohort ~ Gender), stringsAsFactors = FALSE)
  
  if (isTRUE(!("Male" %in% colnames(DF_tabulate)))) {DF_tabulate$Male <- as.numeric("0")}
  if (isTRUE(!("Female" %in% colnames(DF_tabulate)))) {DF_tabulate$Female <- as.numeric("0")} 
  DF_tabulate$Total <- DF_tabulate$Female + DF_tabulate$Male
  
  DF_tabulate <- DF_tabulate[,c("Age.Cohort","Female","Male","Total")]
  DF_tabulate$Female[which(is.na(DF_tabulate$Female))] <- 0
  DF_tabulate$Male[which(is.na(DF_tabulate$Male))] <- 0
  DF_tabulate$Total[which(is.na(DF_tabulate$Total))] <- 0
  
  cohort.missing <- setdiff(c("Adult (18-64yo)","Child (< 18yo)","Older Adult (>= 65yo)"),
                            unique(DF_tabulate$Age.Cohort))
  if (length(cohort.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Age.Cohort=cohort.missing,
                                    Female=0, Male=0, Total=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate$Age.Cohort <- factor(DF_tabulate$Age.Cohort,
                                   levels = c("Child (< 18yo)","Adult (18-64yo)","Older Adult (>= 65yo)"))
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Age.Cohort),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  plot <- ggplot(DF_tabulate_full, aes(x=PatientAge, y=No.Patients, fill=Gender)) +
    geom_bar(stat="identity") +
    
    labs(title = "Age and Gender Distribution",
         subtitle = paste("N = ", nrow(DF_Fxn), " ", comment, sep="")) +
    
    scale_x_continuous(name="Age", breaks = seq(0,100,5), limits=c(0, 100)) +
    scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    scale_fill_discrete(name = "Gender") +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.955, 0.88),
          legend.text=element_text(size=10),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  # Save to local computer
  file_id = paste("gender_age_distribution_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    DF_tabulate_full_02 <- data.frame(DF_Fxn %>% group_by(PatientAge) %>% tally())
    
    # Customize hover text of domains 
    for (i in 1:2) {
      for (age_No in 1:length(plot_dynamic_int$x$data[[i]]$text)) {
        age_elem <- gsub("(^PatientAge:[[:blank:]]+)([[:digit:]]+)(.*)", "\\2", plot_dynamic_int$x$data[[i]]$text[[age_No]])
        text_add <- paste("<br />Total.No.Patients: ", DF_tabulate_full_02$n[DF_tabulate_full_02$PatientAge == age_elem], sep="")
        
        plot_dynamic_int$x$data[[i]]$text[[age_No]] <- paste(plot_dynamic_int$x$data[[i]]$text[[age_No]], text_add, sep="")
      }  
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Gender_Age_Distribution/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

pt_mutation_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","VariantPathogenicityStatus")])
  
  # TABLE - all variants
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
  DF_tabulate$No.Patients <- as.numeric("0")
  for (row_No in 1:nrow(DF_tabulate)) {
    DF_tabulate$No.Patients[row_No] <- 
      length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
  }
  DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 12))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM - all variants
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Patients)/10)*10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else {y_increment = 10
    }
  } else {y_increment = 50
  }
  
  # Plot parameters
  height = 10
  if (nrow(DF_tabulate) == 1) {width = 5.5
  } else {width = 5 + nrow(DF_tabulate) -1
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Patients) == 1)) {comment1 = "patient"
  } else {comment1 = "patients"
  }
  
  if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Patients[[1]] == 1)) {comment2 = "entry"
  } else {comment2 = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Patients)) +
    geom_bar(stat="identity") +
    
    labs(title = "Mutation Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                          nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
    
    xlab("Mutation Count") + 
    scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("pt_mutation_count_allvariants_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Patient_Mutation_Count/All_Variants/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
  
  # TABLE - pathogenic
  #----------------------------------------------
  DF_Fxn <- unique(DF[which(DF$VariantPathogenicityStatus %in% c("Likely Pathogenic","Pathogenic")),
                      c("PatientID","VariantPathogenicityStatus")])
  
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  
  if (nrow(DF_tabulate_pre) > 0) {
    DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
    DF_tabulate$No.Patients <- as.numeric("0")
    for (row_No in 1:nrow(DF_tabulate)) {
      DF_tabulate$No.Patients[row_No] <- 
        length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
    }
    DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
    
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                                     base_size = 12))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 1, l = 1, r = ncol(Output.table))
    
    # HISTOGRAM - pathogenic
    #----------------------------------------------
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Patients)/10)*10
    if (isTRUE(ymax < 200)) {
      if (isTRUE(ymax <= 20)) {y_increment = 1
      } else if (isTRUE(ymax <= 30)) {y_increment = 2
      } else if (isTRUE(ymax <= 100)) {y_increment = 5
      } else {y_increment = 10
      }
    } else {y_increment = 50
    }
    
    # Plot parameters
    height = 10
    if (nrow(DF_tabulate) == 1) {width = 5.5
    } else {width = 5 + nrow(DF_tabulate) -1
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$No.Patients) == 1)) {comment1 = "patient"
    } else {comment1 = "patients"
    }
    
    if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Patients[[1]] == 1)) {comment2 = "entry"
    } else {comment2 = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Patients)) +
      geom_bar(stat="identity") +
      
      labs(title = "Mutation Distribution (Pathogenic)",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                            nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
      
      xlab("Mutation Count") + 
      scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))
    
    # Save to local computer
    file_id = paste("pt_mutation_count_pathovariants_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Patient_Mutation_Count/Pathogenic_Variants/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
  
  # TABLE - VUS
  #----------------------------------------------
  DF_Fxn <- unique(DF[which(DF$VariantPathogenicityStatus %in% c("Unknown significance","Unknown")),
                      c("PatientID","VariantPathogenicityStatus")])
  
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  
  if (nrow(DF_tabulate_pre) > 0) {
    DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
    DF_tabulate$No.Patients <- as.numeric("0")
    for (row_No in 1:nrow(DF_tabulate)) {
      DF_tabulate$No.Patients[row_No] <- 
        length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
    }
    DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
    
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                                     base_size = 12))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 1, l = 1, r = ncol(Output.table))
    
    # HISTOGRAM - pathogenic
    #----------------------------------------------
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Patients)/10)*10
    if (isTRUE(ymax < 200)) {
      if (isTRUE(ymax <= 20)) {y_increment = 1
      } else if (isTRUE(ymax <= 30)) {y_increment = 2
      } else if (isTRUE(ymax <= 100)) {y_increment = 5
      } else {y_increment = 10
      }
    } else {y_increment = 50
    }
    
    # Plot parameters
    height = 10
    if (nrow(DF_tabulate) == 1) {width = 5.5
    } else {width = 5 + nrow(DF_tabulate) -1
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$No.Patients) == 1)) {comment1 = "patient"
    } else {comment1 = "patients"
    }
    
    if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Patients[[1]] == 1)) {comment2 = "entry"
    } else {comment2 = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Patients)) +
      geom_bar(stat="identity") +
      
      labs(title = "Mutation Distribution (VUS)",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                            nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
      
      xlab("Mutation Count") + 
      scale_y_continuous(name="Number of Individuals", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))
    
    # Save to local computer
    file_id = paste("pt_mutation_count_VUS_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Patient_Mutation_Count/UnknownSignificance_Variants/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
}

specimen_type_stacked_fxn <- function(DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","Specimen_Type")])
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Specimen_Type) %>% tally())
  colnames(DF_tabulate) <- c("Specimen_Type","Frequency")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Frequency, decreasing = TRUE),]
  
  specimen.missing <- setdiff(c("formalin-fixed paraffin embedded tissue (FFPE)",
                                "bone marrow (BM)"),
                              unique(DF_tabulate$Specimen_Type))
  if (length(specimen.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Specimen_Type=specimen.missing,
                                    Frequency=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate$Specimen_Type <- factor(DF_tabulate$Specimen_Type,
                                      levels = c("formalin-fixed paraffin embedded tissue (FFPE)",
                                                 "bone marrow (BM)"))
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Specimen_Type),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.025)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 11))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$Assay <- "STAMP_v2"
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$Frequency)/10) * 10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else {y_increment = 50
    }
  } else {y_increment = 100
  }
  
  # Plot parameters
  height = 10
  width = 5
  
  plot <- ggplot(DF_tabulate, aes(x=Assay, y=Frequency, fill=Specimen_Type)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    
    labs(title = "Specimen Type Distribution",
         subtitle = paste("N = ", length(unique(DF$PatientID)), " specimen", sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Number of Specimen", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="bottom",
          legend.title=element_blank()) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  # Save to local computer
  file_id = paste("specimen_type_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,"_stacked.tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.5), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    # Customize hover text of domains 
    for (i in 1:length(plot_dynamic_int$x$data)) {
      for (elem_No in 1:length(plot_dynamic_int$x$data[[i]]$text)) {
        plot_dynamic_int$x$data[[i]]$text[[elem_No]] <- 
          gsub("^Assay: STAMP_v2<br />","", plot_dynamic_int$x$data[[i]]$text[[elem_No]])
        
        text_add <- paste("<br />Total.No.Samples: ", sum(DF_tabulate$Frequency), sep="")
        
        plot_dynamic_int$x$data[[i]]$text[[elem_No]] <- paste(plot_dynamic_int$x$data[[i]]$text[[elem_No]], text_add, sep="")
      }  
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Specimen_Type_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

tumor_purity_count_fxn <- function (DF, cohort, outdir, width, height) {
  DF_Fxn <- unique(DF[,c("PatientID","TUMOR_PURITY")])
  
  # Round up tumor purity values 
  DF_Fxn$Tumor.Purity = 5*ceiling(DF_Fxn$TUMOR_PURITY/5) 
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Tumor.Purity) %>% tally())
  colnames(DF_tabulate) <- c("Tumor.Purity","No.Patients")
  
  percent.missing <- setdiff(seq(0,100,5), unique(DF_tabulate$Tumor.Purity))
  
  if (length(percent.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Tumor.Purity=percent.missing,
                                    No.Patients=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Tumor.Purity, decreasing = FALSE),]
  
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Patients)/10) * 10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else {y_increment = 50
    }
  } else {y_increment = 100
  }
  
  # Plot parameters
  height = 10
  width = 15
  
  plot <- ggplot(DF_tabulate, aes(x=Tumor.Purity, y=No.Patients)) +
    geom_bar(stat="identity") +
    
    labs(title = "Tumor Purity Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " specimen", sep="")) +
    
    scale_x_continuous(name="Tumor Purity (rounded up)", breaks = seq(0,100,5)) + 
    scale_y_continuous(name="Number of Specimen", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("tumor_purity_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    print(plot)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Tumor_Purity_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

test_volume_timeline_fxn <- function (DF, cohort, outdir, width, height, PerSite) {
  DF_Fxn <- unique(DF[,c("PatientID","AssayReportDateReviewed")])
  # Convert to month/year
  DF_Fxn$AssayReportDateReviewed <- gsub("-[[:digit:]]{2}$","", DF_Fxn$AssayReportDateReviewed)
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(AssayReportDateReviewed) %>% tally())
  colnames(DF_tabulate) <- c("AssayReportDateReviewed","No.Patients")
  
  months.missing <- setdiff(c("2016-08","2016-09","2016-10","2016-11","2016-12",
                              "2017-01","2017-02","2017-03","2017-04","2017-05","2017-06",
                              "2017-07","2017-08","2017-09","2017-10","2017-11","2017-12",
                              "2018-01","2018-02","2018-03","2018-04","2018-05","2018-06",
                              "2018-07","2018-08","2018-09","2018-10"), 
                            unique(DF_tabulate$AssayReportDateReviewed))
  if (length(months.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(AssayReportDateReviewed=months.missing,No.Patients=0, stringsAsFactors = FALSE))
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$AssayReportDateReviewed, decreasing = FALSE),]
  
  for (row_No in 1:nrow(DF_tabulate)) {
    if (isTRUE(row_No == 1)) {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Patients[row_No]
    } else {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Patients[row_No] + DF_tabulate$CumulativeCount[row_No -1]
    }
  }
  
  Output.table <- tableGrob(DF_tabulate[,c("AssayReportDateReviewed","No.Patients")], rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 12))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$CumulativeCount)/10)*10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {
      y_increment = 1
    } else if (isTRUE(ymax <= 30)) {
      y_increment = 2
    } else if (isTRUE(ymax <= 100)) {
      y_increment = 5
    } else {
      y_increment = 50
    }
  } else {y_increment = 100
  }
  
  # Plot parameters
  height = 10
  width = 20
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Patients) == 1)) {comment = "order"
  } else {comment = "orders"
  }
  
  # Reformat "AssayReportDateReviewed" 
  Month_Key <- data.frame(month=c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                                  "Aug","Sep","Oct","Nov","Dec"),
                          digit=seq(1,12,1), 
                          stringsAsFactors = FALSE)
  for (row_No in 1:nrow(DF_tabulate)) {
    row_Month <- Month_Key$month[which(Month_Key$digit == as.numeric(gsub("(^[[:digit:]]{4}[-])([[:digit:]]{2})","\\2",DF_tabulate$AssayReportDateReviewed[row_No])))]
    row_Year <- gsub("(^[[:digit:]]{4})(.*)","\\1",DF_tabulate$AssayReportDateReviewed[row_No])
    DF_tabulate$DateReviewed[row_No] <- paste(row_Month,row_Year,sep=" ")
  }
  DF_tabulate$DateReviewed <- factor(DF_tabulate$DateReviewed,
                                     levels = c("Aug 2016","Sep 2016","Oct 2016","Nov 2016","Dec 2016",
                                                "Jan 2017","Feb 2017","Mar 2017","Apr 2017","May 2017","Jun 2017","Jul 2017","Aug 2017","Sep 2017","Oct 2017","Nov 2017","Dec 2017",
                                                "Jan 2018","Feb 2018","Mar 2018","Apr 2018","May 2018","Jun 2018","Jul 2018","Aug 2018","Sep 2018","Oct 2018"))
  
  plot <- ggplot(DF_tabulate[,c(1:3)], aes(x=AssayReportDateReviewed, y=No.Patients)) +
    geom_bar(stat="identity") +
    
    geom_point(aes(y=CumulativeCount, color="red")) +
    geom_line(aes(y=CumulativeCount, group=1), color="red") +
    
    labs(title = "Order Volume Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment, sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Frequency of Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="none")
  
  plot_dynamic_bar <- ggplot(DF_tabulate[,c(4,2,3)], aes(x=DateReviewed, y=No.Patients)) +
    geom_bar(stat="identity") +
    
    geom_point(aes(y=CumulativeCount, color="red")) +
    geom_line(aes(y=CumulativeCount, group=1), color="red") +
    
    labs(title = "Order Volume Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " orders", sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Frequency of Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="none")
  
  # Save to local computer
  file_id = paste("test_volume_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot_dynamic_bar)
    
    # Customize hover text of domains 
    plot_dynamic_int$x$data[[2]]$text <- 
      gsub("<br />No.Patients:[[:blank:]]+.*","",gsub("<br />colour: red","",plot_dynamic_int$x$data[[2]]$text))
    
    plot_dynamic_int$x$data[[3]]$text <- 
      gsub("<br />No.Patients:[[:blank:]]+.*","",plot_dynamic_int$x$data[[3]]$text)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    
    if (isTRUE(PerSite)) {
      filename_Full = paste("STAMPEDE/PerTumor_OrderVolume/", file_id, sep="")
    } else {
      filename_Full = paste("STAMPEDE/", file_id, sep="")  
    }
    
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

variant_type_distribution_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","var.type")])
  
  DF_Fxn$VariantType <- NA
  for (row_No in 1:nrow(DF_Fxn)) {
    if (isTRUE(DF_Fxn$var.type[row_No] == "SNV")) {
      DF_Fxn$VariantType[row_No] <- "SNV"
    } else if (isTRUE(DF_Fxn$var.type[row_No] %in% 
                      c("Frameshift_Deletion","Frameshift_Delins","Frameshift_Duplication","Frameshift_Insertion"))) {
      DF_Fxn$VariantType[row_No] <- "Frameshift Indel"
    } else if (isTRUE(DF_Fxn$var.type[row_No] %in% 
                      c("Deletion","Delins","Duplication","Insertion"))) {
      DF_Fxn$VariantType[row_No] <- "In-Frame Indel"
    }
  }
  DF_Fxn <- DF_Fxn[,c("PatientID","VariantType")]
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency = age of patient
  DF_tabulate_full <- data.frame(DF_Fxn %>% group_by(VariantType) %>% tally())
  colnames(DF_tabulate_full) <- c("VariantType","No.Entries")
  
  vartype.missing <- setdiff(c("SNV","In-Frame Indel","Frameshift Indel"),
                             unique(DF_tabulate_full$VariantType))
  if (length(vartype.missing) > 0) {
    DF_tabulate_full <- rbind(DF_tabulate_full,
                         data.frame(VariantType=vartype.missing,No.Entries=0, stringsAsFactors = FALSE))
  }
  DF_tabulate_full$VariantType <- factor(DF_tabulate_full$VariantType,
                                         levels = c("SNV","In-Frame Indel","Frameshift Indel"))
  DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$VariantType, decreasing = FALSE),]
  
  Output.table <- tableGrob(DF_tabulate_full, rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate_full$No.Entries)/10)*10
  if (isTRUE(ymax < 200)) {
    if (isTRUE(ymax <= 20)) {
      y_increment = 1
    } else if (isTRUE(ymax <= 30)) {
      y_increment = 2
    } else if (isTRUE(ymax <= 100)) {
      y_increment = 5
    } else {
      y_increment = 25
    }
  } else {y_increment = 50
  }
  
  # Plot parameters
  height = 7.5
  width = 10
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate_full$No.Entries) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate_full, aes(x=VariantType, y=No.Entries)) +
    geom_bar(stat="identity") +
    
    labs(title = "Variant Type Distribution",
         subtitle = paste("N = ", sum(DF_tabulate_full$No.Entries), " ", comment, sep="")) +
    
    xlab("Variant Types") +
    scale_y_continuous(name="Number of Entries", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.955, 0.88),
          legend.text=element_text(size=10),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  # Save to local computer
  file_id = paste("var_type_distribution_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
    dev.off()
  }
  
# Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Var_Type_Distribution/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

#################################
## PrimaryTumorSite: Iterate through list
#################################
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_PerSite/"
if (!dir.exists(outdir)){dir.create(outdir)} 
cat(paste("Outdirectory: ", outdir, sep=""),"\n")

# Across all primary tumor sites
#----------------------------------------------
# Plots whose x-axis labels need to be manually edited
# site_count_fxn | test_volume_timeline_fxn | top_gene_count_fxn | top_variant_count_fxn
site_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
top_gene_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
top_variant_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
gender_age_distribution_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
pt_mutation_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
specimen_type_stacked_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
tumor_purity_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
test_volume_timeline_fxn (DF = STAMP_DF, cohort="all", outdir = outdir, PerSite = FALSE)

# Per unique primary tumor site
#----------------------------------------------
# Specify key table
site.list <- data.frame(PrimaryTumorSite=unique(STAMP_DF$PrimaryTumorSite), stringsAsFactors = FALSE)
site.list$CohortName <- gsub("[(].*", "", site.list$PrimaryTumorSite)
site.list$CohortName <- gsub("[[:blank:]]and[[:blank:]]", " ", site.list$CohortName)

for (row_No in 1:nrow(site.list)) {
  s <- strsplit(site.list$CohortName[row_No], " ")[[1]] 
  site.list$CohortName[row_No] <- 
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
  
  remove(s)
}

site.list$CohortName <- gsub("[[:blank:]]", "", site.list$CohortName)
site.list <- site.list[order(site.list$PrimaryTumorSite),]

# Plots whose x-axis labels need to be manually edited: 
# top_gene_count_fxn | top_variant_count_fxn | test_volume_timeline_fxn
for (site_num in 1:nrow(site.list)) {
  cohort_id = site.list$CohortName[site_num]
  site_DF <- STAMP_DF[which(STAMP_DF$PrimaryTumorSite == site.list$PrimaryTumorSite[site_num]),]
  
  cat(paste(site_num,": ", cohort_id, sep=""),"\n")
  top_gene_count_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir)
  top_variant_count_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir)
  gender_age_distribution_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir)
  pt_mutation_count_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir)
  specimen_type_stacked_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir)
  tumor_purity_count_fxn (DF = site_DF, cohort=cohort_id,outdir = outdir)
  test_volume_timeline_fxn (DF = site_DF, cohort=cohort_id, outdir = outdir, PerSite = TRUE)
  
  remove(cohort_id,site_DF)
}
remove(row_No,s,site_num)

#################################
## Gene: Iterate through list
#################################
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_PerGene/"
if (!dir.exists(outdir)){dir.create(outdir)} 
cat("\n","\n",paste("Outdirectory: ", outdir, sep=""),"\n")

# Across all genes
#----------------------------------------------
# Plots whose x-axis labels need to be manually edited
# gene_count_fxn | top_site_count_fxn | variant_type_distribution_fxn
gene_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
top_site_count_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)
variant_type_distribution_fxn (DF = STAMP_DF, cohort="all", outdir = outdir)

# Per unique gene
#----------------------------------------------
# Specify key table
gene.list <- data.frame(Gene=unique(STAMP_DF$VariantGene), stringsAsFactors = FALSE)
gene.list <- gene.list[order(gene.list$Gene),]

# Plots whose x-axis labels need to be manually edited: 
# top_site_count_fxn | top_variant_count_fxn | variant_type_distribution_fxn
for (gene_num in 1:length(gene.list)) {
  gene_id = gene.list[gene_num]
  gene_DF <- STAMP_DF[which(STAMP_DF$VariantGene == gene_id),]
  
  cat(paste(gene_num,": ", gene_id, sep=""),"\n")
  top_site_count_fxn (DF = gene_DF, cohort=gene_id, outdir = outdir)
  top_variant_count_fxn (DF = gene_DF, cohort=gene_id, outdir = outdir)
  variant_type_distribution_fxn (DF = gene_DF, cohort=gene_id, outdir = outdir)
 
  remove(gene_id,gene_DF)
}
remove(gene_num)

#################################
## Mutation Plots 
#################################
outdir = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/TIFF_Lollipop/"
if (!dir.exists(outdir)){dir.create(outdir)} 
cat("\n","\n",paste("Outdirectory: ", outdir, sep=""),"\n")

# Classify pathogenicity statuses
pathogenic = c("Pathogenic", "Likely Pathogenic")
vus = c("Unknown significance", "Unknown")
benign = c("Likely Benign")

setwd("~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/")
source("mutation_hotspot.R")
remove(pathogenic,vus,benign)

sink()

# Delete temporary directory
if (dir.exists(tempdir)){unlink(tempdir, recursive = TRUE)} 

#################################
## Reference TABLES
#################################
outdir="~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/"

# TABLE: PrimaryTumorSite_No.Individuals
#----------------------------------------------
Site_List <- data.frame(PrimaryTumorSite = sort(unique(STAMP_DF$PrimaryTumorSite)),
                        stringsAsFactors = FALSE)
Site_List$No.Cases <- NA
for (row_No in 1:nrow(Site_List)) {
  site_ID = Site_List$PrimaryTumorSite[row_No]
  Site_List$No.Cases[row_No] = length(unique(STAMP_DF$PatientID[which(STAMP_DF$PrimaryTumorSite == site_ID)]))
  
  remove(site_ID)
}

write.table(Site_List, file=paste(outdir,"list_cancers.csv",sep=""),
            append = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# TABLE: Gene_No.Occurrences
#----------------------------------------------
Gene_Summary <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/STAMPEDE_GeneName_Description.csv", 
                         header = FALSE, stringsAsFactors = FALSE, sep = ",")
colnames(Gene_Summary) <- c("VariantGene","Summary")

GeneName_List <- data.frame(VariantGene = sort(unique(STAMP_DF$VariantGene)),
                            stringsAsFactors = FALSE)
GeneName_List$No.Occurrence <- NA
for (row_No in 1:nrow(GeneName_List)) {
  gene_ID = GeneName_List$VariantGene[row_No]
  GeneName_List$No.Occurrence[row_No] = nrow(STAMP_DF[which(STAMP_DF$VariantGene == gene_ID),])
  
  remove(gene_ID)
}
remove(row_No)

GeneName_List <- left_join(GeneName_List, Gene_Summary, by = "VariantGene")

write.table(GeneName_List, file=paste(outdir,"list_genes.csv",sep=""),
            append = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

remove(Site_List,GeneName_List,Gene_Summary)

#################################
## iFRAME CODE TABLE
#################################
Folder_root = "STAMPEDE"

# ALL
#----------------------------------------------
DF_misc <- data.frame(Plot_Name = c("site_count","test_volume","gene_count"),
                      Folder = Folder_root,
                      Label = "all",
                      iframe = c("2549","3322","3328"),
                      stringsAsFactors = FALSE)

# Per Gene
#----------------------------------------------
gene.uniqlist <- sort(unique(STAMP_DF$VariantGene))

DF_SiteCount <- data.frame(Plot_Name = "top_site_count", 
                           Folder = paste(Folder_root,"/Top_Site_Count",sep=""),
                           Label = gene.uniqlist,
                           iframe = NA,
                           iter=1,
                           stringsAsFactors = FALSE)

DF_VariantCount_Gene <- data.frame(Plot_Name = "top_variant_count", 
                                   Folder = paste(Folder_root,"/Top_Variant_Count",sep=""),
                                   Label = gene.uniqlist,
                                   iframe = NA,
                                   iter=2,
                                   stringsAsFactors = FALSE)

DF_PerGene <- rbind(DF_SiteCount,DF_VariantCount_Gene)
DF_PerGene <- DF_PerGene[order(DF_PerGene$iter, decreasing = FALSE),]
DF_PerGene <- DF_PerGene[order(DF_PerGene$Label),]

gene_missing <- c("CCND1","CCND2","CDK12","CREBBP","CUL3","DDR2","EP300","EZH2",
                  "FLT3","HNF1A","MED12","NTRK2","PIK3R1","POLD1","POLE","RAF1",
                  "RIT1","RPS4Y1","SETBP1","SETD2","SOX2","SPOP","SRC")
DF_PerGene$iframe[which(DF_PerGene$Label %in% gene_missing & DF_PerGene$Plot_Name == "top_variant_count")] <- "0"

DF_PerGene <- rbind(data.frame(Plot_Name = "top_site_count", 
                 Folder = paste(Folder_root,"/Top_Site_Count",sep=""),
                 Label = "all",
                 iframe = "3331",
                 iter=1,
                 stringsAsFactors = FALSE), DF_PerGene)

for (row_No in 2:nrow(DF_PerGene)) {
  # Input iframe only if plot exists
  if (is.na(DF_PerGene$iframe[row_No])) {
    
    # Indicate start of iteration
    if (isTRUE(row_No == 2)) {
      DF_PerGene$iframe[row_No] <- as.numeric("3333")
      
    } else {
      if (isTRUE(DF_PerGene$iframe[row_No -1] == "0")) {
        DF_PerGene$iframe[row_No] <- as.numeric(DF_PerGene$iframe[row_No -2]) +2
      } else {
        DF_PerGene$iframe[row_No] <- as.numeric(DF_PerGene$iframe[row_No -1]) +2  
      }
    }
  }
}
DF_PerGene[DF_PerGene == 0] <- NA

DF_PerGene <- DF_PerGene[order(DF_PerGene$iter, decreasing = FALSE),]
DF_PerGene <- DF_PerGene[,1:4]

DF_lollipop <- data.frame(Plot_Name = "gene_lollipop", 
                          Folder = paste(Folder_root,"/Gene_Lollipop",sep=""),
                          Label = gene.uniqlist,
                          iframe = NA,
                          stringsAsFactors = FALSE)
for (row_No in 1:nrow(DF_lollipop)) {DF_lollipop$iframe[row_No] <- seq(3789,4037,2)[row_No] }

DF_VariantType_Gene <- data.frame(Plot_Name = "var_type_distribution", 
                                  Folder = paste(Folder_root,"/Var_Type_Distribution",sep=""),
                                  Label = c("all", gene.uniqlist),
                                  iframe = NA,
                                  stringsAsFactors = FALSE)
end_No = which(gene.uniqlist == "ALK")
end_No_2 = which(gene.uniqlist == "AR")
for (row_No in 1:nrow(DF_VariantType_Gene)) {
  if (isTRUE(DF_VariantType_Gene$Label[row_No] == "all")) {
    DF_VariantType_Gene$iframe[row_No] <- "4245"
  } else if (isTRUE(row_No <= end_No +1)) {
    DF_VariantType_Gene$iframe[row_No] <- seq(4247,4251,2)[row_No -1]
  } else if (isTRUE(DF_VariantType_Gene$Label[row_No] == "APC")) {
    DF_VariantType_Gene$iframe[row_No] <- "4253"
  } else if (isTRUE(DF_VariantType_Gene$Label[row_No] == "AR")) {
    DF_VariantType_Gene$iframe[row_No] <- "4256"
  } else {
    DF_VariantType_Gene$iframe[row_No] <- seq(4259,4497,2)[row_No -end_No_2 -1]  
  }
}
  
DF_PerGene <- rbind(DF_PerGene,DF_lollipop,DF_VariantType_Gene)
remove(DF_SiteCount,DF_VariantCount_Gene,DF_lollipop,DF_VariantType_Gene)

# Per PrimaryTumorSite
#----------------------------------------------
DF_GeneCount <- data.frame(Plot_Name = "top_gene_count", 
                           Folder = paste(Folder_root,"/Top_Gene_Count",sep=""),
                           Label = append("all", unlist(site.list$PrimaryTumorSite)),
                           iframe = NA,
                           stringsAsFactors = FALSE)
DF_GeneCount$iframe[1] <- "3324"


DF_VariantCount <- data.frame(Plot_Name = "top_variant_count", 
                              Folder = paste(Folder_root,"/Top_Variant_Count",sep=""),
                              Label = append("all", unlist(site.list$PrimaryTumorSite)),
                              iframe = NA,
                              stringsAsFactors = FALSE)
DF_VariantCount$iframe[1] <- "2557"
site_missing <- 
  site.list$PrimaryTumorSite[!(site.list$CohortName %in% 
                                 c("AmpullaOfVater","Breast","CentralNervousSystem","ColonRectum",
                                   "Esophagus","Gallbladder","IntrahepaticBileDucts","LiverHepatocellular",
                                   "Lung","LymphNode","Pancreas","Peritoneum","ProstateGland","Skin",
                                   "SmallIntestine","SoftTissue","Stomach","Thyroid","Uvea"))]
DF_VariantCount$iframe[which(DF_VariantCount$Label %in% site_missing)] <- "0"


DF_GenderAge <- data.frame(Plot_Name = "gender_age_distribution", 
                           Folder = paste(Folder_root,"/Gender_Age_Distribution",sep=""),
                           Label = append("all", unlist(site.list$PrimaryTumorSite)),
                           iframe = NA,
                           stringsAsFactors = FALSE)
DF_GenderAge$iframe[1] <- "2560"


DF_Mutation.all <- data.frame(Plot_Name = "pt_mutation_count_allvariants", 
                              Folder = paste(Folder_root,"/Patient_Mutation_Count/All_Variants",sep=""),
                              Label = append("all", unlist(site.list$PrimaryTumorSite)),
                              iframe = NA,
                              stringsAsFactors = FALSE)
DF_Mutation.all$iframe[1] <- "2564"



DF_Mutation.patho <- data.frame(Plot_Name = "pt_mutation_count_pathovariants", 
                                Folder = paste(Folder_root,"/Patient_Mutation_Count/Pathogenic_Variants",sep=""),
                                Label = append("all", unlist(site.list$PrimaryTumorSite)),
                                iframe = NA,
                                stringsAsFactors = FALSE)
DF_Mutation.patho$iframe[1] <- "2567"
site_missing <- c("back","foot","trachea")
DF_Mutation.patho$iframe[which(DF_Mutation.patho$Label %in% site_missing)] <- "0"


DF_Mutation.VUS <- data.frame(Plot_Name = "pt_mutation_count_VUS", 
                              Folder = paste(Folder_root,"/Patient_Mutation_Count/UnknownSignificance_Variants",sep=""),
                              Label = append("all", unlist(site.list$PrimaryTumorSite)),
                              iframe = NA,
                              stringsAsFactors = FALSE)
DF_Mutation.VUS$iframe[1] <- "2570"
site_missing <- "mesenteric mass"
DF_Mutation.VUS$iframe[which(DF_Mutation.VUS$Label %in% site_missing)] <- "0"


DF_SpecimenCount <- data.frame(Plot_Name = "specimen_type_count", 
                               Folder = paste(Folder_root,"/Specimen_Type_Count",sep=""),
                               Label = append("all", unlist(site.list$PrimaryTumorSite)),
                               iframe = NA,
                               stringsAsFactors = FALSE)
DF_SpecimenCount$iframe[1] <- "2573"


DF_TumorPurity <- data.frame(Plot_Name = "tumor_purity_count", 
                             Folder = paste(Folder_root,"/Tumor_Purity_Count",sep=""),
                             Label = append("all", unlist(site.list$PrimaryTumorSite)),
                             iframe = NA,
                             stringsAsFactors = FALSE)
DF_TumorPurity$iframe[1] <- "2576"

DF_iframe <- data.frame(Label = DF_GeneCount$Label,
                        iframe_GeneCt = DF_GeneCount$iframe,
                        iframe_VariantCt = DF_VariantCount$iframe,
                        iframe_GenderAge = DF_GenderAge$iframe,
                        iframe_All = DF_Mutation.all$iframe,
                        iframe_Patho = DF_Mutation.patho$iframe,
                        iframe_VUS = DF_Mutation.VUS$iframe,
                        iframe_Specimen = DF_SpecimenCount$iframe,
                        iframe_Purity = DF_TumorPurity$iframe,
                        stringsAsFactors = FALSE)

for (row_No in 2:nrow(DF_iframe)) {
  for (col_No in 2:ncol(DF_iframe)) {
    
    # Input iframe only if plot exists
    if (is.na(DF_iframe[row_No,col_No])) {
      
      # Indicate start of iteration
      if (isTRUE(row_No == 2 & col_No == 2)) {
        DF_iframe[row_No, col_No] <- as.numeric("2580")
        
      } else if (isTRUE(col_No == 2)) {
        DF_iframe[row_No, col_No] <- as.numeric(DF_iframe[row_No -1, ncol(DF_iframe)]) +2
        
      } else {
        if (DF_iframe[row_No, col_No -1] == "0") {
          DF_iframe[row_No, col_No] <- as.numeric(DF_iframe[row_No, col_No -2]) +2
        } else {
          DF_iframe[row_No, col_No] <- as.numeric(DF_iframe[row_No, col_No -1]) +2  
        }
      }
    }
    remove(col_No)
  }
}
DF_iframe[DF_iframe == 0] <- NA

for (row_No in 1:nrow(DF_iframe)) {
  DF_GeneCount$iframe[row_No] <- as.numeric(DF_iframe$iframe_GeneCt[row_No])
  DF_VariantCount$iframe[row_No] <- as.numeric(DF_iframe$iframe_VariantCt[row_No])
  DF_GenderAge$iframe[row_No] <- as.numeric(DF_iframe$iframe_GenderAge[row_No])
  DF_Mutation.all$iframe[row_No] <- as.numeric(DF_iframe$iframe_All[row_No])
  DF_Mutation.patho$iframe[row_No] <- as.numeric(DF_iframe$iframe_Patho[row_No])
  DF_Mutation.VUS$iframe[row_No] <- as.numeric(DF_iframe$iframe_VUS[row_No])
  DF_SpecimenCount$iframe[row_No] <- as.numeric(DF_iframe$iframe_Specimen[row_No])
  DF_TumorPurity$iframe[row_No] <- as.numeric(DF_iframe$iframe_Purity[row_No])
}

DF_Tumor_TestVolume <- data.frame(Plot_Name = "test_volume", 
                             Folder = paste(Folder_root,"/PerTumor_OrderVolume",sep=""),
                             Label = unlist(site.list$PrimaryTumorSite),
                             iframe = NA,
                             stringsAsFactors = FALSE)
for (row_No in 1:nrow(DF_Tumor_TestVolume)) {DF_Tumor_TestVolume$iframe[row_No] <- seq(4144,4242,2)[row_No] }

DF_iframe_FULL <- rbind(DF_misc,DF_GeneCount,DF_VariantCount,DF_GenderAge,
                        DF_Mutation.all,DF_Mutation.patho,DF_Mutation.VUS,
                        DF_SpecimenCount,DF_TumorPurity,DF_PerGene,DF_Tumor_TestVolume)

DF_tumor <- DF_iframe_FULL[which(DF_iframe_FULL$Label %in% unlist(site.list$PrimaryTumorSite)),]
DF_tumor <- DF_tumor[order(DF_tumor$Label, decreasing = FALSE),]

DF_gene <- DF_iframe_FULL[which(DF_iframe_FULL$Label %in% gene.list),]
DF_gene <- DF_gene[order(DF_gene$Label, decreasing = FALSE),]

DF_all <- DF_iframe_FULL[which(DF_iframe_FULL$Label == "all"),]
DF_all <- DF_all[order(DF_all$Folder, decreasing = FALSE),]

DF_iframe_FULL_ordered <- rbind(DF_all,DF_tumor,DF_gene)

## Write to local computer
#----------------------------------------------
write.table(DF_iframe_FULL_ordered, 
            file="~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/site_mappings.csv",
            append = FALSE, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

remove(DF_misc,DF_GeneCount,DF_VariantCount,DF_GenderAge,DF_Mutation.all,DF_Mutation.patho,
       DF_Mutation.VUS,DF_SpecimenCount,DF_TumorPurity,DF_PerGene,DF_iframe,DF_Tumor_TestVolume,
       DF_iframe_FULL,site.list,Folder_root,gene.list,gene.uniqlist,row_No,site_missing,gene_missing,
       end_No,end_No_2,DF_all,DF_tumor,DF_gene)
