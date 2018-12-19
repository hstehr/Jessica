## Generate lollipop plots for SNV, Frameshift, Indels, Insertions, Deletions, Duplications
## Input: "Mutation_Hotspot/20181114_syapse_export_DF_STAMP_4Map.csv" = 9592 total STAMP entries

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/STAMP/")

# Load Library
#----------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)

# Lollipop Plot
#----------------------------------------------
Mutation_Plot <- function(protein, domains, variant_data, 
                          yaxis_max, y_tick, y_legend, x_tick,
                          variant_type, col="gray33"){
  
  # Visualization options
  voff = 0
  seg_size = c(6, 10)
  
  plot <- ggplot(data=variant_data) +
    # Add points
    geom_segment(aes(x=x, xend=x, y=voff, yend=y), color="gray88", alpha=0.8) +
    geom_point(aes(x=x, y=y), fill=col, colour="gray88", size=3, shape=21, alpha=0.8) +
    
    # Plot protein
    geom_segment(data=protein, aes(x=start, xend=end, y=voff-y_legend, yend=voff-y_legend), 
                 color="gray88", size=seg_size[1]) +
    # Add domain
    geom_segment(data=domains, aes(x=start, xend=end, y=voff-y_legend, yend=voff-y_legend, 
                                   color=name), size=seg_size[1]) +
    
    # Labels
    ylab("# Mutations") +
    labs(title = variant_type, 
         subtitle = paste("N = ", sum(variant_data$y), sep=""),
         color = "Domains") +
    
    # Scaling
    scale_fill_brewer(palette="Set2") +
    scale_x_continuous(breaks=seq(0, protein$end, x_tick)) +
    scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c((-1*y_legend), yaxis_max)) +
    
    # Theme
    theme_light() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust=0, vjust=0, face="bold", size=14),
          plot.subtitle = element_text(hjust=1, face="bold")
    )
}

# Lollipop PIPELINE
#----------------------------------------------
Mutation_Pipeline <- function (gene.list, DF, assay, assayName,
                               DF_genes, DF_domain) {
  for (gene_num in 1:length(gene.list)) {
    
    # Specify gene INFO
    gene_id = gene.list[gene_num]
    
    protein = data.frame(symbol=gene_id, start=1, end=DF_genes[DF_genes$Name == gene_id, "AA.Length"])
    
    row_id = which(DF_domain$Name == gene_id)
    domain_start = DF_domain[row_id, "Domain_Start"]
    domain_end = DF_domain[row_id, "Domain_End"]
    domain_name = DF_domain[row_id, "Domain"]
    domains = data.frame(start=domain_start, end=domain_end, name=domain_name)
    
    remove(row_id)
    remove(domain_name, domain_start, domain_end)
    
    # Subset genes of interest from STAMP DF
    Curr_gene <- DF[DF$base.gene == gene_id,]
    
    # Remove cases with variant position > AA length of transcript
    var.remove <- which(Curr_gene$var.position > protein$end)
    DF_remove <- Curr_gene[var.remove,]
    if (nrow(DF_remove) != 0){
      Curr_gene <- Curr_gene[!Curr_gene$sys.label %in% DF_remove$sys.label,]
    }
    remove(DF_remove, var.remove)
    
    # Extract variant INFO from smpl.hgvsProtein column
    variant.loc <- Curr_gene$var.position
    
    # Extract pathogenicity status from smpl.pathogenicityStatus column
    variant.pathogenicity = Curr_gene$smpl.pathogenicityStatus
    
    # Create variant location-pathogenicity DF
    Curr_variant = data.frame(allele=as.numeric(variant.loc), significance=variant.pathogenicity)
    
    # Pathogenicity statuses present in DF
    index.sig = unique(variant.pathogenicity)
    pList = list()
    
    # Generate lollipop plots for each pathogenicity status
    loc.freq.max = max(as.numeric(as.vector(sort(table(as.factor(variant.loc))))))
    
    # Specify y-axis parameters 
    if (loc.freq.max <= 10) {
      yaxis_max = 10
      y_tick = 1
      y_legend = 1
    } else if (loc.freq.max <= 20) {
      yaxis_max = 20
      y_tick = 2
      y_legend = 1.25
    } else if (loc.freq.max <= 50) {
      yaxis_max = 50
      y_tick = 5
      y_legend = 5
    } else if (loc.freq.max <= 100) {
      yaxis_max = 100
      y_tick = 10
      y_legend = 10
    } else if (loc.freq.max <= 250) {
      yaxis_max = 250
      y_tick = 25
      y_legend = 25
    } else {
      yaxis_max = loc.freq.max
      y_tick = 50
      y_legend = 50
    }
    
    # Specify x-axis parameters 
    if (as.numeric(protein$end) <= 800) {
      x_tick = 25
    } else if (as.numeric(protein$end) <= 1400) {
      x_tick = 50
    } else if (as.numeric(protein$end) <= 2400) {
      x_tick = 100
    } else if (as.numeric(protein$end) > 2400) {
      x_tick = 200
    }
    
    # Generate plots for each pathogenicity status
    if (as.numeric(table(DF$base.gene == gene_id)["TRUE"]) > 0) {
      allele.count = table(variant.loc)
      pAll <- Mutation_Plot(protein = protein, domains = domains, 
                            variant_data = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count)), 
                            yaxis_max = yaxis_max, y_tick = y_tick, y_legend = y_legend, x_tick = x_tick, 
                            variant_type = "All Variants", 
                            col = "black")
      pNum = 1
      pList[[pNum]] <- pAll
      
      if ("Pathogenic" %in% index.sig | "Likely Pathogenic" %in% index.sig) {
        allele.count = table(variant.loc[which(variant.pathogenicity %in% pathogenic)])
        pPathogenic <- Mutation_Plot(protein = protein, domains = domains, 
                                     variant_data = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count)),  
                                     yaxis_max = yaxis_max, y_tick = y_tick, y_legend = y_legend, x_tick = x_tick, 
                                     variant_type = "Pathogenic", 
                                     col = "firebrick4")
      }
      if (exists("pPathogenic")) {
        pNum = pNum + 1
        pList[[pNum]] <- pPathogenic
      }
      
      if ("Unknown significance" %in% index.sig | "Unknown" %in% index.sig) {
        allele.count = table(variant.loc[which(variant.pathogenicity %in% vus)])
        pVUS <- Mutation_Plot(protein = protein, domains = domains, 
                              variant_data = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count)),  
                              yaxis_max = yaxis_max, y_tick = y_tick, y_legend = y_legend, x_tick = x_tick, 
                              variant_type = "Unknown Significance", 
                              col = "darkslategrey")
      }
      if (exists("pVUS")) {
        pNum = pNum + 1
        pList[[pNum]] <- pVUS
      }
      
      if ("Likely Benign" %in% index.sig) {
        allele.count = table(variant.loc[which(variant.pathogenicity %in% benign)])
        pBenign <- Mutation_Plot(protein = protein, domains = domains, 
                                 variant_data = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count)), 
                                 yaxis_max = yaxis_max, y_tick = y_tick, y_legend = y_legend, x_tick = x_tick, 
                                 variant_type = "Benign", 
                                 col = "darkgreen")
      }
      if (exists("pBenign")) {
        pNum = pNum + 1
        pList[[pNum]] <- pBenign
      }
      
      # Merge all plots into single frame
      gpanels <- ggarrange(plotlist = pList,
                           nrow = 4, ncol = 1,
                           legend = "bottom", common.legend = TRUE)
      #Save plot to local computer
      gpanels <- annotate_figure(gpanels,
                                 top = text_grob(paste(gene_id, " variants in ", assayName, " (", protein$end, " aa)", sep=""), 
                                                 face = "bold", size = 18))
      file_id = paste("Mutation_Hotspot/LollipopPlots/", gene_id, "_", assay, ".jpg", sep="")
      ggexport(gpanels, filename=file_id, width = 4000, height = 4000, res=350)
      
      remove(pAll, gpanels)
      if (exists("pPathogenic")) {remove(pPathogenic)}
      if (exists("pVUS")) {remove(pVUS)}
      if (exists("pBenign")) {remove(pBenign)}
    }
    
    remove(Curr_gene, Curr_variant, domains, protein)
    remove(allele.count, variant.loc, variant.pathogenicity)
    remove(x_tick, y_legend, y_tick, yaxis_max, loc.freq.max)
    remove(file_id, gene_id, gene_num, index.sig)
    remove(pList, pNum)
  }
}

# Classify pathogenicity statuses
#----------------------------------------------
pathogenic = c("Pathogenic", "Likely Pathogenic")
vus = c("Unknown significance", "Unknown")
benign = c("Likely Benign")

# Read relevant files
#----------------------------------------------
# Protein length INFO
DF_genes <- read.csv(file = "Mutation_Hotspot/STAMPv2_Annotation_v2_genes.csv",
                     header = TRUE,
                     na.strings = c(""," ","NA"),
                     stringsAsFactors = FALSE,
                     sep = ",")
DF_genes <- DF_genes[complete.cases(DF_genes$Name), ]
DF_genes <- DF_genes[complete.cases(DF_genes$AA.Length), ]

# Protein domain INFO
DF_domain <- read.csv(file = "Mutation_Hotspot/STAMPv2_Annotation_v2_domains.csv",
                      header = TRUE,
                      na.strings = c(""," ","NA"),
                      stringsAsFactors = FALSE,
                      sep = ",")
DF_domain <- DF_domain[complete.cases(DF_domain$Name), ]
DF_domain <- DF_domain[complete.cases(DF_domain$Domain), ]

DF <- read.csv(file = "Mutation_Hotspot/20181114_syapse_export_DF_STAMP_4Map.csv",
               header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")

# Structure DF columns 
#----------------------------------------------
personnel <- c("smpl.hasOrderingPhysician", "sys.uniqueId", "smpl.gender")
assay <- c("smpl.assayName", "base.chromosome")
specimen <- c("smpl.specimenType", "smpl.specimenSite")
col.change <- c(personnel, assay, specimen)
DF[ ,col.change] <- lapply(DF[ ,col.change], factor)
remove(personnel, assay, specimen, col.change)

# Check missing gene INFO
#----------------------------------------------
gene.DF = data.frame(gene=sort(unique(DF_genes$Name)), stringsAsFactors=FALSE)
domain.DF = data.frame(gene=sort(unique(DF_domain$Name)), stringsAsFactors=FALSE)

gene.list <- data.frame(gene=unique(DF$base.gene), stringsAsFactors=FALSE)
print(paste("Transcripts with missing amino acid length:", sep=""))
print(anti_join(gene.list, gene.DF, by = "gene"))
print(paste("Transcripts with missing domain information:", sep=""))
print(anti_join(gene.list, domain.DF, by = "gene"))

remove(gene.DF, domain.DF)

# Categorize by assay data
DF_STAMP_v1 <- DF[DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" |
                    DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)", ]
Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_v1$base.gene)), 
                  DF = DF_STAMP_v1, assay = "STAMP_v1", 
                  assayName = "STAMP v1 Panel",
                  DF_genes = DF_genes, DF_domain = DF_domain)'
'

DF_STAMP_v2 <- DF[DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)", ]
Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_v2$base.gene)), 
                  DF = DF_STAMP_v2, assay = "STAMP_v2", 
                  assayName = "STAMP v2 Panel",
                  DF_genes = DF_genes, DF_domain = DF_domain)

DF_STAMP_all <- DF
Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_all$base.gene)), 
                  DF = DF, assay = "STAMP_all", 
                  assayName = "STAMP v1 and v2 Panels",
                  DF_genes = DF_genes, DF_domain = DF_domain)

remove(DF_STAMP_v1,DF_STAMP_v2,DF_STAMP_all)
