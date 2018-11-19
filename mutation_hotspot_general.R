###################################
# Load Library
###################################
library(rentrez)
library(ggplot2)
library(stringr)
library(ggpubr)
library(MASS)
library(viridis)
library(dplyr)

###################################
# Functions
###################################
# Extract variant position
#----------------------------------------------
PositionExtract <- function(DF,DF_sub,var_type) {
  
  DF_inter <- DF[DF$sys.label %in% DF_sub$variant,]
  
  # Extract variant positions
  #----------------------------------------------
  # Subset cases missing smpl.hgvsProtein field
  phgvs.na <- which(is.na(DF_inter$smpl.hgvsProtein))
  DF_NAprotein <- DF_inter[phgvs.na,]
  
  if (nrow(DF_NAprotein) != 0){
    DF_inter <- DF_inter[!DF_inter$smpl.hgvsCoding %in% DF_NAprotein$smpl.hgvsCoding,]
    assign("DF_NAprotein", DF_NAprotein, envir = .GlobalEnv)
  }
  
  DF_inter$smpl.hgvsProtein <- gsub("[[:space:]]*$", "", DF_inter$smpl.hgvsProtein)
  
  DF_inter$var.position <- gsub("^p.[[:alpha:]]{3}", "", DF_inter$smpl.hgvsProtein)
  if (var_type == "SNV") {
    DF_inter$var.position <- gsub("[[:alpha:]]{3}$", "", DF_inter$var.position)
  } else {
    DF_inter$var.position <- gsub("([[:digit:]]+)((_)*.*)", "\\1", DF_inter$var.position)
  }
  DF_inter$var.position <- as.numeric(DF_inter$var.position)
  
  DF_inter$aa.start <- gsub("^p.", "", DF_inter$smpl.hgvsProtein)
  if (var_type == "SNV") {
    DF_inter$aa.start <- gsub("[[:digit:]]+[[:alpha:]]{,3}$", "", DF_inter$aa.start)
  } else if (var_type == "Duplications" | var_type == "Frameshift") {
    DF_inter$aa.start <- gsub("[[:digit:]]+[[:alnum:]]+", "", DF_inter$aa.start)
    if (var_type == "Duplications") {
      DF_inter$aa.start <- gsub("_[[:alpha:]]+", "", DF_inter$aa.start)
    }
  } else {
    DF_inter$aa.start <- gsub("[[:digit:]]+.*", "", DF_inter$aa.start)
  }
  
  DF_inter$aa.end <- gsub("^p.[[:alpha:]]{3}[[:digit:]]+", "", DF_inter$smpl.hgvsProtein)
  if (var_type != "SNV") {
    DF_inter$aa.end <- gsub("(^_)*([[:alpha:]]{3})([[:alnum:]]+)", "\\2", DF_inter$aa.end)
    if (var_type == "Indels") {
      DF_inter$aa.end <- gsub("(^_[[:digit:]]+)*(del)(.*)" , "\\2", DF_inter$aa.end)
      }}
  
  assign("DF_inter", DF_inter, envir = .GlobalEnv)
}

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
    )}

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

# Read gene INFO files
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

#########################################################################################################
# STAMP - Solid Tumor Actionable Mutation Panel
###################################
DF_Full <- read.csv(file = "20181018_syapse_export_all_variants_patientNameAndMrnRemoved.csv",
                    header = TRUE,
                    na.strings = c(""," ","NA"),
                    stringsAsFactors = FALSE,
                    sep = ",")
# Shorten column names of DF
colnames(DF_Full) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_Full))

DF_STAMP_Full <- DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)", ]

DF = DF_STAMP_Full
remove(DF_STAMP_Full, DF_Full)

# Structure DF columns 
#----------------------------------------------
personnel <- c("smpl.hasOrderingPhysician", "sys.uniqueId", "smpl.gender")
assay <- c("smpl.assayName", "base.chromosome")
specimen <- c("smpl.specimenType", "smpl.specimenSite")
col.change <- c(personnel, assay, specimen)
DF[ ,col.change] <- lapply(DF[ ,col.change], factor)
remove(personnel, assay, specimen, col.change)

# Categorize protein sequence variants
#----------------------------------------------
# Subset columns of interest for simplification
# Beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron
# End of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron
DF_var <- data.frame(variant=DF$sys.label, gene=DF$base.gene,
                     coding=DF$smpl.hgvsCoding, protein=DF$smpl.hgvsProtein,
                     
                     intron=grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+.*[[:blank:]]*$", DF$smpl.hgvsCoding),
                     
                     synonymous=grepl("=",DF$smpl.hgvsProtein),
                     upstream=grepl("^c.-[[:digit:]]+[ATCG]>[ATCG][[:blank:]]*$", DF$smpl.hgvsCoding),
                     
                     SNVprotein=grepl("^p.[[:upper:]]{1}[[:lower:]]{2}[[:digit:]]+[[:upper:]]{1}[[:lower:]]{2}[[:blank:]]*$", 
                                      DF$smpl.hgvsProtein),
                     SNVcoding=grepl("^c.[[:digit:]]+[ATCG]>[ATCG][[:blank:]]*$", DF$smpl.hgvsCoding),
                     
                     fs=grepl("fs", DF$smpl.hgvsProtein),
                     indel=grepl("del.*ins", DF$sys.label),
                     ins=grepl("ins", DF$smpl.hgvsCoding),
                     del=grepl("del", DF$smpl.hgvsCoding),
                     dup=grepl("dup", DF$smpl.hgvsCoding)
)
DF_var$gene <- as.character(DF_var$gene)
DF_var$variant <- as.character(DF_var$variant)

DF_var$remain <- NA
for (i in 1:nrow(DF_var)) {
  if (length(grep(FALSE, DF_var[i,5:14])) == 10) {
    DF_var$remain[i] <- TRUE
  } else {
    DF_var$remain[i] <- FALSE
  }}
remove(i)

# Subset based on matching of strings in smpl.hgvsCoding or smpl.hgvsProtein
# Classification order: SNV > Frameshift > indel > insertion > deletion > duplication

# "p.=" indicates protein has not been analysed, RNA was, but no change is expected
DF_synonymous <- DF_var[DF_var$synonymous == TRUE, 1:4]
DF_upstream <- DF_var[DF_var$upstream == TRUE, 1:4]
DF_intron <- DF_var[DF_var$intron == TRUE, 1:4]

DF_SNV <- DF_var[DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE &
                   DF_var$synonymous == FALSE, 1:4]
DF_Frameshift <- DF_var[DF_var$fs == TRUE & 
                          DF_var$synonymous == FALSE & DF_var$upstream == FALSE & 
                          DF_var$intron == FALSE, 1:4]

DF_indel <- DF_var[DF_var$indel == TRUE & 
                     DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                     DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE & 
                     DF_var$fs == FALSE , 1:4]
DF_ins <- DF_var[DF_var$ins == TRUE & 
                   DF_var$indel == FALSE & 
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE & 
                   DF_var$fs == FALSE , 1:4]
DF_del <- DF_var[DF_var$del == TRUE & 
                   DF_var$ins == FALSE & DF_var$indel == FALSE & 
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE & 
                   DF_var$fs == FALSE , 1:4]
DF_dup <- DF_var[DF_var$dup == TRUE & 
                   DF_var$del == FALSE & DF_var$ins == FALSE & DF_var$indel == FALSE & 
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE & 
                   DF_var$fs == FALSE , 1:4]

DF_remain <- DF_var[DF_var$remain == TRUE, 1:4]
remove(DF_var)

print(paste("All cases have been categorized based on type of nucleotide change, as indicated:", 
            nrow(DF) == (nrow(DF_synonymous) + nrow(DF_upstream) + nrow(DF_intron) +
                           nrow(DF_SNV) + nrow(DF_Frameshift) + nrow(DF_indel) + 
                           nrow(DF_ins) + nrow(DF_del) + nrow(DF_dup) + nrow(DF_remain))), 
      sep="")

# Unique genes per variant category
#----------------------------------------------
gene.no <- function(DF, gene_type) {
  print(paste(gene_type, " mutations: n=", 
              length(unique(DF$gene)), " genes and n=", length(DF$gene), " cases)", sep=""))
}

print("Synonymous, upstream & intronic mutations are not mapped to lollipop plots")
remove(DF_synonymous,DF_upstream, DF_intron,DF_remain)

DF_NAprotein_Full <- data.frame()
DF_STAMP_4Map <- data.frame()

PositionExtract(DF = DF, DF_sub = DF_SNV, var_type = "SNV")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)

PositionExtract(DF = DF, DF_sub = DF_Frameshift, var_type = "Frameshift")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)

PositionExtract(DF = DF, DF_sub = DF_indel, var_type = "Indels")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)
DF_NAprotein_Full <- rbind(DF_NAprotein_Full, DF_NAprotein)

PositionExtract(DF = DF, DF_sub = DF_ins, var_type = "Insertions")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)
DF_NAprotein_Full <- rbind(DF_NAprotein_Full, DF_NAprotein)

PositionExtract(DF = DF, DF_sub = DF_del, var_type = "Deletions")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)
DF_NAprotein_Full <- rbind(DF_NAprotein_Full, DF_NAprotein)

PositionExtract(DF = DF, DF_sub = DF_dup, var_type = "Duplications")
DF_STAMP_4Map <- rbind(DF_STAMP_4Map, DF_inter)
DF_NAprotein_Full <- rbind(DF_NAprotein_Full, DF_NAprotein)

remove(DF_SNV, DF_Frameshift,DF_indel,DF_ins,DF_del,DF_dup)
remove(DF_inter,DF_NAprotein)

DF <- DF_STAMP_4Map

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
                  DF_genes = DF_genes, DF_domain = DF_domain)

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

#########################################################################################################
