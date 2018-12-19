## Subset 20181113_syapse_export_all_variants_QC.csv based on mutation type = 9918 total STAMP entries
## Specific categories of mutations
## Synonymous, Upstream, Intronic, SNV, Frameshift, Indels, Insertions, Deletions, Duplications
## Synonymous, upstream & intronic mutations are not mapped to lollipop plots = 305 entries
## Output: "Mutation_Hotspot/20181114_syapse_export_DF_STAMP_4Map.csv" = 9592 entries
## Output: "Mutation_Hotspot/20181114_syapse_export_DF_NAprotein.csv" >> missing smpl.hgvsProtein = 21 entries

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/")

# Load Library
#----------------------------------------------
library("dplyr")
library("eeptools")

# Extract variant positions FUNCTION
#----------------------------------------------
Extract_VarPosition <- function(DF) {
  DF$var.position <- gsub("^p.[[:alpha:]]{3}", "", DF$smpl.hgvsProtein)
  
  for (row_No in 1:nrow(DF)) {
    
    if (DF$var.type[row_No] == "SNV") {
      DF$var.position[row_No] <- gsub("[[:alpha:]]{3}$", "", DF$var.position[row_No])
      DF$var.position[row_No] <- gsub("(ThrextTer37)$", "", DF$var.position[row_No])
    } else if (DF$var.type[row_No] == "Synonymous") {
      DF$var.position[row_No] <- gsub("([[:digit:]]+)([[:punct:]]+.*)", "\\1", DF$var.position[row_No])
    } else {
      DF$var.position[row_No] <- gsub("([[:digit:]]+)((_)*.*)", "\\1", DF$var.position[row_No])
    }
    DF$var.position[row_No] <- as.numeric(DF$var.position[row_No])
    
    DF$aa.start[row_No] <- gsub("^p.", "", DF$smpl.hgvsProtein[row_No])
    if (DF$var.type[row_No] == "SNV") {
      DF$aa.start[row_No] <- gsub("[[:digit:]]+[[:alpha:]]{,3}$", "", DF$aa.start[row_No])
      DF$aa.start[row_No] <- gsub("(840ThrextTer)$", "", DF$aa.start[row_No])
    } else if (DF$var.type[row_No] == "Duplication" | DF$var.type[row_No] == "Frameshift") {
      DF$aa.start[row_No] <- gsub("[[:digit:]]+[[:alnum:]]+", "", DF$aa.start[row_No])
      if (DF$var.type[row_No] == "Duplication") {
        DF$aa.start[row_No] <- gsub("_[[:alpha:]]+", "", DF$aa.start[row_No])
      }
    } else {
      DF$aa.start[row_No] <- gsub("([[:alpha:]]+)([[:digit:]]+.*)", "\\1", DF$aa.start[row_No])
    }
    
    DF$aa.end[row_No] <- gsub("^p.[[:alpha:]]{3}[[:digit:]]+", "", DF$smpl.hgvsProtein[row_No])
    if (DF$var.type[row_No] == "SNV") {
      DF$aa.end[row_No] <- gsub("(extTer37)$", "", DF$aa.end[row_No])
    } else if (DF$var.type[row_No] == "Frameshift") {
      DF$aa.end[row_No] <- gsub("([[:alpha:]]{3})(.*)", "\\1", DF$aa.end[row_No])
    } else if (DF$var.type[row_No] == "Indel" | DF$var.type[row_No] == "Insertion" | 
               DF$var.type[row_No] == "Deletion" |DF$var.type[row_No] == "Duplication") {
      DF$aa.end[row_No] <- gsub("(^_)([[:alpha:]]{3})(.*)", "\\2",DF$aa.end[row_No])
      
      if (DF$var.type[row_No] == "Indel") {
        DF$aa.end[row_No] <- gsub("(.*delins)([[:alpha:]]{3})(.*)", "\\2",DF$aa.end[row_No])
      } else if (DF$var.type[row_No] == "Insertion") {
        DF$aa.end[row_No] <- gsub("(^_774ins)([[:alpha:]]{3})(.*)", "\\2",DF$aa.end[row_No])
      } else if (DF$var.type[row_No] == "Deletion") {
        DF$aa.end[row_No] <- gsub("([[:alpha:]]{3})(.*)", "\\1",DF$aa.end[row_No])
        DF$aa.end[row_No] <- gsub("^_759", "",DF$aa.end[row_No])
      } 
    }
  }
  
  assign("DF", DF, envir = .GlobalEnv)
}  

# Load relevant file: STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
DF_Full <- read.csv(file = "STAMP/20181113_syapse_export_all_variants_QC.csv",
                    header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")

DF_STAMP_Full <- DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)", ]

DF = DF_STAMP_Full
remove(DF_STAMP_Full, DF_Full)
DF$smpl.hgvsProtein <- gsub("[[:space:]]*$", "", DF$smpl.hgvsProtein)

# Structure DF columns
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
DF_var <- data.frame(patientid=DF$sys.uniqueId,
                     variant=DF$sys.label, gene=DF$base.gene,
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
  if (length(grep(FALSE, DF_var[i,6:15])) == 10) {
    DF_var$remain[i] <- TRUE
  } else {
    DF_var$remain[i] <- FALSE
  }}
remove(i)

# Subset based on matching of strings in smpl.hgvsCoding or smpl.hgvsProtein
# Classification order: SNV > Frameshift > indel > insertion > deletion > duplication
#----------------------------------------------
colname_sub <- c("sys.uniqueId", "sys.label", "base.gene", "smpl.hgvsCoding", "smpl.hgvsProtein")

# "p.=" indicates protein has not been analysed, RNA was, but no change is expected
DF_synonymous <- DF_var[DF_var$synonymous == TRUE, 1:5]
DF_synonymous$var.type <- "Synonymous"
colnames(DF_synonymous)[1:5] <- colname_sub
DF_synonymous$smpl.hgvsCoding <- as.character(DF_synonymous$smpl.hgvsCoding)
DF_synonymous$smpl.hgvsProtein <- as.character(DF_synonymous$smpl.hgvsProtein)
DF_synonymous_FULL <- inner_join(DF[DF$sys.label %in% DF_synonymous$sys.label, ], DF_synonymous, 
                                 by = colname_sub)

DF_upstream <- DF_var[DF_var$upstream == TRUE, 1:5]
DF_upstream$var.type <- NA
for (row_No in 1:nrow(DF_upstream)) {
  if (isTRUE(grepl("^c.(-)[[:digit:]]+[ATCG]>[ATCG][[:blank:]]*$", DF_upstream$coding[row_No]))) {
    DF_upstream$var.type[row_No] <- "SNV"
  }
}
colnames(DF_upstream)[1:5] <- colname_sub
DF_upstream$smpl.hgvsCoding <- as.character(DF_upstream$smpl.hgvsCoding)
DF_upstream$smpl.hgvsProtein <- as.character(DF_upstream$smpl.hgvsProtein)
DF_upstream_FULL <- inner_join(DF[DF$sys.label %in% DF_upstream$sys.label, ], DF_upstream, 
                                 by = colname_sub)
DF_upstream_FULL$var.position <- NA
DF_upstream_FULL$aa.start <- NA
DF_upstream_FULL$aa.end <- NA

DF_intron <- DF_var[DF_var$intron == TRUE, 1:5]
DF_intron$var.type <- NA
for (row_No in 1:nrow(DF_intron)) {
  if (isTRUE(grepl("del.*ins", DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "Indel"
  } else if (isTRUE(grepl("del", DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "Deletion"
  } else if (isTRUE(grepl("dup", DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "Duplication"
  } else if (isTRUE(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+[ATCG]>[ATCG][[:blank:]]*$", DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "SNV"
  }
}
colnames(DF_intron)[1:5] <- colname_sub
DF_intron$smpl.hgvsCoding <- as.character(DF_intron$smpl.hgvsCoding)
DF_intron$smpl.hgvsProtein <- as.character(DF_intron$smpl.hgvsProtein)
DF_intron_FULL <- inner_join(DF[DF$sys.label %in% DF_intron$sys.label, ], DF_intron, 
                               by = colname_sub)
DF_intron_FULL$var.position <- NA
DF_intron_FULL$aa.start <- NA
DF_intron_FULL$aa.end <- NA

DF_SNV <- DF_var[DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE &
                   DF_var$synonymous == FALSE, 1:5]
DF_SNV$var.type <- "SNV"

DF_Frameshift <- DF_var[DF_var$fs == TRUE &
                          DF_var$synonymous == FALSE & DF_var$upstream == FALSE &
                          DF_var$intron == FALSE, 1:5]
DF_Frameshift$var.type <- "Frameshift"

DF_indel <- DF_var[DF_var$indel == TRUE &
                     DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                     DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE &
                     DF_var$fs == FALSE , 1:5]
DF_indel$var.type <- "Indel"

DF_ins <- DF_var[DF_var$ins == TRUE &
                   DF_var$indel == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_ins$var.type <- "Insertion"

DF_del <- DF_var[DF_var$del == TRUE &
                   DF_var$ins == FALSE & DF_var$indel == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_del$var.type <- "Deletion"

DF_dup <- DF_var[DF_var$dup == TRUE &
                   DF_var$del == FALSE & DF_var$ins == FALSE & DF_var$indel == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$SNVprotein == FALSE & DF_var$SNVcoding == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_dup$var.type <- "Duplication"

DF_remain <- DF_var[DF_var$remain == TRUE, 1:5]
remove(DF_var, row_No)

print(paste("All cases have been categorized based on type of nucleotide change, as indicated:",
            nrow(DF) == (nrow(DF_synonymous) + nrow(DF_upstream) + nrow(DF_intron) + nrow(DF_SNV) + 
                           nrow(DF_Frameshift) + nrow(DF_indel) + nrow(DF_ins) + nrow(DF_del) + 
                           nrow(DF_dup) + nrow(DF_remain))),sep="")

# Unique genes per variant category
#----------------------------------------------
gene.no <- function(DF, gene_type) {
  print(paste(gene_type, " mutations: n=",
              length(unique(DF[[3]])), " genes and n=", length(DF[[3]]), " cases)", sep=""))
}

# gene.no(DF = DF_synonymous, gene_type = "Synonymous")
# gene.no(DF = DF_upstream, gene_type = "Upstream")
# gene.no(DF = DF_intron, gene_type = "Intronic")
# gene.no(DF = DF_SNV, gene_type = "SNV")
# gene.no(DF = DF_Frameshift, gene_type = "Frameshift")
# gene.no(DF = DF_indel, gene_type = "Indels")
# gene.no(DF = DF_ins, gene_type = "Insertions")
# gene.no(DF = DF_del, gene_type = "Deletions")
# gene.no(DF = DF_dup, gene_type = "Duplications")
# gene.no(DF = DF_remain, gene_type = "Unknown")

print("Synonymous, upstream & intronic mutations are not mapped to lollipop plots")
remove(DF_synonymous,DF_upstream, DF_intron,DF_remain)

DF_Map <- rbind(DF_SNV, DF_Frameshift,DF_indel,DF_ins,DF_del,DF_dup)
colnames(DF_Map)[1:5] <- colname_sub
DF_Map$smpl.hgvsCoding <- as.character(DF_Map$smpl.hgvsCoding)
DF_Map$smpl.hgvsProtein <- as.character(DF_Map$smpl.hgvsProtein)
DF_Map_All <- inner_join(DF[DF$sys.label %in% DF_Map$sys.label,], DF_Map, by = colname_sub)

# Subset cases missing smpl.hgvsProtein field
#----------------------------------------------
DF_NAprotein <- DF_Map_All[is.na(DF_Map_All$smpl.hgvsProtein),]
DF_NAprotein$var.position <- NA
DF_NAprotein$aa.start <- NA
DF_NAprotein$aa.end <- NA

# Extract variant positions
#----------------------------------------------
DF_STAMP_4Map <- DF_Map_All[!is.na(DF_Map_All$smpl.hgvsProtein),]
Extract_VarPosition(DF = DF_STAMP_4Map)
DF_STAMP_4Map <- DF

remove(DF_SNV,DF_Frameshift,DF_indel,DF_ins,DF_del,DF_dup,DF_Map,DF_Map_All,DF)

# ## Simple review of data
# #----------------------------------------------
# print(paste("Total number of unique patient_id: ", length(unique(DF_STAMP_4Map$sys.uniqueId)), sep=""))
# print(paste("Total number of unique specimen sites: ", length(sort(unique(DF_STAMP_4Map$smpl.specimenSite))), sep=""))
# sort(unique(DF_STAMP_4Map$smpl.specimenSite))
# print(paste("Total number of unique Dx: ", length(sort(unique(DF_STAMP_4Map$smpl.ppDiagnosticSummary))), sep=""))
# sort(unique(DF_STAMP_4Map$smpl.ppDiagnosticSummary))
# print(paste("Total number of unique genes: ", length(sort(unique(DF_STAMP_4Map$base.gene))), sep=""))
# sort(unique(DF_STAMP_4Map$base.gene))
# print(paste("Total number of unique pathogenicity statuses: ", length(sort(unique(DF_STAMP_4Map$smpl.pathogenicityStatus))), sep=""))
# table(sort(DF_STAMP_4Map$smpl.pathogenicityStatus))
# print(paste("Total number of unique variant categories: ", length(sort(unique(DF_STAMP_4Map$var.type))), sep=""))
# table(sort(DF_STAMP_4Map$var.type))
# print(paste("Total number of unique start codons: ", length(sort(unique(DF_STAMP_4Map$aa.start))), sep=""))
# table(sort(DF_STAMP_4Map$aa.start))
# print(paste("Total number of unique end codons: ", length(sort(unique(DF_STAMP_4Map$aa.end))), sep=""))
# table(sort(DF_STAMP_4Map$aa.end))

## Write to local computer
#----------------------------------------------
write.csv(DF_STAMP_4Map, file = "STAMP/Mutation_Hotspot/20181114_syapse_export_DF_STAMP_4Map.csv",
          na = "NA", row.names = FALSE)
write.csv(DF_NAprotein[,1:64], file = "STAMP/Mutation_Hotspot/20181114_syapse_export_DF_NAprotein.csv",
          na = "NA", row.names = FALSE)

## Subset 20181113_syapse_export_all_variants_QC.csv based on mutation type = 9918 total STAMP entries
## Classification of mutations into "AMPLIFICATION", "DELETION","FUSION","MUTATION"
## Output: "ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv"

# Extract variant positions
#----------------------------------------------
DF_STAMP_VariantAnno <- rbind(DF_STAMP_4Map[,1:64], DF_synonymous_FULL)
Extract_VarPosition(DF = DF_STAMP_VariantAnno)
DF_STAMP_VariantAnno <- DF

DF_STAMP_VariantAnno <- rbind(DF_STAMP_VariantAnno, DF_NAprotein,
                              DF_intron_FULL, DF_upstream_FULL)
remove(DF_synonymous_FULL,DF_STAMP_4Map,DF,DF_NAprotein,colname_sub,
       DF_intron_FULL, DF_upstream_FULL)

# Classification of variants for clinical trial matching
#----------------------------------------------
# Frameshift are a special type of amino acid deletion/insertion
# Source: http://varnomen.hgvs.org/recommendations/protein/variant/frameshift/
DF_STAMP_VariantAnno$var.anno <- NA
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  if (DF_STAMP_VariantAnno$var.type[row_No] == "SNV" | DF_STAMP_VariantAnno$var.type[row_No] == "Indel" |
      DF_STAMP_VariantAnno$var.type[row_No] == "Frameshift") {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "MUTATION"
  } else if (DF_STAMP_VariantAnno$var.type[row_No] == "Duplication" | 
             DF_STAMP_VariantAnno$var.type[row_No] == "Insertion") {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "AMPLIFICATION"
  } else if (DF_STAMP_VariantAnno$var.type[row_No] == "Deletion") {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "DELETION"
  } else {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "OTHER"
  }
}
remove(gene.no, Extract_VarPosition,row_No)

# Remove patients without DOB = 4 entries
DF_STAMP_VariantAnno <- DF_STAMP_VariantAnno[!is.na(DF_STAMP_VariantAnno$base.dob),]
ncol_STAMP <- as.numeric(ncol(DF_STAMP_VariantAnno))

## Structure patient DOB and input current age
#----------------------------------------------
curr_year <- as.numeric(gsub("^([[:digit:]]{2})([[:digit:]]{2})", "\\2", 
                             format(as.Date(Sys.Date(), format="%d/%m/%Y"),"%Y")))

DF_STAMP_VariantAnno$month = gsub("(^[[:digit:]]{,2})(.*)", "\\1", DF_STAMP_VariantAnno$base.dob)
DF_STAMP_VariantAnno$day=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})(.*)", "\\3", DF_STAMP_VariantAnno$base.dob)
DF_STAMP_VariantAnno$year=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})([/])([[:digit:]]{2})[[:blank:]]*$", "\\5", 
                                DF_STAMP_VariantAnno$base.dob)

# Assume no individual is >= 100yo
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  if (DF_STAMP_VariantAnno$year[row_No] > curr_year) {
    DF_STAMP_VariantAnno$year[row_No] <- paste("19", DF_STAMP_VariantAnno$year[row_No], sep="")
  } else if (DF_STAMP_VariantAnno$year[row_No] <= curr_year) {
    DF_STAMP_VariantAnno$year[row_No] <- paste("20", DF_STAMP_VariantAnno$year[row_No], sep="")
  }
}

DF_STAMP_VariantAnno$patient.dob <- NA
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  DF_STAMP_VariantAnno$patient.dob <- paste(DF_STAMP_VariantAnno$month, 
                                             DF_STAMP_VariantAnno$day, 
                                             DF_STAMP_VariantAnno$year, sep="/")
}

DF_STAMP_VariantAnno$patient.dob <- as.Date(DF_STAMP_VariantAnno$patient.dob, "%m/%d/%Y")
DF_STAMP_VariantAnno <- DF_STAMP_VariantAnno[,c(1:ncol_STAMP,ncol(DF_STAMP_VariantAnno))]

# Age rounded down to nearest integer
DF_STAMP_VariantAnno$current.age <- as.numeric(floor(age_calc(dob = DF_STAMP_VariantAnno$patient.dob, 
                                                               enddate = Sys.Date(), units = "years")))

# Plot age distribution
DF <- data.frame(patient.id = DF_STAMP_VariantAnno$sys.uniqueId, 
                 age = DF_STAMP_VariantAnno$current.age,
                 assayName = DF_STAMP_VariantAnno$smpl.assayName)
patient.list <- as.character(sort(unique(DF$patient.id)))
DF_new <- data.frame()
for (id_num in 1:length(patient.list)) {
  DF_new <- rbind(DF_new, DF[min(which(DF$patient.id == patient.list[id_num])),])
}
DF <- DF_new

tiff(filename = "STAMP/AgeDistribution.tiff",
     width = 15, height = 7, units = "in", res = 200)

ggplot(DF, aes(DF$age, fill=DF$assayName)) +
  geom_histogram(bins = 97, col="gray") +
  labs(subtitle = paste("N = ", nrow(DF), sep="")) +
  
  scale_x_continuous(name="Current Age", breaks = seq(0, 100, 5)) +
  scale_y_continuous(name="Number of Individuals", breaks = seq(0,150,10)) + 
  
  scale_fill_discrete(name = "NGS Assay") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust=1, face="bold",size=15),
        legend.position = "bottom",
        legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
        legend.text=element_text(size=10),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

dev.off()
remove(DF_new,DF,curr_year, row_No,ncol_STAMP,id_num,patient.list)

# ## Simple review of data
# #----------------------------------------------
# print(paste("Total number of variant categories: ", length(sort(unique(DF_STAMP_VariantAnno$var.type))), sep=""))
# table(sort(DF_STAMP_VariantAnno$var.type))
# print(paste("Total number of variant annotations for clinical trials: ", length(sort(unique(DF_STAMP_VariantAnno$var.anno))), sep=""))
# table(sort(DF_STAMP_VariantAnno$var.anno))
# print("Variants with an annotation == AMPLIFICATION have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "AMPLIFICATION"]))
# print("Variants with an annotation == DELETION have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "DELETION"]))
# print("Variants with an annotation == MUTATION have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "MUTATION"]))
# print("Variants with an annotation == OTHER have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "OTHER"]))

## Write to local computer
#----------------------------------------------
write.csv(DF_STAMP_VariantAnno, file = "ClinicalTrialMatching/20181114_syapse_export_DF_STAMP_VariantAnno.csv",
          na = "NA", row.names = FALSE)
