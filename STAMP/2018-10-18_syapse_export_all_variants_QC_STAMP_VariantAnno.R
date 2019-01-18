setwd("~/Documents/ClinicalDataScience_Fellowship/")

# FUNCTION
#----------------------------------------------
# Extract variant positions 
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
    } else if (DF$var.type[row_No] == "Duplication" | grepl("Frameshift", DF$var.type[row_No])) {
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
    } else if (grepl("Frameshift", DF$var.type[row_No])) {
      DF$aa.end[row_No] <- gsub("([[:alpha:]]{3})(.*)", "\\1", DF$aa.end[row_No])
    } else if (DF$var.type[row_No] == "Delins" | DF$var.type[row_No] == "Insertion" | 
               DF$var.type[row_No] == "Deletion" |DF$var.type[row_No] == "Duplication") {
      DF$aa.end[row_No] <- gsub("(^_)([[:alpha:]]{3})(.*)", "\\2",DF$aa.end[row_No])
      
      if (DF$var.type[row_No] == "Delins") {
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

# Unique genes per variant category
gene.no <- function(DF, gene_type) {
  print(paste(gene_type, " mutations: n=",
              length(DF[[3]]), " entries and n=", length(unique(DF[[3]])), " genes)", sep=""))
}

# Load relevant file: STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
DF_Full <- read.csv(file = paste("STAMP/", Syapse_Export_timestamp, "_syapse_export_all_variants_QC.csv",sep=""),
                    header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")

DF_STAMP_Full <- DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" |
                           DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)", ]

DF = DF_STAMP_Full
print(paste("Total STAMP entries: n=",(nrow(DF)), sep=""))
print(paste("Total number of unique patient_id: ", length(unique(DF$sys.uniqueId)), sep=""))

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
                     delins=grepl("del.*ins", DF$sys.label),
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

# Classification correction
DF_var$SNVprotein[which(DF_var$SNVprotein == TRUE & 
                          (DF_var$delins == TRUE | DF_var$ins == TRUE | 
                             DF_var$del == TRUE | DF_var$dup == TRUE))] <- FALSE

# Subset based on matching of strings in smpl.hgvsCoding or smpl.hgvsProtein
# Classification order: SNV > Frameshift > delins > insertion > deletion > duplication
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
# sort(unique(DF_synonymous$base.gene))

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
# sort(unique(DF_upstream$base.gene))

DF_intron <- DF_var[DF_var$intron == TRUE, ]
DF_intron$var.type <- NA
for (row_No in 1:nrow(DF_intron)) {
  if (isTRUE(DF_intron$fs[row_No])) {
    DF_intron$var.type[row_No] <- "Frameshift"
  } else if (isTRUE(DF_intron$delins[row_No])) {
    DF_intron$var.type[row_No] <- "Delins"
  } else if (isTRUE(DF_intron$del[row_No])) {
    DF_intron$var.type[row_No] <- "Deletion"
  } else if (isTRUE(DF_intron$dup[row_No])) {
    DF_intron$var.type[row_No] <- "Duplication"
  } else if (isTRUE(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+[ATCG]>[ATCG][[:blank:]]*$", 
                          DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "SNV"
  }
}
DF_intron <- DF_intron[,c(1:5,17)]
colnames(DF_intron)[1:5] <- colname_sub
DF_intron$smpl.hgvsCoding <- as.character(DF_intron$smpl.hgvsCoding)
DF_intron$smpl.hgvsProtein <- as.character(DF_intron$smpl.hgvsProtein)
DF_intron_FULL <- inner_join(DF[DF$sys.label %in% DF_intron$sys.label, ], DF_intron, 
                             by = colname_sub)
DF_intron_FULL$var.position <- NA
DF_intron_FULL$aa.start <- NA
DF_intron_FULL$aa.end <- NA
# sort(unique(DF_intron$base.gene))

DF_SNV <- DF_var[DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE &
                   DF_var$synonymous == FALSE, 1:5]
DF_SNV$var.type <- "SNV"

DF_Frameshift <- DF_var[DF_var$fs == TRUE &
                          DF_var$synonymous == FALSE & DF_var$upstream == FALSE &
                          DF_var$intron == FALSE,]
DF_Frameshift$var.type <- "Frameshift"
for (row_No in 1:nrow(DF_Frameshift)) {
  if (isTRUE(DF_Frameshift$delins[row_No])) {
    DF_Frameshift$var.type[row_No] <- "Frameshift_Delins"
  } else if (isTRUE(DF_Frameshift$del[row_No])) {
    DF_Frameshift$var.type[row_No] <- "Frameshift_Deletion"
  } else if (isTRUE(DF_Frameshift$ins[row_No])) {
    DF_Frameshift$var.type[row_No] <- "Frameshift_Insertion"
  } else if (isTRUE(DF_Frameshift$dup[row_No])) {
    DF_Frameshift$var.type[row_No] <- "Frameshift_Duplication"
  }
}
DF_Frameshift <- DF_Frameshift[,c(1:5,17)]

DF_delins <- DF_var[DF_var$delins == TRUE &
                      DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                      DF_var$fs == FALSE , 1:5]
DF_delins$var.type <- "Delins"

DF_ins <- DF_var[DF_var$ins == TRUE &
                   DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_ins$var.type <- "Insertion"

DF_del <- DF_var[DF_var$del == TRUE &
                   DF_var$ins == FALSE & DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_del$var.type <- "Deletion"

DF_dup <- DF_var[DF_var$dup == TRUE &
                   DF_var$del == FALSE & DF_var$ins == FALSE & DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_dup$var.type <- "Duplication"

DF_remain <- DF_var[DF_var$remain == TRUE, 1:5]
remove(DF_var, row_No)

## Data Review
#----------------------------------------------
print(paste("All cases have been categorized based on type of nucleotide change, as indicated:",
            nrow(DF) == (nrow(DF_synonymous) + nrow(DF_upstream) + nrow(DF_intron) + nrow(DF_SNV) + 
                           nrow(DF_Frameshift) + nrow(DF_delins) + nrow(DF_ins) + nrow(DF_del) + 
                           nrow(DF_dup) + nrow(DF_remain))),sep="")

cat("\n")

gene.no(DF = DF_synonymous, gene_type = "Synonymous")
cat("\n")

gene.no(DF = DF_upstream, gene_type = "Upstream")
Temp_FS <- DF_upstream[DF_upstream$var.type == "SNV",]
print(paste("Upstream_SNV mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
cat("\n")

gene.no(DF = DF_intron, gene_type = "Intronic")
Temp_FS <- DF_intron[DF_intron$var.type == "SNV",]
print(paste("Intronic_SNV mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
Temp_FS <- DF_intron[DF_intron$var.type == "Delins",]
print(paste("Intronic_Delins mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
Temp_FS <- DF_intron[DF_intron$var.type == "Deletion",]
print(paste("Intronic_Deletion mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
Temp_FS <- DF_intron[DF_intron$var.type == "Duplication",]
print(paste("Intronic_Duplication mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
Temp_FS <- DF_intron[DF_intron$var.type == "Frameshift",]
print(paste("Intronic_Frameshift mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$base.gene)), " genes)", sep=""))
cat("\n")

gene.no(DF = DF_SNV, gene_type = "SNV")
cat("\n")

Temp_InFrame <- rbind(DF_delins, DF_ins, DF_del, DF_dup)
print(paste("In-Frame mutations: n=", nrow(Temp_InFrame), " entries and n=",
            length(unique(Temp_InFrame$gene)), " genes)", sep=""))
gene.no(DF = DF_delins, gene_type = "Delins")
gene.no(DF = DF_ins, gene_type = "Insertions")
gene.no(DF = DF_del, gene_type = "Deletions")
gene.no(DF = DF_dup, gene_type = "Duplications")
cat("\n")

gene.no(DF = DF_Frameshift, gene_type = "Frameshift")
Temp_FS <- DF_Frameshift[DF_Frameshift$var.type == "Frameshift_Delins",]
print(paste("Frameshift_Delins mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$gene)), " genes)", sep=""))
Temp_FS <- DF_Frameshift[DF_Frameshift$var.type == "Frameshift_Insertion",]
print(paste("Frameshift_Insertion mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$gene)), " genes)", sep=""))
Temp_FS <- DF_Frameshift[DF_Frameshift$var.type == "Frameshift_Deletion",]
print(paste("Frameshift_Deletion mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$gene)), " genes)", sep=""))
Temp_FS <- DF_Frameshift[DF_Frameshift$var.type == "Frameshift_Duplication",]
print(paste("Frameshift_Duplication mutations: n=", nrow(Temp_FS), " entries and n=",
            length(unique(Temp_FS$gene)), " genes)", sep=""))
remove(Temp_FS,Temp_InFrame)
cat("\n")

print(paste("Synonymous, upstream & intronic mutations are not mapped to lollipop plots (n=",
            (nrow(DF_synonymous) + nrow(DF_upstream) + nrow(DF_intron)), " entries)", sep=""))
remove(DF_synonymous,DF_upstream, DF_intron,DF_remain)

# Merge DFs of interest
DF_Map <- rbind(DF_SNV, DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup)
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
print(paste("STAMP entries missing smpl.hgvsProtein: n=",(nrow(DF_NAprotein)), sep=""))

# Extract variant positions
#----------------------------------------------
DF_STAMP_4Map <- DF_Map_All[!is.na(DF_Map_All$smpl.hgvsProtein),]
Extract_VarPosition(DF = DF_STAMP_4Map)
DF_STAMP_4Map <- DF
print(paste("STAMP entries mapped to lollipop plots: n=",(nrow(DF_STAMP_4Map)), sep=""))

remove(DF_SNV,DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup,DF_Map,DF_Map_All,DF)

## Data Review
#----------------------------------------------
cat("\n")
print("General overview for DF_STAMP_4Map entries")
print(paste("Total number of unique patient_id: ", length(unique(DF_STAMP_4Map$sys.uniqueId)), sep=""))
print(paste("Total number of unique specimen sites: ", length(sort(unique(DF_STAMP_4Map$smpl.specimenSite))), sep=""))
sort(unique(DF_STAMP_4Map$smpl.specimenSite))
print(paste("Total number of unique Dx: ", length(sort(unique(DF_STAMP_4Map$smpl.ppDiagnosticSummary))), sep=""))
sort(unique(DF_STAMP_4Map$smpl.ppDiagnosticSummary))
print(paste("Total number of unique genes: ", length(sort(unique(DF_STAMP_4Map$base.gene))), sep=""))
sort(unique(DF_STAMP_4Map$base.gene))
print(paste("Total number of unique pathogenicity statuses: ", length(sort(unique(DF_STAMP_4Map$smpl.pathogenicityStatus))), sep=""))
print(table(sort(DF_STAMP_4Map$smpl.pathogenicityStatus)))
print(paste("Total number of unique variant categories: ", length(sort(unique(DF_STAMP_4Map$var.type))), sep=""))
print(table(sort(DF_STAMP_4Map$var.type)))
print(paste("Total number of unique start codons: ", length(sort(unique(DF_STAMP_4Map$aa.start))), sep=""))
print(table(sort(DF_STAMP_4Map$aa.start)))
print(paste("Total number of unique end codons: ", length(sort(unique(DF_STAMP_4Map$aa.end))), sep=""))
print(table(sort(DF_STAMP_4Map$aa.end)))
cat("\n")

## Write to local computer
#----------------------------------------------
write.csv(DF_STAMP_4Map, file = paste("Mutation_Hotspot/", Syapse_Export_timestamp, "_syapse_export_DF_STAMP_4Map.csv", sep=""),
          na = "NA", row.names = FALSE)
write.csv(DF_NAprotein[,1:64], file = paste("Mutation_Hotspot/", Syapse_Export_timestamp, "_syapse_export_DF_NAprotein.csv", sep=""),
          na = "NA", row.names = FALSE)

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
  if (DF_STAMP_VariantAnno$var.type[row_No] == "Synonymous") {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "OTHER"
  } else {
    DF_STAMP_VariantAnno$var.anno[row_No] <- "MUTATION"
  }
}

# Remove patients without DOB = 4 entries
DF_STAMP_VariantAnno <- DF_STAMP_VariantAnno[!is.na(DF_STAMP_VariantAnno$base.dob),]
ncol_STAMP <- as.numeric(ncol(DF_STAMP_VariantAnno))

print(paste("STAMP entries with valid DOB: n=",(nrow(DF_STAMP_VariantAnno)), sep=""))
print(paste("Unique patient ID with valid DOB: n=",(length(unique(DF_STAMP_VariantAnno$sys.uniqueId))), sep=""))

## Structure patient DOB and input current age
#----------------------------------------------
curr_year <- as.numeric(gsub("^([[:digit:]]{2})([[:digit:]]{2})", "\\2", 
                             format(as.Date(Sys.Date(), format="%Y-%m-%d"),"%Y")))

DF_STAMP_VariantAnno$month = gsub("(^[[:digit:]]{,2})(.*)", "\\1", DF_STAMP_VariantAnno$base.dob)
DF_STAMP_VariantAnno$day=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})(.*)", "\\3", DF_STAMP_VariantAnno$base.dob)
DF_STAMP_VariantAnno$year=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})([/])([[:digit:]]{2})[[:blank:]]*$", "\\5", 
                               DF_STAMP_VariantAnno$base.dob)

# Assume no individual is >= 100yo
# sort(as.numeric(unique(DF_STAMP_VariantAnno$year)))
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  if (DF_STAMP_VariantAnno$year[row_No] > curr_year) {
    DF_STAMP_VariantAnno$year[row_No] <- paste("19", DF_STAMP_VariantAnno$year[row_No], sep="")
  } else if (DF_STAMP_VariantAnno$year[row_No] <= curr_year) {
    DF_STAMP_VariantAnno$year[row_No] <- paste("20", DF_STAMP_VariantAnno$year[row_No], sep="")
  }
}
# sort(as.numeric(unique(DF_STAMP_VariantAnno$year)))

DF_STAMP_VariantAnno$patient.dob <- NA
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  DF_STAMP_VariantAnno$patient.dob <- paste(DF_STAMP_VariantAnno$month, 
                                            DF_STAMP_VariantAnno$day, 
                                            DF_STAMP_VariantAnno$year, sep="/")
}

DF_STAMP_VariantAnno$patient.dob <- as.Date(DF_STAMP_VariantAnno$patient.dob, "%m/%d/%Y")
DF_STAMP_VariantAnno <- DF_STAMP_VariantAnno[,c(1:ncol_STAMP,ncol(DF_STAMP_VariantAnno))]

# Age rounded down to nearest integer -- relative to smpl.dateReceived
DF_STAMP_VariantAnno$dateReceived <- gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                          DF_STAMP_VariantAnno$smpl.dateReceived)

# Use sys.date_created as alternative for missing smpl.dateReceived
row_NA <- which(is.na(DF_STAMP_VariantAnno$dateReceived))
if (length(row_NA) > 0) {
  DF_STAMP_VariantAnno$dateReceived[row_NA] <- gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                                    DF_STAMP_VariantAnno$sys.date_created[row_NA])
}

DF_STAMP_VariantAnno$dateReceived <- as.Date(DF_STAMP_VariantAnno$dateReceived, format = "%Y-%m-%d")

DF_STAMP_VariantAnno$current.age <- as.numeric(floor(age_calc(dob = DF_STAMP_VariantAnno$patient.dob, 
                                                              enddate = DF_STAMP_VariantAnno$dateReceived,
                                                              units = "years")))

# Remove dateReceived column
DF_STAMP_VariantAnno <- DF_STAMP_VariantAnno[, !(colnames(DF_STAMP_VariantAnno) %in% c("dateReceived"))]

## Data Review
#----------------------------------------------
Temp_Age <- unique(DF_STAMP_VariantAnno$sys.uniqueId[DF_STAMP_VariantAnno$current.age < 18])
print(paste("No. child patients (birth-17yo): n=",(length(Temp_Age)), sep=""))
Temp_Age <- unique(DF_STAMP_VariantAnno$sys.uniqueId[DF_STAMP_VariantAnno$current.age > 65])
print(paste("No. Older Adult patients (65yo+): n=",(length(Temp_Age)), sep=""))
Temp_Age <- DF_STAMP_VariantAnno[DF_STAMP_VariantAnno$current.age >= 18 &
                                   DF_STAMP_VariantAnno$current.age < 65,]
print(paste("No. Adult patients (18-64yo): n=",(length(unique(Temp_Age$sys.uniqueId))), sep=""))
cat("\n")
print("General overview for adult patient entries")
print(paste("Total number of STAMP entries: ", nrow(Temp_Age), sep=""))
print(paste("Total number of unique specimen sites: ", length(sort(unique(Temp_Age$smpl.specimenSite))), sep=""))
print(paste("Total number of unique Dx: ", length(sort(unique(Temp_Age$smpl.ppDiagnosticSummary))), sep=""))
print(paste("Total number of unique genes: ", length(sort(unique(Temp_Age$base.gene))), sep=""))
print(paste("Total number of unique pathogenicity statuses: ", length(sort(unique(Temp_Age$smpl.pathogenicityStatus))), sep=""))
print(table(sort(Temp_Age$smpl.pathogenicityStatus)))
print(paste("Total number of unique variant categories: ", length(sort(unique(Temp_Age$var.type))), sep=""))
print(table(sort(Temp_Age$var.type)))
cat("\n")
remove(Temp_Age)

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

tiff(filename = paste("STAMP/Syapse_", Syapse_Export_timestamp, "_AgeDistribution.tiff", sep=""),
     width = 15, height = 7, units = "in", res = 200)

plot <- ggplot(DF, aes(DF$age, fill=DF$assayName)) +
  geom_histogram(bins = (max(DF$age) - min(DF$age) +1), col="gray") +
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

print(plot)
dev.off()
remove(DF_new,DF,curr_year, row_No,ncol_STAMP,id_num,patient.list)

## Data Review
#----------------------------------------------
# print(paste("Total number of variant categories: ", length(sort(unique(DF_STAMP_VariantAnno$var.type))), sep=""))
# print(paste("Total number of variant annotations for clinical trials: ", length(sort(unique(DF_STAMP_VariantAnno$var.anno))), sep=""))
# table(sort(DF_STAMP_VariantAnno$var.anno))
# print("Variants with an annotation == MUTATION have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "MUTATION"]))
# print("Variants with an annotation == OTHER have the following categories")
# table(sort(DF_STAMP_VariantAnno$var.type[DF_STAMP_VariantAnno$var.anno == "OTHER"]))

## Merge with primary tumor site data
#----------------------------------------------
DF_TumorSite <- read.csv(file = "STAMP/2018-12-17_syapse_primary_tumor_site.csv",
                         header = TRUE, na.strings = c("NA","None"), stringsAsFactors = FALSE,sep = ",")
# Remove extraneous columns
DF_TumorSite <- DF_TumorSite[,c(1,2,7)]
# Shorten column names of DF
colnames(DF_TumorSite) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_TumorSite))

# Confirm smpl.primaryTumorSite and smpl.histologicalDiagnosis are same for all entries of each unique patient
patient.list <- sort(unique(DF_TumorSite$sys.uniqueId))
for (row_No in 1:length(patient.list)) {
  DF <- DF_TumorSite[DF_TumorSite$sys.uniqueId == patient.list[row_No],c(1:3)]
  DF <- DF[complete.cases(DF$smpl.primaryTumorSite),]
  if (nrow(DF) > 1) {
    for (rep_No in 2:nrow(DF)) {
      if (isTRUE(DF$smpl.primaryTumorSite[1] == DF$smpl.primaryTumorSite[rep_No]) &
          isTRUE(DF$smpl.histologicalDiagnosis[1] == DF$smpl.histologicalDiagnosis[rep_No])) {
        DF[rep_No,c(1:2)] <- NA
      } else if (isTRUE(DF$smpl.primaryTumorSite[1] == DF$smpl.primaryTumorSite[rep_No]) &
                 isTRUE(is.na(DF$smpl.histologicalDiagnosis[1]) & is.na(DF$smpl.histologicalDiagnosis[rep_No]))) {
        DF[rep_No,c(1:2)] <- NA
      }
    }
    DF <- DF[complete.cases(DF$smpl.primaryTumorSite),]
    if (nrow(DF) > 1) {
      print(DF)
    }
  }
}

DF_TumorSite <- DF_TumorSite %>% dplyr::distinct(sys.uniqueId, .keep_all = TRUE)
DF_STAMP_VariantAnno$sys.uniqueId <- as.character(DF_STAMP_VariantAnno$sys.uniqueId)
DF_STAMP_VariantAnno <- left_join(DF_STAMP_VariantAnno, DF_TumorSite,
                                  by = c("sys.uniqueId"))

# primaryTumorSite grouped according to medical specializations 
#----------------------------------------------
DF_STAMP_VariantAnno$smpl.primaryTumorSite <- tolower(DF_STAMP_VariantAnno$smpl.primaryTumorSite)

# Confirm smpl.primaryTumorSite is classified 
primaryTumorSite.STAMP <- sort(unique(DF_STAMP_VariantAnno$smpl.primaryTumorSite))
primaryTumorSite.key <- sort(unique(DiseaseGroupCategory_LongFormat$primaryTumorSite))

for (elem_No in 1:length(primaryTumorSite.STAMP)) {
  if (isTRUE(is.element(primaryTumorSite.STAMP[elem_No], primaryTumorSite.key) == FALSE)) {
    cat("\n")
    print(paste(primaryTumorSite.STAMP[elem_No], 
                ": primaryTumorSite has not been classified into DiseaseGroupCategory", sep=""))
    cat("\n")
  }
}
remove(primaryTumorSite.STAMP,primaryTumorSite.key)

DF_STAMP_VariantAnno$primaryTumorSite.category <- NA
for (row_No in 1:nrow(DF_STAMP_VariantAnno)) {
  DiseaseGroupCategory.name <- 
    DiseaseGroupCategory_LongFormat$Disease.Group.category[which(DiseaseGroupCategory_LongFormat$primaryTumorSite == 
                                                                   DF_STAMP_VariantAnno$smpl.primaryTumorSite[row_No])] 
  
  DF_STAMP_VariantAnno$primaryTumorSite.category[row_No] <- paste(as.character(DiseaseGroupCategory.name), collapse=", ")
}

# If missing smpl.primaryTumorSite, primaryTumorSite.category = "Unknown"
DF_STAMP_VariantAnno$primaryTumorSite.category[which(DF_STAMP_VariantAnno$primaryTumorSite.category == "")] <- "unknown"

print(paste("Total number of unique primaryTumorSite.category: ", length(sort(unique(DF_STAMP_VariantAnno$primaryTumorSite.category))), sep=""))
print(table(DF_STAMP_VariantAnno$primaryTumorSite.category))
# breast          dermatology dermatology, sarcoma        endocrinology     gastroenterology 
# 189                  593                   11                    7                 1100 
# genitourinary           gynecology       hematolymphoid            neurology        ophthalmology 
# 148                   94                  112                  227                   23 
# otolaryngology          pulmonology              sarcoma              unknown 

remove(DF,DF_TumorSite,patient.list,row_No,rep_No,Extract_VarPosition,gene.no,plot,
       elem_No,row_NA,DiseaseGroupCategory.name)

## Data Review
#----------------------------------------------
Temp_Age <- DF_STAMP_VariantAnno[DF_STAMP_VariantAnno$current.age >= 18 &
                                   DF_STAMP_VariantAnno$current.age < 65,]
Temp_col <- sapply(Temp_Age$smpl.primaryTumorSite, tolower)
sort(unique(Temp_col))
Temp_col <- sapply(Temp_Age$primaryTumorSite.category, tolower)
sort(unique(Temp_col))
Temp_col <- sapply(Temp_Age$smpl.specimenSite, tolower)
sort(unique(Temp_col))
Temp_col <- sapply(Temp_Age$smpl.ppDiagnosticSummary, tolower)
sort(unique(Temp_col))
remove(Temp_Age,Temp_col)

## Write to local computer
#----------------------------------------------
write.csv(DF_STAMP_VariantAnno, 
          file = paste("ClinicalTrialMatching/", Syapse_Export_timestamp, "_syapse_export_DF_STAMP_VariantAnno.csv", sep=""),
          na = "NA", row.names = FALSE)
