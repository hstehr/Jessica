## Classification #1 (var.type): Synonymous, Upstream, Intronic, SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Classification #2 (var.anno): "MUTATION","OTHER"
## Classification #3 (PrimaryTumorSite.Category)
## Classification #4 (Exon_Number)
## Output: "syapse_export_all_variants_QC.csv" [overwrite]

# Categorize protein sequence variants (var.type)
#----------------------------------------------
# Subset columns of interest for simplification
# Beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron
# End of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron

DF <- STAMP_DF
DF_var <- data.frame(patientid=DF$PatientID,
                     variant=DF$VariantLabel, gene=DF$VariantGene,
                     coding=DF$VariantHGVSCoding, protein=DF$VariantHGVSProtein,
                     
                     # Categorize based on HGVS.Protein
                     synonymous=grepl("=",DF$VariantHGVSProtein),
                     SNVprotein=grepl("^p.[[:upper:]]{1}[[:lower:]]{2}[[:digit:]]+[[:upper:]]{1}[[:lower:]]{2}",
                                      DF$VariantHGVSProtein),
                     fs=grepl("fs", DF$VariantHGVSProtein),
                     
                     # Categorize based on HGVS.Coding
                     intron=grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+.*", DF$VariantHGVSCoding),
                     upstream=grepl("^c.-[[:digit:]]+[ATCG]>[ATCG]", DF$VariantHGVSCoding),
                     SNVcoding=grepl("^c.[[:digit:]]+[ATCG]>[ATCG]", DF$VariantHGVSCoding),
                     ins=grepl("ins", DF$VariantHGVSCoding),
                     del=grepl("del", DF$VariantHGVSCoding),
                     dup=grepl("dup", DF$VariantHGVSCoding),
                     
                     # Categorize based on Variant.Label
                     delins=grepl("del.*ins", DF$VariantLabel),
                     
                     # Empty variable for unclassified variant
                     remain=NA,
                     stringsAsFactors = FALSE)

for (i in 1:nrow(DF_var)) {
  if (length(grep(FALSE, DF_var[i,6:15])) == 10) { DF_var$remain[i] <- TRUE
  } else { DF_var$remain[i] <- FALSE
  }
}

# Classification correction
DF_var$SNVprotein[which(DF_var$SNVprotein == TRUE & 
                          (DF_var$delins == TRUE | DF_var$ins == TRUE | 
                             DF_var$del == TRUE | DF_var$dup == TRUE))] <- FALSE

# Subset based on matching of strings in smpl.hgvsCoding or smpl.hgvsProtein
# Classification order: SNV > Frameshift > In-Frame (delins > insertion > deletion > duplication)
#----------------------------------------------
colname_subset <- c("PatientID", "VariantLabel", "VariantGene", "VariantHGVSCoding", "VariantHGVSProtein")

# Synonymous Variants
# "p.=" indicates protein has not been analysed, RNA was, but no change is expected
DF_synonymous <- DF_var[DF_var$synonymous == TRUE, 1:5]
DF_synonymous$var.type <- "Synonymous"

# Upstream Variants
DF_upstream <- DF_var[DF_var$upstream == TRUE, ]
DF_upstream$var.type <- NA
for (row_No in 1:nrow(DF_upstream)) {
  if (isTRUE(grepl("^c.(-)[[:digit:]]+[ATCG]>[ATCG]", DF_upstream$coding[row_No]))) {
    DF_upstream$var.type[row_No] <- "SNV"
  } else if (isTRUE(DF_upstream$fs[row_No])) { DF_upstream$var.type[row_No] <- "Frameshift"
  } else if (isTRUE(DF_upstream$delins[row_No])) { DF_upstream$var.type[row_No] <- "Delins"
  } else if (isTRUE(DF_upstream$ins[row_No])) { DF_upstream$var.type[row_No] <- "Insertion"
  } else if (isTRUE(DF_upstream$del[row_No])) { DF_upstream$var.type[row_No] <- "Deletion"
  } else if (isTRUE(DF_upstream$dup[row_No])) { DF_upstream$var.type[row_No] <- "Duplication"
  }
}
DF_upstream <- DF_upstream[,c(1:5,17)]
colnames(DF_upstream)[1:5] <- colname_subset
DF_upstream <- inner_join(DF[DF$VariantLabel %in% DF_upstream$VariantLabel, ], 
                          DF_upstream, by = colname_subset)
DF_upstream <- cbind(DF_upstream, 
                     data.frame(aa.start = NA, var.position = NA, aa.end = NA))
DF_upstream$VariantHGVSProtein <- DF_upstream$VariantHGVSCoding
DF_upstream$var.position <- gsub("(^c.)([-][[:digit:]]{,3})(.*)","\\2",DF_upstream$VariantHGVSCoding)
DF_upstream$aa.start <- gsub("(^c.[-][[:digit:]]{,3})([[:alpha:]]{1})(.*)","\\2",DF_upstream$VariantHGVSCoding)
DF_upstream$aa.end <- gsub("(^c..*[>}])(.*$)","\\2",DF_upstream$VariantHGVSCoding)

# Intronic Variants
DF_intron <- DF_var[DF_var$intron == TRUE, ]
DF_intron$var.type <- NA
for (row_No in 1:nrow(DF_intron)) {
  if (isTRUE(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+[ATCG]>[ATCG]", DF_intron$coding[row_No]))) {
    DF_intron$var.type[row_No] <- "SNV"
  } else if (isTRUE(DF_intron$fs[row_No])) { DF_intron$var.type[row_No] <- "Frameshift"
  } else if (isTRUE(DF_intron$delins[row_No])) { DF_intron$var.type[row_No] <- "Delins"
  } else if (isTRUE(DF_intron$ins[row_No])) { DF_intron$var.type[row_No] <- "Insertion"
  } else if (isTRUE(DF_intron$del[row_No])) { DF_intron$var.type[row_No] <- "Deletion"
  } else if (isTRUE(DF_intron$dup[row_No])) { DF_intron$var.type[row_No] <- "Duplication"
  }
}
DF_intron <- DF_intron[,c(1:5,17)]
colnames(DF_intron)[1:5] <- colname_subset
DF_intron <- inner_join(DF[DF$VariantLabel %in% DF_intron$VariantLabel, ], 
                        DF_intron, by = colname_subset)
DF_intron <- cbind(DF_intron, 
                   data.frame(aa.start = NA, var.position = NA, aa.end = NA))

# SNV Variants
DF_SNV <- DF_var[(DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE) &
                   DF_var$synonymous == FALSE &
                   DF_var$intron == FALSE, 1:5]
DF_SNV$var.type <- "SNV"

# Frameshift Variants
DF_Frameshift <- DF_var[DF_var$fs == TRUE &
                          DF_var$synonymous == FALSE & DF_var$upstream == FALSE &
                          DF_var$intron == FALSE, ]
DF_Frameshift$var.type <- "Frameshift"
for (row_No in 1:nrow(DF_Frameshift)) {
  if (isTRUE(DF_Frameshift$delins[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Delins"
  } else if (isTRUE(DF_Frameshift$ins[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Insertion"
  } else if (isTRUE(DF_Frameshift$del[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Deletion"
  } else if (isTRUE(DF_Frameshift$dup[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Duplication"
  }
}
DF_Frameshift <- DF_Frameshift[,c(1:5,17)]

# Delins Variants
DF_delins <- DF_var[DF_var$delins == TRUE &
                      DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                      DF_var$fs == FALSE , 1:5]
DF_delins$var.type <- "Delins"

# Insertion Variants
DF_ins <- DF_var[DF_var$ins == TRUE &
                   DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_ins$var.type <- "Insertion"

# Deletion Variants
DF_del <- DF_var[DF_var$del == TRUE &
                   DF_var$ins == FALSE & DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_del$var.type <- "Deletion"

# Duplication Variants
DF_dup <- DF_var[DF_var$dup == TRUE &
                   DF_var$del == FALSE & DF_var$ins == FALSE & DF_var$delins == FALSE &
                   DF_var$intron == FALSE & DF_var$synonymous == FALSE &
                   DF_var$fs == FALSE , 1:5]
DF_dup$var.type <- "Duplication"

# Uncategorized Variants
DF_remain <- DF_var[DF_var$remain == TRUE, 1:5]

cat(paste("All cases have been categorized based on type of nucleotide change, as indicated:",
          nrow(DF) == (nrow(DF_synonymous) + nrow(DF_upstream) + nrow(DF_intron) + nrow(DF_SNV) + 
                         nrow(DF_Frameshift) + nrow(DF_delins) + nrow(DF_ins) + nrow(DF_del) + 
                         nrow(DF_dup)),sep=""),"\n")

# A <- rbind(DF_upstream[,colname_subset],
#            DF_intron[,colname_subset])
# colnames(A) <- c("patientid","variant","gene","coding","protein")
# A$var.type <- NA
# A <- rbind(A,DF_synonymous,DF_SNV,DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup)
# A[duplicated(A[,c("patientid","variant","gene","coding")]),]
# remove(A)

cat(paste("No. entries not categorized based on type of nucleotide change: n=", nrow(DF_remain), sep=""),"\n")

# Merge entries from Synonymous, SNV, Frameshift & In-Frame mutations 
#----------------------------------------------
DF_Map <- rbind(DF_synonymous,DF_SNV, DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup)
colnames(DF_Map)[1:5] <- colname_subset
DF_Map <- inner_join(DF[DF$VariantLabel %in% DF_Map$VariantLabel,], 
                     DF_Map, by = colname_subset)

# Subset entries with missing VariantHGVSProtein field
#----------------------------------------------
DF_NAprotein <- DF_Map[is.na(DF_Map$VariantHGVSProtein),]
DF_NAprotein <- cbind(DF_NAprotein, 
                      data.frame(aa.start = NA, var.position = NA, aa.end = NA))

# Extract STAMP genes with coverage in hotspots 
snv.hotspot.list <- sort(unique(Genes$Name[which(Genes$Coverage == "promoter")]))
assign("snv.hotspot.list", snv.hotspot.list, envir = .GlobalEnv)

# Promoter
row.change = which(DF_NAprotein$VariantGene %in% snv.hotspot.list)
DF_NAprotein$VariantHGVSProtein[row.change] <- "Promoter"

# Subsitute with coding nomenclature
row.change = which(is.na(DF_NAprotein$VariantHGVSProtein))
DF_NAprotein$VariantHGVSProtein[row.change] <- DF_NAprotein$VariantHGVSCoding[row.change]

cat(paste("STAMP entries missing VariantHGVSProtein field: n=",(nrow(DF_NAprotein)), sep=""),"\n")

# FUNCTION: Extract variant positions 
#----------------------------------------------
Extract_VarPosition <- function(DF) {
  DF$aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF$VariantHGVSProtein)
  
  DF$var.position <- as.numeric(gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,4})(.*)", "\\2", 
                                     DF$VariantHGVSProtein))
  
  ## aa.end indicates first amino acid or codon interuption following codon No. 
  ## i.e. frameshift (fs), termination (Ter), deletion (del), duplication (dup)
  for (row_No in 1:nrow(DF)) {
    
    if (isTRUE(DF$var.type[row_No] == "SNV")) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4})([[:alpha:]]{3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
      
    } else if (isTRUE(grepl("Frameshift", DF$var.type[row_No]))) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4})([[:alpha:]]{2,3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
      
    } else if (isTRUE(DF$var.type[row_No] == "Delins")) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4}[_]*[[:digit:]]*[delins]*)([[:alpha:]]{3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
      
    } else if (isTRUE(DF$var.type[row_No] == "Insertion")) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4}[_]*[[:digit:]]*[ins]*)([[:alpha:]]{3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
      
    } else if (isTRUE(DF$var.type[row_No] == "Deletion")) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4}[_]*[[:digit:]]*)([[:alpha:]]{3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
      
    } else if (isTRUE(DF$var.type[row_No] == "Duplication")) {
      DF$aa.end[row_No] <- gsub("(^p.[[:alpha:]]{3}[[:digit:]]{,4}[_]*)([[:alpha:]]{3})(.*)", "\\2", 
                                DF$VariantHGVSProtein[row_No])
    }
  }
  
  assign("DF", DF, envir = .GlobalEnv)
} 

# Extract variant positions for entries with VariantHGVSProtein field (aa.start, var.position, aa.end)
#----------------------------------------------
DF_Map <- DF_Map[!is.na(DF_Map$VariantHGVSProtein),]
Extract_VarPosition(DF = DF_Map)
DF_Map <- DF
cat(paste("STAMP entries mapped to lollipop plots: n=",(nrow(DF_Map[DF_Map$var.type != "Synonymous",])), sep=""),"\n")

row.change = which(is.na(DF_intron$VariantHGVSProtein))
DF_intron$VariantHGVSProtein[row.change] <- DF_intron$VariantHGVSCoding[row.change]

# Merge all entries 
#----------------------------------------------
DF_Full <- rbind(DF_Map, DF_NAprotein, DF_intron, DF_upstream)

# Classification of variants for clinical trial matching (var.anno)
#----------------------------------------------
# Frameshift are a special type of amino acid deletion/insertion
# Source: http://varnomen.hgvs.org/recommendations/protein/variant/frameshift/
DF_Full$var.anno <- NA
for (row_No in 1:nrow(DF_Full)) {
  if (DF_Full$var.type[row_No] == "Synonymous") {
    DF_Full$var.anno[row_No] <- "OTHER"
  } else {
    DF_Full$var.anno[row_No] <- "MUTATION"
  }}

# PrimaryTumorSite grouped according to medical specializations (PrimaryTumorSite.Category)
#----------------------------------------------
# Confirm PrimaryTumorSite is classified 
primaryTumorSite.STAMP <- sort(unique(DF_Full$PrimaryTumorSite))
primaryTumorSite.key <- sort(unique(DiseaseGroupCategory_LongFormat$primaryTumorSite))

# for (elem_No in 1:length(primaryTumorSite.STAMP)) {
#   if (isTRUE(is.element(primaryTumorSite.STAMP[elem_No], primaryTumorSite.key) == FALSE)) {
#     cat("\n")
#     cat(paste(primaryTumorSite.STAMP[elem_No], 
#               ": primaryTumorSite has not been classified into DiseaseGroupCategory", sep=""),"\n","\n")
#   }
# }

DF_Full$PrimaryTumorSite.Category <- NA
for (row_No in 1:nrow(DF_Full)) {
  DiseaseGroupCategory.name <- 
    DiseaseGroupCategory_LongFormat$Disease.Group.category[which(DiseaseGroupCategory_LongFormat$primaryTumorSite == 
                                                                   DF_Full$PrimaryTumorSite[row_No])] 
  
  DF_Full$PrimaryTumorSite.Category[row_No] <- paste(as.character(DiseaseGroupCategory.name), collapse=", ")
}

# If missing smpl.primaryTumorSite, primaryTumorSite.category = "unknown"
DF_Full$PrimaryTumorSite.Category[which(DF_Full$PrimaryTumorSite.Category == "")] <- "unknown"

# Histological Dx grouped according to CTEP categories
#----------------------------------------------
# Confirm HistologicalDx is classified 
HistologicalDx.STAMP <- sort(unique(DF_Full$HistologicalDx))
HistologicalDx.key <- sort(unique(HistologicalDxCategory$histologicalDiagnosis))

# for (elem_No in 1:length(HistologicalDx.STAMP)) {
#   if (isTRUE(is.element(HistologicalDx.STAMP[elem_No], HistologicalDx.key) == FALSE)) {
#     cat("\n")
#     cat(paste(HistologicalDx.STAMP[elem_No], 
#               ": Histological Dx has not been classified into HistologicalDxCategory", sep=""),"\n","\n")
#   }
# }
# 
# # Identify entries with incorrect HGVS genomic nomenclature
# #----------------------------------------------
# for (row_No in 1:nrow(DF_Full)) {
#   if (isTRUE(grepl("^chr[[:alnum:]]{,2}:g.",DF_Full$VariantHGVSGenomic[row_No]) == FALSE |
#              grepl("[[:alpha:]]+$",DF_Full$VariantHGVSGenomic[row_No]) == FALSE)) {
#     cat("\n")
#     cat(paste(DF_Full$PatientID[row_No], " has incorrect HGVS genomic nomenclature: ", sep=""),"\n")
#     print(DF_Full[row_No, 
#                   c("VariantNMAccession","VariantCHR","VariantHGVSGenomic","VariantLabel",
#                     "VariantGene","VariantHGVSCoding","VariantHGVSProtein","VariantPathogenicityStatus")],
#           row.names = FALSE)
#   }
# }

# Input Exon Number
DF_Full$Exon_Number <- NA
for (row_No in 1:nrow(DF_Full)) {
  gene_id <- DF_Full$VariantGene[row_No]
  genomic_pos <- gsub("(^chr[[:digit:]]{,2}:g.)([[:digit:]]+)([_]*[[:digit:]]*[[:alpha:]]+.*)", "\\2", DF_Full$VariantHGVSGenomic[row_No])
  
  Gene.ExonTable.list <- sort(unique(stamp_reference_full$Gene))
  # Remove genes with promoter coverage 
  Gene.ExonTable.list <- Gene.ExonTable.list[!(Gene.ExonTable.list %in% snv.hotspot.list)]
  
  if (isTRUE(gene_id %in% Gene.ExonTable.list &
             # HGVS genomic region CANNOT be NA
             !is.na(genomic_pos))) {
    Gene.ExonTable <- stamp_reference_full[stamp_reference_full$Gene == gene_id,]
    
    DF_Full$Exon_Number[row_No] <- min(Gene.ExonTable$exon_number[which(Gene.ExonTable$start <= genomic_pos)])
  }
  
  remove(Gene.ExonTable.list)
}

cat(paste("No. of genes without annotated exon information: ",
    length(sort(unique(DF_Full$VariantGene[is.na(DF_Full$Exon_Number)]))), sep=""),"\n")
print(sort(unique(DF_Full$VariantGene[is.na(DF_Full$Exon_Number)])))

## Remove entries with missing information
#----------------------------------------------
## Remove entries with Synonymous mutations = 26 entries
DF_Full <- DF_Full[which(DF_Full$var.anno != "OTHER"),]

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_DF", DF_Full, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(DF_Full, file = paste(tempdir, Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Export entries per patient into "_SNVIndel.tsv" file
#----------------------------------------------
patient.list <- sort(unique(STAMP_DF$PatientID))

for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Extract STAMP entries for individual patient
  DF_patient <- STAMP_DF[which(STAMP_DF$PatientID == patient_id),]
  
  # Order alphabetically by gene name
  DF_patient <- DF_patient[order(DF_patient$VariantGene, decreasing = FALSE),]
  
  ## Write to local computer
  write.table(DF_patient, file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

## Output potential matching errors due to structurization
#----------------------------------------------
# Missing information
out.DF <- STAMP_DF[is.na(STAMP_DF$VariantGene),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing gene information:", "\n")
  print(unique(out.DF$PatientID))
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSProtein)),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing HGVS protein information:", "\n")
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantGene, sep=" : ")
  print(sort(unique(out.DF$print)))
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSCoding)),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing HGVS coding information:", "\n")
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantGene, sep=" : ")
  print(sort(unique(out.DF$print)))
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSGenomic)),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing HGVS genomic information:", "\n")
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantGene, sep=" : ")
  print(sort(unique(out.DF$print)))
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$var.anno)),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing variant annotation:", "\n")
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantLabel, sep=" : ")
  print(sort(unique(out.DF$print)))
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$PrimaryTumorSite.Category)),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with missing primary tumor site category:", "\n")
  out.DF$print <- paste(out.DF$PatientID, " - primary tumor site: ", out.DF$PrimaryTumorSite, sep="")
  print(sort(unique(out.DF$print)))
  cat("\n")
}

# Incorrect structure
stamp_reference_gene.list <- sort(unique(stamp_reference_full$Gene))
out.DF <- STAMP_DF[!(STAMP_DF$VariantGene %in% stamp_reference_gene.list),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with gene not listed in stamp reference table:", "\n")
  out.DF <- out.DF[, c("PatientID","VariantGene")]
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantGene, sep=" : ")
  print(out.DF$print , quote = TRUE)
  cat("\n")
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSProtein) &
                           grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*", STAMP_DF$VariantHGVSProtein) == FALSE),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with HGVS protein nomenclature formatted incorrectly:", "\n")
  out.DF <- out.DF[, c("PatientID","VariantHGVSProtein")]
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantHGVSProtein, sep=" : ")
  print(out.DF$print , quote = TRUE)
  cat("\n")
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSCoding) &
                           grepl("^c.[-]*[[:digit:]]+.*", STAMP_DF$VariantHGVSCoding) == FALSE),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with HGVS coding nomenclature formatted incorrectly:", "\n")
  out.DF <- out.DF[, c("PatientID","VariantHGVSCoding")]
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantHGVSCoding, sep=" : ")
  print(out.DF$print , quote = TRUE)
  cat("\n")
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSGenomic) &
                           grepl("^chr[[:alnum:]]{,2}:g.[[:digit:]]+.*", STAMP_DF$VariantHGVSGenomic) == FALSE),]
if (nrow(out.DF) > 0) {
  cat("Patient IDs with HGVS genomic nomenclature formatted incorrectly:", "\n")
  out.DF <- out.DF[, c("PatientID","VariantHGVSGenomic")]
  out.DF$print <- paste(out.DF$PatientID, out.DF$VariantHGVSGenomic, sep=" : ")
  print(out.DF$print , quote = TRUE)
  cat("\n")
}

remove(DF,DF_del,DF_delins,DF_dup,DF_Frameshift,DF_Full,DF_ins,DF_intron,DF_Map,DF_NAprotein,
       DF_remain,DF_SNV,DF_synonymous,DF_upstream,DF_var,colname_subset,DiseaseGroupCategory.name,
       i,primaryTumorSite.key,primaryTumorSite.STAMP,row_No,Extract_VarPosition,
       patient_id,patient_num,patient.list,DF_patient,HistologicalDx.key,HistologicalDx.STAMP,
       gene_id,genomic_pos,Gene.ExonTable,out.DF)

cat("\n")
