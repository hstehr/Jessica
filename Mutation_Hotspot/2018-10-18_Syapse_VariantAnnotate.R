## Classification (var.type): Synonymous, Upstream, Intronic, SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Output: "Syapse_Export_DF_STAMP_4Map.tsv"

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
  }}

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
  }}
DF_upstream <- DF_upstream[,c(1:5,17)]
colnames(DF_upstream)[1:5] <- colname_subset
DF_upstream <- inner_join(DF[DF$VariantLabel %in% DF_upstream$VariantLabel, ], 
                          DF_upstream, by = colname_subset)
DF_upstream <- cbind(DF_upstream, 
                     data.frame(aa.start = NA, var.position = NA, aa.end = NA))

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
  }}
DF_intron <- DF_intron[,c(1:5,17)]
colnames(DF_intron)[1:5] <- colname_subset
DF_intron <- inner_join(DF[DF$VariantLabel %in% DF_intron$VariantLabel, ], 
                        DF_intron, by = colname_subset)
DF_intron <- cbind(DF_intron, 
                   data.frame(aa.start = NA, var.position = NA, aa.end = NA))

# SNV Variants
DF_SNV <- DF_var[DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE &
                   DF_var$synonymous == FALSE, 1:5]
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
  }}
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
DF_Map <- DF_Map[DF_Map$var.type != "Synonymous",]

## Write to local computer
#----------------------------------------------
write.table(DF_Map, file = paste(Syapse_Export_timestamp, "_Syapse_Export_DF_STAMP_4Map.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(DF,DF_del,DF_delins,DF_dup,DF_Frameshift,DF_ins,DF_intron,DF_remain,DF_SNV,
       DF_synonymous,DF_upstream,DF_var,colname_subset,i,row_No,Extract_VarPosition)
