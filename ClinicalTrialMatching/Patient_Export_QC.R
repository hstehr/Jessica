## Classification #1 (var.type): SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Classification #2 (var.anno): "MUTATION"
## Classification #3 (Exon_Number)
## Output: "patient_id_QC.tsv"

# Load relevant file from global environment 
STAMP_DF_structured <- 
  cbind(data.frame(sys.uniqueId = patient.id, stringsAsFactors = FALSE), STAMP_DF)

# Subset for variants with status = "ACCEPT"
STAMP_DF_structured <- STAMP_DF_structured[which(STAMP_DF_structured$Status == "ACCEPT"),]

# Restructure STAMP dataframe
#----------------------------------------------
if (nrow(STAMP_DF_structured) > 0) {
  # Modify nomenclature
  STAMP_DF_structured$CDS.Change <- paste("c.",STAMP_DF_structured$CDS.Change, sep="")
  STAMP_DF_structured$AA.Change <- paste("p.",STAMP_DF_structured$AA.Change, sep="")
  STAMP_DF_structured$hgvsGenomic <- 
    paste("g.",STAMP_DF_structured$Position,STAMP_DF_structured$Ref,">",STAMP_DF_structured$Var,sep="")
  
  # Convert codon nomenclature from 1-letter to 3-letters 
  for (row_No in 1:nrow(STAMP_DF_structured)) {
    chars <- unlist(strsplit(STAMP_DF_structured$AA.Change[row_No], ""))
    for (elem_No in 1:length(chars)) {
      if (isTRUE(grepl("[[:upper:]]", chars[elem_No]))) {
        chars[elem_No] <- AminoAcid_Conversion$Code3[which(AminoAcid_Conversion$Code1 == chars[elem_No])]
      }
    }
    STAMP_DF_structured$AA.Change[row_No] <-  paste0(chars, collapse = "")
  }
  
  # Convert asterisk to "Ter"
  STAMP_DF_structured$AA.Change <- gsub("\\*", "Ter", STAMP_DF_structured$AA.Change)
  
  STAMP_DF_structured$sys.label <- 
    paste(STAMP_DF_structured$Gene, " ", STAMP_DF_structured$CDS.Change, " (", STAMP_DF_structured$AA.Change, ")", sep="")
  
  # Subset columns of interest
  #----------------------------------------------
  colnames_keep <- c("sys.uniqueId","Chr","hgvsGenomic","sys.label","Gene","CDS.Change","AA.Change")
  STAMP_DF_structured <- STAMP_DF_structured[,colnames_keep]
  
  colnames_generic <- c("PatientID","VariantCHR","VariantHGVSGenomic","VariantLabel","VariantGene",
                        "VariantHGVSCoding","VariantHGVSProtein")
  colnames(STAMP_DF_structured) <- colnames_generic
  
  # Categorize protein sequence variants (var.type)
  #----------------------------------------------
  DF_var <- data.frame(PatientID=STAMP_DF_structured$PatientID,
                       VariantCHR=STAMP_DF_structured$VariantCHR,
                       VariantHGVSGenomic=STAMP_DF_structured$VariantHGVSGenomic,
                       VariantLabel=STAMP_DF_structured$VariantLabel,
                       VariantGene=STAMP_DF_structured$VariantGene,
                       VariantHGVSCoding=STAMP_DF_structured$VariantHGVSCoding, 
                       VariantHGVSProtein=STAMP_DF_structured$VariantHGVSProtein,
                       
                       # Categorize based on HGVS.Protein
                       SNVprotein=grepl("^p.[[:upper:]]{1}[[:lower:]]{2}[[:digit:]]+[[:upper:]]{1}[[:lower:]]{2}",
                                        STAMP_DF_structured$VariantHGVSProtein),
                       fs=grepl("fs", STAMP_DF_structured$VariantHGVSProtein),
                       
                       # Categorize based on HGVS.Coding
                       intron=grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+.*", STAMP_DF_structured$VariantHGVSCoding),
                       SNVcoding=grepl("^c.[[:digit:]]+[ATCG]>[ATCG]", STAMP_DF_structured$VariantHGVSCoding),
                       ins=grepl("ins", STAMP_DF_structured$VariantHGVSCoding),
                       del=grepl("del", STAMP_DF_structured$VariantHGVSCoding),
                       dup=grepl("dup", STAMP_DF_structured$VariantHGVSCoding),
                       
                       # Categorize based on Variant.Label
                       delins=grepl("del.*ins", STAMP_DF_structured$VariantLabel),
                       
                       # Empty variable for unclassified variant
                       remain=NA,
                       stringsAsFactors = FALSE)
  
  for (i in 1:nrow(DF_var)) {
    if (length(grep(FALSE, DF_var[i,8:15])) == 8) { DF_var$remain[i] <- TRUE
    } else { DF_var$remain[i] <- FALSE
    }}
  
  # Classification correction
  DF_var$SNVprotein[which(DF_var$SNVprotein == TRUE & 
                            (DF_var$delins == TRUE | DF_var$ins == TRUE | 
                               DF_var$del == TRUE | DF_var$dup == TRUE))] <- FALSE
  
  # Subset based on matching of strings in smpl.hgvsCoding or smpl.hgvsProtein
  # Classification order: SNV > Frameshift > In-Frame (delins > insertion > deletion > duplication)
  #----------------------------------------------
  # Intronic Variants
  DF_intron <- DF_var[DF_var$intron == TRUE, ]
  if (isTRUE(exists("DF_intron") & nrow(DF_intron) > 0)) {
    DF_intron$var.type <- NA
    for (row_No in 1:nrow(DF_intron)) {
      if (isTRUE(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+[ATCG]>[ATCG]", DF_intron$VariantHGVSCoding[row_No]))) {
        DF_intron$var.type[row_No] <- "SNV"
      } else if (isTRUE(DF_intron$fs[row_No])) { DF_intron$var.type[row_No] <- "Frameshift"
      } else if (isTRUE(DF_intron$delins[row_No])) { DF_intron$var.type[row_No] <- "Delins"
      } else if (isTRUE(DF_intron$ins[row_No])) { DF_intron$var.type[row_No] <- "Insertion"
      } else if (isTRUE(DF_intron$del[row_No])) { DF_intron$var.type[row_No] <- "Deletion"
      } else if (isTRUE(DF_intron$dup[row_No])) { DF_intron$var.type[row_No] <- "Duplication"
      }
      
      if (isTRUE(grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*$", DF_intron$VariantHGVSProtein[row_No]) == FALSE)) {
        DF_intron$VariantHGVSProtein[row_No] <- NA
      }
    }
    
    DF_intron <- DF_intron[,c(1:7,17)]
    DF_intron <- cbind(DF_intron, 
                       data.frame(aa.start = NA, var.position = NA, aa.end = NA))
  }
  
  # SNV Variants
  DF_SNV <- DF_var[DF_var$SNVprotein == TRUE | DF_var$SNVcoding == TRUE, 1:7]
  if (isTRUE(exists("DF_SNV") & nrow(DF_SNV) > 0)) {DF_SNV$var.type <- "SNV"}
  
  # Frameshift Variants
  DF_Frameshift <- DF_var[DF_var$fs == TRUE & DF_var$intron == FALSE, ]
  if (isTRUE(exists("DF_Frameshift") & nrow(DF_Frameshift) > 0)) { 
    DF_Frameshift$var.type <- "Frameshift"
    for (row_No in 1:nrow(DF_Frameshift)) {
      if (isTRUE(DF_Frameshift$delins[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Delins"
      } else if (isTRUE(DF_Frameshift$ins[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Insertion"
      } else if (isTRUE(DF_Frameshift$del[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Deletion"
      } else if (isTRUE(DF_Frameshift$dup[row_No])) { DF_Frameshift$var.type[row_No] <- "Frameshift_Duplication"
      }}
    DF_Frameshift <- DF_Frameshift[,c(1:7,17)]
  }
  
  # Delins Variants
  DF_delins <- DF_var[DF_var$delins == TRUE & DF_var$intron == FALSE & 
                        DF_var$fs == FALSE , 1:7]
  if (isTRUE(exists("DF_delins") & nrow(DF_delins) > 0)) { 
    DF_delins$var.type <- "Delins"
  }
  
  # Insertion Variants
  DF_ins <- DF_var[DF_var$ins == TRUE & DF_var$delins == FALSE & 
                     DF_var$intron == FALSE & DF_var$fs == FALSE , 1:7]
  if (isTRUE(exists("DF_ins") & nrow(DF_ins) > 0)) { 
    DF_ins$var.type <- "Insertion"
  }
  
  # Deletion Variants
  DF_del <- DF_var[DF_var$del == TRUE & DF_var$ins == FALSE & 
                     DF_var$intron == FALSE & DF_var$delins == FALSE & 
                     DF_var$fs == FALSE, 1:7]
  if (isTRUE(exists("DF_del") & nrow(DF_del))) { 
    DF_del$var.type <- "Deletion"
  }
  
  # Duplication Variants
  DF_dup <- DF_var[DF_var$dup == TRUE & DF_var$del == FALSE & 
                     DF_var$ins == FALSE & DF_var$delins == FALSE & 
                     DF_var$intron == FALSE & DF_var$fs == FALSE , 1:7]
  if (isTRUE(exists("DF_dup") & nrow(DF_dup) > 0)) { 
    DF_dup$var.type <- "Duplication"
  }
  
  # Uncategorized Variants
  DF_remain <- DF_var[DF_var$remain == TRUE, 1:7]
  cat(paste("No. entries not categorized based on type of nucleotide change: n=", nrow(DF_remain), sep=""),"\n")
  
  # Merge entries from SNV, Frameshift & In-Frame mutations 
  #----------------------------------------------
  STAMP_DF_structured_anno <- rbind(DF_SNV, DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup)
  
  # Subset entries with missing VariantHGVSProtein field
  #----------------------------------------------
  DF_NAprotein <- STAMP_DF_structured_anno[is.na(STAMP_DF_structured_anno$VariantHGVSProtein),]
  cat(paste("STAMP entries missing VariantHGVSProtein field: n=",(nrow(DF_NAprotein)), sep=""),"\n")
  
  # FUNCTION: Extract variant positions 
  #----------------------------------------------
  Extract_VarPosition <- function(DF) {
    DF$aa.start <- gsub("(^p.)([[:alpha:]]{3})(.*)", "\\2", DF$VariantHGVSProtein)
    
    DF$var.position <- gsub("(^p.[[:alpha:]]{3})([[:digit:]]{,4})(.*)", "\\2", 
                            DF$VariantHGVSProtein)
    
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
    
    for (row_No in 1:nrow(DF)) {
      if (isTRUE(grepl("^[[:alpha:]]{,3}$", DF$aa.start[row_No]) == FALSE)) {
        DF$aa.start[row_No] <- NA
      }
      
      if (isTRUE(grepl("^[[:digit:]]+$", DF$var.position[row_No]) == FALSE)) {
        DF$var.position[row_No] <- NA
      } else {
        DF$var.position[row_No] <- as.numeric(DF$var.position[row_No] )
      }
      
      if (isTRUE(grepl("^[[:alpha:]]{,3}$", DF$aa.end[row_No]) == FALSE)) {
        DF$aa.end[row_No] <- NA
      }
      
      if (isTRUE(grepl("^[[:alpha:]]{,3}$", DF$aa.end[row_No]) == FALSE)) {
        DF$aa.end[row_No] <- NA
      }
      
      if (isTRUE(grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*$", DF$VariantHGVSProtein[row_No]) == FALSE)) {
        DF$VariantHGVSProtein[row_No] <- NA
      }
    }
    
    assign("DF", DF, envir = .GlobalEnv)
  } 
  
  # Extract variant positions for entries with VariantHGVSProtein field (aa.start, var.position, aa.end)
  #----------------------------------------------
  STAMP_DF_structured_anno <- STAMP_DF_structured_anno[!is.na(STAMP_DF_structured_anno$VariantHGVSProtein),]
  Extract_VarPosition(DF = STAMP_DF_structured_anno)
  STAMP_DF_structured_anno <- DF
  STAMP_DF_structured_anno <- rbind(STAMP_DF_structured_anno, DF_intron)
  
  # Modify if changes made to nomenclature
  STAMP_DF_structured_anno$VariantLabel <- paste(STAMP_DF_structured_anno$VariantGene, " ",
                                                 STAMP_DF_structured_anno$VariantHGVSCoding, " (",
                                                 STAMP_DF_structured_anno$VariantHGVSProtein, ")",sep="")
  
  # Classification of variants for clinical trial matching (var.anno)
  #----------------------------------------------
  STAMP_DF_structured_anno$var.anno <- "MUTATION"
  
  # Input Exon Number
  #----------------------------------------------
  for (row_No in 1:nrow(STAMP_DF_structured_anno)) {
    gene_id <- STAMP_DF_structured_anno$VariantGene[row_No]
    genomic_pos <- gsub("(^g.)([[:digit:]]+)([_]*[[:digit:]]*[[:alpha:]]+.*)", "\\2", 
                        STAMP_DF_structured_anno$VariantHGVSGenomic[row_No])
    
    Gene.ExonTable <- stamp_reference_full[stamp_reference_full$Gene == gene_id,]
    
    for (exon_row_No in 1:nrow(Gene.ExonTable)) {
      exon_start <- Gene.ExonTable$start[exon_row_No]
      exon_end <- Gene.ExonTable$end[exon_row_No]
      
      if (isTRUE(genomic_pos >= exon_start & genomic_pos <= exon_end)) {
        STAMP_DF_structured_anno$Exon_Number[row_No] <- Gene.ExonTable$exon_number[exon_row_No]
      }
    }
  }
  
  # Missing information
  #----------------------------------------------
  STAMP_DF_structured_anno$PrimaryTumorSite.Category <- "unknown"
  STAMP_DF_structured_anno$PrimaryTumorSite <- "unknown"
  STAMP_DF_structured_anno$VariantPathogenicityStatus <- "NULL"
  
  ## Write to local computer
  #----------------------------------------------
  write.table(STAMP_DF_structured_anno, file = paste(tempdir, patient.id, "_QC.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  remove(DF,DF_del,DF_delins,DF_dup,DF_Frameshift,DF_ins,DF_intron,DF_NAprotein,
         DF_remain,DF_SNV,DF_var,STAMP_DF_structured,chars,colnames_generic,colnames_keep,
         elem_No,exon_end,exon_row_No,exon_start,gene_id,genomic_pos,i,row_No,Extract_VarPosition)
  
} else {
  STAMP_DF_structured_anno <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(STAMP_DF_structured_anno) <- c("PatientID","VariantCHR","VariantHGVSGenomic",
                                          "VariantLabel","VariantGene","VariantHGVSCoding",
                                          "VariantHGVSProtein","var.type","aa.start","var.position",
                                          "aa.end","var.anno","Exon_Number","PrimaryTumorSite.Category","PrimaryTumorSite",
                                          "VariantPathogenicityStatus")
}

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_DF", STAMP_DF_structured_anno, envir = .GlobalEnv)
remove(STAMP_DF_structured_anno)

## Export entries per patient into .tsv file
#----------------------------------------------
# Extract STAMP entries for individual patient
DF_patient <- STAMP_DF[which(STAMP_DF$PatientID == patient.id),]

# Order alphabetically by gene name
DF_patient <- DF_patient[order(DF_patient$VariantGene, decreasing = FALSE),]

## Write to local computer
write.table(DF_patient, file = paste(tempdir, patient.id, ".tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

if (nrow(DF_patient) == 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat(paste(patient.id, " has no entries with variants with ACCEPTED status.", sep=""),"\n","\n")
  sink()
} 
remove(DF_patient)

## Output potential matching errors due to structurization
#----------------------------------------------
# Missing information
out.DF <- STAMP_DF[is.na(STAMP_DF$VariantGene),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Patient has entry with missing gene information", "\n","\n")
  sink()
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSProtein)),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with missing HGVS protein information:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSCoding)),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with missing HGVS coding information:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSGenomic)),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with missing HGVS genomic information:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}

# Incorrect structure
stamp_reference_gene.list <- sort(unique(stamp_reference_full$Gene))
out.DF <- STAMP_DF[!(STAMP_DF$VariantGene %in% stamp_reference_gene.list),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with gene not listed in stamp reference table:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}
remove(stamp_reference_gene.list)

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSProtein) &
                           grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*", STAMP_DF$VariantHGVSProtein) == FALSE),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with HGVS protein nomenclature formatted incorrectly:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSCoding) &
                           grepl("^c.[-]*[[:digit:]]+.*", STAMP_DF$VariantHGVSCoding) == FALSE),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with HGVS coding nomenclature formatted incorrectly:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSGenomic) &
                           grepl("^g.[[:digit:]]+.*", STAMP_DF$VariantHGVSGenomic) == FALSE),]
if (nrow(out.DF) > 0) {
  sink(file = err.output, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Entry with HGVS genomic nomenclature formatted incorrectly:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
  
  sink()
}
remove(out.DF)

cat("\n")
