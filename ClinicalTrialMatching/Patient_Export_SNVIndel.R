## Classification #1 (var.type): SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Classification #2 (var.anno): "MUTATION"
## Classification #3 (Exon_Number)
## Output: "patient_id_SNVIndel.tsv"

cat(paste("Timestamp of patient STAMP processing START: ", Sys.time(), sep=""),"\n","\n")

# Load relevant file from global environment 
STAMP_DF_structured <- 
  cbind(data.frame(sys.uniqueId = patient.id, stringsAsFactors = FALSE), STAMP_DF)

# Subset for variants with status != "NOT REPORTED"
STAMP_DF_structured <- STAMP_DF_structured[which(STAMP_DF_structured$Status != "NOT_REPORTED"),]

# Restructure STAMP dataframe
#----------------------------------------------
if (nrow(STAMP_DF_structured) == 0) {
  
  # Generate empty dataframe
  STAMP_DF_structured_anno <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(STAMP_DF_structured_anno) <- c("PatientID","VariantHGVSGenomic","VariantLabel","VariantGene","VariantHGVSCoding",
                                          "VariantHGVSProtein","var.type","var.anno","Exon_Number",
                                          "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus")
  
  # Output indication to file
  cat(paste(patient.id, " does not have any SNV/Indel entries with statuses that are not \"NOT REPORTED \".", sep=""),"\n","\n")
  
} else if (nrow(STAMP_DF_structured) > 0) {
  
  # Modify nomenclature
  STAMP_DF_structured$CDS.Change <- paste("c.",STAMP_DF_structured$CDS.Change, sep="")
  STAMP_DF_structured$CDS.Change[which(grepl("^c..$",STAMP_DF_structured$CDS.Change))] <- NA
  
  STAMP_DF_structured$AA.Change <- paste("p.",STAMP_DF_structured$AA.Change, sep="")
  STAMP_DF_structured$AA.Change[which(grepl("^p..$",STAMP_DF_structured$AA.Change))] <- NA
  
  STAMP_DF_structured$hgvsGenomic <- 
    paste(STAMP_DF_structured$Chr,":g.",STAMP_DF_structured$Position,STAMP_DF_structured$Ref,">",STAMP_DF_structured$Var,sep="")
  
  STAMP_DF_structured$sys.label <- 
    paste(STAMP_DF_structured$Gene, " ", STAMP_DF_structured$CDS.Change, " (", STAMP_DF_structured$AA.Change, ")", sep="")
  
  # Subset columns of interest
  #----------------------------------------------
  colnames_keep <- c("sys.uniqueId","hgvsGenomic","sys.label","Gene","CDS.Change","AA.Change")
  STAMP_DF_structured <- STAMP_DF_structured[,colnames_keep]
  
  colnames_generic <- c("PatientID","VariantHGVSGenomic","VariantLabel","VariantGene",
                        "VariantHGVSCoding","VariantHGVSProtein")
  colnames(STAMP_DF_structured) <- colnames_generic
  
  # Categorize protein sequence variants (var.type)
  #----------------------------------------------
  DF_var <- data.frame(PatientID=STAMP_DF_structured$PatientID,
                       VariantHGVSGenomic=STAMP_DF_structured$VariantHGVSGenomic,
                       VariantLabel=STAMP_DF_structured$VariantLabel,
                       VariantGene=STAMP_DF_structured$VariantGene,
                       VariantHGVSCoding=STAMP_DF_structured$VariantHGVSCoding, 
                       VariantHGVSProtein=STAMP_DF_structured$VariantHGVSProtein,
                       
                       # Categorize based on HGVS.Protein
                       SNV_Protein_Genomic=(grepl("^p.[[:upper:]]{1}[[:digit:]]+[[:upper:]]{1}",
                                                  STAMP_DF_structured$VariantHGVSProtein) | 
                                              grepl("^chr[[:alnum:]]+:g.[[:digit:]]+[[:upper:]]{1}>[[:upper:]]{1}",
                                                    STAMP_DF_structured$VariantHGVSGenomic)),
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
    if (length(grep(FALSE, DF_var[i,7:14])) == 8) { DF_var$remain[i] <- TRUE
    } else { DF_var$remain[i] <- FALSE
    }
  }
  
  # Classification correction
  DF_var$SNV_Protein_Genomic[which(DF_var$SNV_Protein_Genomic == TRUE & 
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
    
    DF_intron <- DF_intron[,c(1:6,16)]
    DF_intron <- cbind(DF_intron, 
                       data.frame(aa.start = NA, var.position = NA, aa.end = NA))
  }
  
  # SNV Variants
  DF_SNV <- DF_var[DF_var$SNV_Protein_Genomic == TRUE | DF_var$SNVcoding == TRUE, 1:6]
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
      }
    }
    DF_Frameshift <- DF_Frameshift[,c(1:6,16)]
  }
  
  # Delins Variants
  DF_delins <- DF_var[DF_var$delins == TRUE & DF_var$intron == FALSE & 
                        DF_var$fs == FALSE , 1:6]
  if (isTRUE(exists("DF_delins") & nrow(DF_delins) > 0)) { 
    DF_delins$var.type <- "Delins"
  }
  
  # Insertion Variants
  DF_ins <- DF_var[DF_var$ins == TRUE & DF_var$delins == FALSE & 
                     DF_var$intron == FALSE & DF_var$fs == FALSE , 1:6]
  if (isTRUE(exists("DF_ins") & nrow(DF_ins) > 0)) { 
    DF_ins$var.type <- "Insertion"
  }
  
  # Deletion Variants
  DF_del <- DF_var[DF_var$del == TRUE & DF_var$ins == FALSE & 
                     DF_var$intron == FALSE & DF_var$delins == FALSE & 
                     DF_var$fs == FALSE, 1:6]
  if (isTRUE(exists("DF_del") & nrow(DF_del))) { 
    DF_del$var.type <- "Deletion"
  }
  
  # Duplication Variants
  DF_dup <- DF_var[DF_var$dup == TRUE & DF_var$del == FALSE & 
                     DF_var$ins == FALSE & DF_var$delins == FALSE & 
                     DF_var$intron == FALSE & DF_var$fs == FALSE , 1:6]
  if (isTRUE(exists("DF_dup") & nrow(DF_dup) > 0)) { 
    DF_dup$var.type <- "Duplication"
  }
  
  # Uncategorized Variants
  DF_remain <- DF_var[DF_var$remain == TRUE, 1:6]
  cat(paste("No. entries not categorized based on type of nucleotide change: n=", nrow(DF_remain), sep=""),"\n","\n")
  
  # Merge entries from SNV, Frameshift & In-Frame mutations 
  #----------------------------------------------
  STAMP_DF_structured_anno <- rbind(DF_SNV, DF_Frameshift,DF_delins,DF_ins,DF_del,DF_dup)
  
  # Classification of variants for clinical trial matching (var.anno)
  #----------------------------------------------
  STAMP_DF_structured_anno$var.anno <- "MUTATION"
  
  # Input Exon Number
  #----------------------------------------------
  for (row_No in 1:nrow(STAMP_DF_structured_anno)) {
    gene_id <- STAMP_DF_structured_anno$VariantGene[row_No]
    genomic_pos <- gsub("(^chr[[:digit:]]+[:]g.)([[:digit:]]+)([_]*[[:digit:]]*[[:alpha:]]+.*)", "\\2", 
                        STAMP_DF_structured_anno$VariantHGVSGenomic[row_No])
    
    Gene.ExonTable <- stamp_reference_full[stamp_reference_full$Gene == gene_id,]
    
    for (exon_row_No in 1:nrow(Gene.ExonTable)) {
      exon_start <- Gene.ExonTable$start[exon_row_No]
      exon_end <- Gene.ExonTable$end[exon_row_No]
      
      if (isTRUE(genomic_pos >= exon_start & genomic_pos <= exon_end)) {
        STAMP_DF_structured_anno$Exon_Number[row_No] <- Gene.ExonTable$exon_number[exon_row_No]
      }
    }
    remove(Gene.ExonTable)
  }
  
  # Missing information
  #----------------------------------------------
  STAMP_DF_structured_anno$PrimaryTumorSite.Category <- "unknown"
  STAMP_DF_structured_anno$PrimaryTumorSite <- "unknown"
  STAMP_DF_structured_anno$VariantPathogenicityStatus <- "NULL"
  
  remove(DF_del,DF_delins,DF_dup,DF_Frameshift,DF_ins,DF_intron,
         DF_remain,DF_SNV,DF_var,STAMP_DF_structured,colnames_generic,colnames_keep,
         exon_end,exon_row_No,exon_start,gene_id,genomic_pos,i,row_No)
} 

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_DF", STAMP_DF_structured_anno, envir = .GlobalEnv)

## Export entries per patient into .tsv file
#----------------------------------------------
# Extract STAMP entries for individual patient
DF_patient <- STAMP_DF[which(STAMP_DF$PatientID == patient.id),]

# Order alphabetically by gene name
DF_patient <- DF_patient[order(DF_patient$VariantGene, decreasing = FALSE),]

## Write to local computer
write.table(DF_patient, file = paste(tempdir, patient.id, "_SNVIndel.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

## Output potential matching errors due to structurization
#----------------------------------------------
# Missing information
out.DF <- STAMP_DF[is.na(STAMP_DF$VariantGene),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries missing VariantGene field:", "\n","\n")
  print(out.DF, row.names = FALSE, quote = FALSE)
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSProtein)),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries missing VariantHGVSProtein field:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSCoding)),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries missing VariantHGVSCoding field:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

out.DF <- STAMP_DF[which(is.na(STAMP_DF$VariantHGVSGenomic)),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries missing VariantHGVSGenomic field:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

# Incorrect structure
stamp_reference_gene.list <- sort(unique(stamp_reference_full$Gene))
out.DF <- STAMP_DF[!(STAMP_DF$VariantGene %in% stamp_reference_gene.list),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries with Gene not listed in stamp reference table:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSCoding) &
                           grepl("^c.[-]*[[:digit:]]+.*", STAMP_DF$VariantHGVSCoding) == FALSE),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries with VariantHGVSCoding field formatted incorrectly:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

out.DF <- STAMP_DF[which(!is.na(STAMP_DF$VariantHGVSGenomic) &
                           grepl("^chr[[:alnum:]]{1,}:g.[[:digit:]]+.*", STAMP_DF$VariantHGVSGenomic) == FALSE),]
if (nrow(out.DF) > 0) {
  cat("STAMP entries with VariantHGVSGenomic field formatted incorrectly:", "\n")
  print(unique(out.DF[,c("VariantLabel","VariantHGVSGenomic")]), row.names = FALSE, quote = FALSE)
  cat("\n")
}

remove(out.DF,DF_patient,STAMP_DF_structured_anno,stamp_reference_gene.list)
cat("\n")
