## Confirmed for timestamp = c("2019-02-01","2018-12-11","2018-09-06","2018-01-01","2019-02-25")

## Parse info across ARMs into separate dataframes
## Classification #1 (Variant_Type): SNV, Frameshift/In-frame (i.e. Delins, Insertions, Deletions, Duplications)
## Output: "Patient_Variant_Report_QC_.xlsx"

cat(paste("Timestamp of Patient_Variant_Report: ", Patient_Variant_Report_timestamp, sep=""),"\n","\n")

################################
## Manual edit
################################
# Removal of ARMs based on information in "ARM-GENE-LOOK-UP-TABLE"
if (isTRUE(Patient_Variant_Report_timestamp %in% c("2019-02-01","2018-12-11"))) {
  sink(file = out.ouput, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("ARM-Z1C and ARM-Z1F have been removed based on comments in ARM-GENE-LOOK-UP-TABLE", "\n","\n")
  
  sink()
  
  # names(PATIENT_VARIANT_REPORT)
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-21] # `ARM-Z1C`
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-22] # `ARM-Z1F`
  
} else if (isTRUE(Patient_Variant_Report_timestamp == "2019-02-25")) {
  sink(file = out.ouput, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("ARM-Z1C, ARM-Z1F, ARM-K1 and ARM-M have been removed based on comments in ARM-GENE-LOOK-UP-TABLE", "\n")
  
  sink()
  
  # names(PATIENT_VARIANT_REPORT)
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-21] # `ARM-Z1C`
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-22] # `ARM-Z1F`
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-14] # `ARM-K1`
  PATIENT_VARIANT_REPORT <- PATIENT_VARIANT_REPORT[-16] # `ARM-M`
}

## Specify parameters of output files 
#----------------------------------------------
DF_Inclusion_NonHotspot_Rules <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Inclusion_NonHotspot_Rules) <- c("Arm_Name", "Description", "oncominevariantclass",
                                             "Gene_Name","Exon","Function",	"Level_of_Evidence",
                                             "PRESENT", "Protein","VARIANT_TYPE","Chromosome",
                                             "position","Ref","Alt","ALLELE_FREQ","Variant_ID")

DF_Exclusion_NonHotspot_Rules <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Exclusion_NonHotspot_Rules) <- colnames(DF_Inclusion_NonHotspot_Rules)

DF_Exclusion_Variants <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Exclusion_Variants) <- c("Arm_Name", "Gene_Name","Variant_ID","Variant_Type",
                                     "Protein","Level_of_Evidence","HGVS","Chromosome",
                                     "Position","Ref","Alt","PRESENT","ALLELE_FREQ",
                                     "CN_value","FUSION_READ_DEPTH","Variant Source")

DF_Inclusion_Variants <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Inclusion_Variants) <- colnames(DF_Exclusion_Variants)

DF_IHC_Results <- data.frame(matrix(NA, ncol = 7))
colnames(DF_IHC_Results) <- c("Arm_Name","Gene","Status_POSITIVE_NEGATIVE_INDETERMINATE",
                              "Variant_PRESENT_NEGATIVE_EMPTY","Description","LOE","IHC_RESULT")

DF_Comments <- data.frame(matrix(NA, ncol = 16))
colnames(DF_Comments) <- c("Arm_Name","Gene","X__1","X__2","X__3",
                           "X__4","X__5","X__6","X__7","X__8","X__9",
                           "X__10","X__11","X__12","X__13","X__14")

DF_Histologic_Disease_Exclusion_Codes <- data.frame(matrix(NA, ncol = 6))

## Parse info across ARMs from "Disease Exclusion LOOK-UP Table"
#----------------------------------------------
for (tab_No in which(names(PATIENT_VARIANT_REPORT) == "Disease Exclusion LOOK-UP Table")) {
  Disease_Exclusion_file <- PATIENT_VARIANT_REPORT[[tab_No]]
  
  ## Remove rows that are all empty
  Disease_Exclusion_file <- 
    Disease_Exclusion_file[rowSums(is.na(Disease_Exclusion_file)) != 
                             ncol(Disease_Exclusion_file),]  
  
  row_start_list = (which(Disease_Exclusion_file[[1]] == "Histologic Disease Exclusion Codes"))
  
  for (Histo_No in 1:length(row_start_list)) {
    row_start = row_start_list[Histo_No] +2
    
    # Extract ARM_Name
    if (isTRUE(row_start_list[Histo_No] == 1)) { Arm_Name <- colnames(Disease_Exclusion_file)[2]
    } else { Arm_Name <- as.character(Disease_Exclusion_file[row_start_list[Histo_No] -1,2])
    }
    Arm_Name <- gsub("^(EAY131)(.*)", "ARM\\2", Arm_Name)
    
    if (isTRUE(Histo_No == length(row_start_list))) { row_end = as.numeric(nrow(Disease_Exclusion_file))
    } else { row_end = row_start_list[Histo_No +1] -2 
    }
    
    # Extract rows per ARM_Name
    DF_Histologic_Disease_Exclusion_pre <- Disease_Exclusion_file[c(row_start:row_end),1:5]
    DF_Histologic_Disease_Exclusion_pre$Arm_Name <- Arm_Name
    
    # Reorder columns for consistency
    DF_Histologic_Disease_Exclusion_pre <- data.frame(DF_Histologic_Disease_Exclusion_pre[,c(6,1:5)])
    
    if (isTRUE(row_start_list[Histo_No] == 1)) { DF_Histologic_Disease_Exclusion_Codes <- DF_Histologic_Disease_Exclusion_pre
    } else { DF_Histologic_Disease_Exclusion_Codes <- rbind(DF_Histologic_Disease_Exclusion_Codes,
                                                            DF_Histologic_Disease_Exclusion_pre)
    }}
  remove(DF_Histologic_Disease_Exclusion_pre, Disease_Exclusion_file,Arm_Name,Histo_No,row_end,row_start,row_start_list,tab_No)
}

if (isTRUE(Patient_Variant_Report_timestamp == "2018-01-01")) {
  colnames(DF_Histologic_Disease_Exclusion_Codes) <- c("Arm_Name","SHORT.NAME","CTEP.CATEGORY",
                                                       "CTEP.TERM","CTEP.SUBCATEGORY","MedDRA.CODE")
  
  DF_Histologic_Disease_Exclusion_Codes <- 
    DF_Histologic_Disease_Exclusion_Codes[, c("Arm_Name","CTEP.CATEGORY","CTEP.SUBCATEGORY","CTEP.TERM","SHORT.NAME","MedDRA.CODE")]
  
} else {
  colnames(DF_Histologic_Disease_Exclusion_Codes) <- c("Arm_Name","CTEP.CATEGORY","CTEP.SUBCATEGORY",
                                                       "CTEP.TERM","SHORT.NAME","MedDRA.CODE")
}

## Remove rows that are all empty
DF_Histologic_Disease_Exclusion_Codes <- 
  DF_Histologic_Disease_Exclusion_Codes[which(DF_Histologic_Disease_Exclusion_Codes$MedDRA.CODE != "None"),]

## Parse info across ARMs from original file into corresponding output files
#----------------------------------------------
if (isTRUE(Patient_Variant_Report_timestamp == "2018-01-01")) {
  arm_start = which(names(PATIENT_VARIANT_REPORT) == "GENE-STRAND Look-up table") +1
  arm_end = length(names(PATIENT_VARIANT_REPORT))
  
} else {
  arm_start = which(names(PATIENT_VARIANT_REPORT) == "Patient ID information") +1
  arm_end = which(names(PATIENT_VARIANT_REPORT) == "MOIs") -1
  
}

for (Arm_No in arm_start:arm_end) {
  DF <- PATIENT_VARIANT_REPORT[[Arm_No]]
  Arm_Name <- colnames(DF)[1]
  
  ## Remove rows that are all empty
  DF <- DF[rowSums(is.na(DF)) != ncol(DF),]  
  
  row_end = (which(DF[[2]] == "Exclusion Variants")) -1
  
  if (isTRUE(row_end > 2)) {
    if (isTRUE(which(DF[[2]] == "Exclusion Non-Hotspot Rules") > 0)) {
      row_int_end = (which(DF[[2]] == "Exclusion Non-Hotspot Rules")) -1
      row_int_start = row_int_end +3
      
      DF_Inclusion_NonHotspot_Rules_pre <- DF[c(2:row_int_end),1:16]
      colnames(DF_Inclusion_NonHotspot_Rules_pre) <- colnames(DF_Inclusion_NonHotspot_Rules) 
      DF_Inclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Inclusion_NonHotspot_Rules <- rbind(DF_Inclusion_NonHotspot_Rules, DF_Inclusion_NonHotspot_Rules_pre)
      
      DF_Exclusion_NonHotspot_Rules_pre <- DF[c(row_int_start:row_end),1:16]
      colnames(DF_Exclusion_NonHotspot_Rules_pre) <- colnames(DF_Exclusion_NonHotspot_Rules) 
      DF_Exclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Exclusion_NonHotspot_Rules <- rbind(DF_Exclusion_NonHotspot_Rules, DF_Exclusion_NonHotspot_Rules_pre)
      
    } else {
      DF_Inclusion_NonHotspot_Rules_pre <- DF[c(2:row_end),1:16]
      colnames(DF_Inclusion_NonHotspot_Rules_pre) <- colnames(DF_Inclusion_NonHotspot_Rules) 
      DF_Inclusion_NonHotspot_Rules_pre$Arm_Name <- Arm_Name
      DF_Inclusion_NonHotspot_Rules <- rbind(DF_Inclusion_NonHotspot_Rules, DF_Inclusion_NonHotspot_Rules_pre)
    }}
  
  row_start = (which(DF[[2]] == "Exclusion Variants")) +2
  row_end = (which(DF[[2]] == "Inclusion Variants")) -1
  DF_Exclusion_Variants_pre <- DF[c(row_start:row_end),1:16]
  colnames(DF_Exclusion_Variants_pre) <- colnames(DF_Exclusion_Variants) 
  DF_Exclusion_Variants_pre$Arm_Name <- Arm_Name
  DF_Exclusion_Variants <- rbind(DF_Exclusion_Variants, DF_Exclusion_Variants_pre)
  
  row_start = (which(DF[[2]] == "Inclusion Variants")) +2
  row_end = (which(DF[[2]] == "COMMENTS:")) -1
  
  if (isTRUE(which(DF[[2]] == "IHC Results") > 0)) {
    row_int_end = (which(DF[[2]] == "IHC Results")) -1
    row_int_start = row_int_end +3
    
    if (row_start <= row_int_end) {
      DF_Inclusion_Variants_pre <- DF[c(row_start:row_int_end),1:16]
      colnames(DF_Inclusion_Variants_pre) <- colnames(DF_Inclusion_Variants) 
      DF_Inclusion_Variants_pre$Arm_Name <- Arm_Name
      DF_Inclusion_Variants <- rbind(DF_Inclusion_Variants, DF_Inclusion_Variants_pre)
    }
    
    DF_IHC_Results_pre <- DF[c(row_int_start:row_end),1:7]
    colnames(DF_IHC_Results_pre) <- colnames(DF_IHC_Results) 
    DF_IHC_Results_pre$Arm_Name <- Arm_Name
    DF_IHC_Results <- rbind(DF_IHC_Results, DF_IHC_Results_pre)
    
  } else {
    DF_Inclusion_Variants_pre <- DF[c(row_start:row_end),1:16]
    colnames(DF_Inclusion_Variants_pre) <- colnames(DF_Inclusion_Variants) 
    DF_Inclusion_Variants_pre$Arm_Name <- Arm_Name
    DF_Inclusion_Variants <- rbind(DF_Inclusion_Variants, DF_Inclusion_Variants_pre)
  }
  
  row_start = (which(DF[[2]] == "COMMENTS:")) +1
  row_end = nrow(DF)
  if (row_start < row_end) {
    DF_Comments_pre <- DF[c(row_start:row_end),1:16]
    colnames(DF_Comments_pre) <- colnames(DF_Comments) 
    DF_Comments_pre$Arm_Name <- Arm_Name
    DF_Comments<- rbind(DF_Comments, DF_Comments_pre)
  }
  
  if (exists("DF_Inclusion_NonHotspot_Rules_pre")) {remove(DF_Inclusion_NonHotspot_Rules_pre)}
  if (exists("DF_Exclusion_NonHotspot_Rules_pre")) {remove(DF_Exclusion_NonHotspot_Rules_pre)}
  if (exists("DF_Exclusion_Variants_pre")) {remove(DF_Exclusion_Variants_pre)}
  if (exists("DF_Inclusion_Variants_pre")) {remove(DF_Inclusion_Variants_pre)}
  if (exists("DF_IHC_Results_pre")) {remove(DF_IHC_Results_pre)}
  if (exists("DF_Comments_pre")) {remove(DF_Comments_pre)}
  remove(Arm_Name,row_start,row_end)
}

## Remove rows that are all empty or where column 2 == c("none", "None")
DF_Inclusion_NonHotspot_Rules <- 
  DF_Inclusion_NonHotspot_Rules[rowSums(is.na(DF_Inclusion_NonHotspot_Rules)) != ncol(DF_Inclusion_NonHotspot_Rules),]

DF_Exclusion_NonHotspot_Rules <- 
  DF_Exclusion_NonHotspot_Rules[rowSums(is.na(DF_Exclusion_NonHotspot_Rules)) != 
                                  ncol(DF_Exclusion_NonHotspot_Rules),]

DF_Exclusion_Variants <- 
  DF_Exclusion_Variants[rowSums(is.na(DF_Exclusion_Variants)) != ncol(DF_Exclusion_Variants),]
DF_Exclusion_Variants <- DF_Exclusion_Variants[DF_Exclusion_Variants$Gene_Name != "None", ] 

DF_Inclusion_Variants <- 
  DF_Inclusion_Variants[rowSums(is.na(DF_Inclusion_Variants)) != ncol(DF_Inclusion_Variants),]  
DF_Inclusion_Variants <- DF_Inclusion_Variants[DF_Inclusion_Variants$Gene_Name != "None", ] 

DF_IHC_Results <- DF_IHC_Results[rowSums(is.na(DF_IHC_Results)) != ncol(DF_IHC_Results),]
DF_IHC_Results <- DF_IHC_Results[DF_IHC_Results$Gene != "none", ]  

DF_Comments <- DF_Comments[rowSums(is.na(DF_Comments)) != ncol(DF_Comments),]

################################
## Manual edit = DF_Comments
################################
if (isTRUE(Patient_Variant_Report_timestamp %in% c("2019-02-01","2018-12-11","2018-09-06","2019-02-25"))) {
  DF_Comments$Note[max(which(DF_Comments$Arm_Name == "ARM-T"))] <- 
    paste(DF_Comments$Gene[1], DF_Comments$Gene[2], DF_Comments$Gene[3], sep=" ")
  DF_Comments <- DF_Comments[DF_Comments$Arm_Name == "ARM-T" & 
                               !is.na(DF_Comments$Note),]
  DF_Comments <- DF_Comments[,c(1:11,ncol(DF_Comments))]
  colnames(DF_Comments)[1:11] <- colnames(DF_Exclusion_Variants)[1:11]
}

################################
## Manual edit = DF_Inclusion_Variants
################################
# Remove blank space after "p."
for (row_No in 1:nrow(DF_Inclusion_Variants)) {
  if (isTRUE(DF_Inclusion_Variants$Variant_Type[row_No] == "SNV")) {
    DF_Inclusion_Variants$Protein[row_No] <- gsub("(^p.[[:blank:]]*)(.*)", "\\2", DF_Inclusion_Variants$Protein[row_No])
    DF_Inclusion_Variants$Protein[row_No] <- paste("p.", DF_Inclusion_Variants$Protein[row_No], sep="")
  }}

# Spelling correction
DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Variant_Type == "Del")] <- "Indel"
DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Variant_Type == "MNV")] <- "SNV"

if (Patient_Variant_Report_timestamp  == "2018-09-06") {
  # Spelling correction
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Protein == "p.G719C")] <- "p.Gly719Cys"
  
  # Editing of "MET Exon 14 skipping"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082G>C")] <- "p.Asp1028His"
  # Input coding info for intronic mutations
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.2:c.3082+1G>T")] <- "c.3082+1G>T"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Variant_ID == "MVAR27")] <- "c.3082_3082+26del"
  
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Protein == "p.G983_G1054del")] <- "p.Gly983_Gly1054del"
  
  # NM_000368.4(TSC1):c.1907_1908del (p.Glu636Glyfs)
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Protein == "p.Glu636fs51")] <- "p.Glu636Glyfs"
  
  # NM_000142.4(FGFR3):c.1949A>C (p.Lys650Thr)
  DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Protein == "p.Lys650Thr")] <- "SNV"
  
} else if (Patient_Variant_Report_timestamp  == "2018-12-11") {
  # Editing of "MET Exon 14 skipping"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082G>C")] <- "p.Asp1028His"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.2:c.3082+1G>T")] <- "c.3082+1G>T"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082+1G>A")] <- "c.3082+1G>A"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_000245.2:c.2888-18_2888-2del")] <- "c.2888-18_2888-2del"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Variant_ID == "MVAR27")] <- "c.3082_3082+26del"
  
  # NM_000142.4(FGFR3):c.1949A>C (p.Lys650Thr)
  DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Protein == "p.Lys650Thr")] <- "SNV"
  
  # NM_000368.4(TSC1):c.1907_1908del (p.Glu636Glyfs)
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Protein == "p.Glu636fs51")] <- "p.Glu636Glyfs"
  
} else if (Patient_Variant_Report_timestamp  == "2019-02-01") {
  # Editing of "MET Exon 14 skipping"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082G>C")] <- "p.Asp1028His"
  # Input coding info for intronic mutations
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Variant_ID == "MVAR27")] <- "c.3082_3082+26del"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082+1G>A")] <- "c.3082+1G>A"
  
  # NM_000368.4(TSC1):c.1907_1908del (p.Glu636Glyfs)
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Protein == "p.Glu636fs51")] <- "p.Glu636Glyfs"
  
  # NM_000142.4(FGFR3):c.1949A>C (p.Lys650Thr)
  DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Protein == "p.Lys650Thr")] <- "SNV"
  
} else if (Patient_Variant_Report_timestamp  == "2018-01-01") {
  row.change <- which(grepl("Va", DF_Inclusion_Variants$Protein) == TRUE & DF_Inclusion_Variants$Variant_Type == "Indel")
  DF_Inclusion_Variants$Protein[row.change] <- gsub("Va", "VA", DF_Inclusion_Variants$Protein[row.change])
  
  row.change <- which(grepl("Va", DF_Exclusion_Variants$Protein) == TRUE & DF_Exclusion_Variants$Variant_Type == "Indel")
  DF_Exclusion_Variants$Protein[row.change] <- gsub("Va", "VA", DF_Exclusion_Variants$Protein[row.change])
  
} else if (Patient_Variant_Report_timestamp  == "2019-02-25") {
  # Editing of "MET Exon 14 skipping"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082G>C")] <- "p.Asp1028His"
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$HGVS == "NM_001127500.1:c.3082+1G>A")] <- "c.3082+1G>A"
  # Input coding info for intronic mutations
  DF_Inclusion_Variants$Protein[which(DF_Inclusion_Variants$Variant_ID == "MVAR27")] <- "c.3082_3082+26del"
  
  # NM_000142.4(FGFR3):c.1949A>C (p.Lys650Thr)
  DF_Inclusion_Variants$Variant_Type[which(DF_Inclusion_Variants$Protein == "p.Lys650Thr")] <- "SNV"
}

################################
## Manual edit = DF_Exclusion_Variants
################################
# Remove blank space after "p." and/or include "p."
for (row_No in 1:nrow(DF_Exclusion_Variants)) {
  if (DF_Exclusion_Variants$Variant_Type[row_No] == "SNV") {
    DF_Exclusion_Variants$Protein[row_No] <- gsub("(^p.[[:blank:]]*)(.*)", "\\2", DF_Exclusion_Variants$Protein[row_No])
    DF_Exclusion_Variants$Protein[row_No] <- paste("p.", DF_Exclusion_Variants$Protein[row_No], sep="")
  }}

# Spelling correction
DF_Exclusion_Variants$Variant_Type[which(DF_Exclusion_Variants$Variant_Type == "Del")] <- "Indel"
DF_Exclusion_Variants$Variant_Type[which(DF_Exclusion_Variants$Variant_Type == "MNV")] <- "SNV"

################################
## Manual edit = DF_Inclusion_Variants & DF_Exclusion_Variants
################################
if (isTRUE(Patient_Variant_Report_timestamp == "2018-01-01")) {
  ## Convert codon nomenclature from 1-letter to 3-letters 
  for (row_No in 1:nrow(DF_Inclusion_Variants)) {
    if (isTRUE(DF_Inclusion_Variants$Variant_Type[row_No] %in% c("SNV","Indel"))) {
      chars <- unlist(strsplit(DF_Inclusion_Variants$Protein[row_No], ""))
      for (elem_No in 1:length(chars)) {
        if (isTRUE(grepl("[[:upper:]]", chars[elem_No]))) {
          chars[elem_No] <- AminoAcid_Conversion$Code3[which(AminoAcid_Conversion$Code1 == chars[elem_No])]
        }
      }
      DF_Inclusion_Variants$Protein[row_No] <-  paste0(chars, collapse = "")
    }
  }
  
  for (row_No in 1:nrow(DF_Exclusion_Variants)) {
    if (isTRUE(DF_Exclusion_Variants$Variant_Type[row_No] %in% c("SNV","Indel"))) {
      chars <- unlist(strsplit(DF_Exclusion_Variants$Protein[row_No], ""))
      for (elem_No in 1:length(chars)) {
        if (isTRUE(grepl("[[:upper:]]", chars[elem_No]))) {
          chars[elem_No] <- AminoAcid_Conversion$Code3[which(AminoAcid_Conversion$Code1 == chars[elem_No])]
        }
      }
      DF_Exclusion_Variants$Protein[row_No] <-  paste0(chars, collapse = "")
    }
  }
  
  # Convert asterisk to "Ter"
  DF_Inclusion_Variants$Protein <- gsub("\\*", "Ter", DF_Inclusion_Variants$Protein)
  DF_Exclusion_Variants$Protein <- gsub("\\*", "Ter", DF_Exclusion_Variants$Protein)
  
  remove(chars,elem_No)
}

# Annotate Variant_Type for consistency
for (row_No in 1:nrow(DF_Inclusion_Variants)) {
  if (isTRUE(DF_Inclusion_Variants$Variant_Type[row_No] == "Indel")) {
    if (grepl("fs", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Frameshift"
    } else if (grepl("del.*ins", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Delins"
    } else if (grepl("ins", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Insertion"
    } else if (grepl("del", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Deletion"
    } else if (grepl("dup", DF_Inclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Inclusion_Variants$Variant_Type[row_No] <- "Duplication"
    }}}

for (row_No in 1:nrow(DF_Exclusion_Variants)) {
  if (isTRUE(DF_Exclusion_Variants$Variant_Type[row_No] == "Indel")) {
    if (grepl("fs", DF_Exclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Exclusion_Variants$Variant_Type[row_No] <- "Frameshift"
    } else if (grepl("del.*ins", DF_Exclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Exclusion_Variants$Variant_Type[row_No] <- "Delins"
    } else if (grepl("ins", DF_Exclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Exclusion_Variants$Variant_Type[row_No] <- "Insertion"
    } else if (grepl("del", DF_Exclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Exclusion_Variants$Variant_Type[row_No] <- "Deletion"
    } else if (grepl("dup", DF_Exclusion_Variants$Protein[row_No]) == TRUE) {
      DF_Exclusion_Variants$Variant_Type[row_No] <- "Duplication"
    }}}

## Write variables in global environment
#----------------------------------------------
assign("Inclusion_Variants", DF_Inclusion_Variants, envir = .GlobalEnv)
assign("Exclusion_Variants", DF_Exclusion_Variants, envir = .GlobalEnv)
assign("Inclusion_NonHotspot_Rules", DF_Inclusion_NonHotspot_Rules, envir = .GlobalEnv)
assign("Exclusion_NonHotspot_Rules", DF_Exclusion_NonHotspot_Rules, envir = .GlobalEnv)
assign("IHC_Results", DF_IHC_Results, envir = .GlobalEnv)
assign("Comments", DF_Comments, envir = .GlobalEnv)
assign("Disease_Exclusion_Codes", DF_Histologic_Disease_Exclusion_Codes, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
list_of_datasets <- list("Inclusion_Variants" = DF_Inclusion_Variants, 
                         "Exclusion_Variants" = DF_Exclusion_Variants,
                         "Inclusion_NonHotspot_Rules" = DF_Inclusion_NonHotspot_Rules,
                         "Exclusion_NonHotspot_Rules" = DF_Exclusion_NonHotspot_Rules,
                         "IHC_Results" = DF_IHC_Results,
                         "Comments" = DF_Comments,
                         "Disease_Exclusion_Codes" = DF_Histologic_Disease_Exclusion_Codes)
write.xlsx(list_of_datasets, file = paste(tempdir, Patient_Variant_Report_timestamp,
                                          "_Patient_Variant_Report_QC.xlsx", sep=""))

remove(DF,list_of_datasets,arm_end,Arm_No,arm_start,row_int_end,row_int_start,row_No,
       DF_Inclusion_Variants, DF_Exclusion_Variants, DF_Inclusion_NonHotspot_Rules,
       DF_Exclusion_NonHotspot_Rules, DF_IHC_Results, DF_Comments,
       DF_Histologic_Disease_Exclusion_Codes)

out.DF <- Inclusion_Variants[which(!is.na(Inclusion_Variants$Protein) &
                                     grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*", Inclusion_Variants$Protein) == FALSE &
                                     Inclusion_Variants$Variant_Type != "Fusion" &
                                     Inclusion_Variants$Variant_Type != "CNV"),]
if (nrow(out.DF) > 0) {
  sink(file = out.ouput, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Inclusion_Variants", "\n")
  cat("Clinical trials that cannot be matched using current pipeline due to incorrectly formatted HGVS protein nomenclature:", "\n")
  out.DF$HGVSGenomic <- paste("g.",out.DF$Position,out.DF$Ref,">",out.DF$Alt, sep="")
  out.DF$Protein <- gsub("^p.","",out.DF$Protein)
  out.DF <- out.DF[, c("Arm_Name","Gene_Name","Protein","HGVSGenomic")]
  print(out.DF , quote = FALSE, row.names = FALSE)
  cat("\n")
  
  sink()
}

out.DF <- Exclusion_Variants[which(!is.na(Exclusion_Variants$Protein) &
                                     grepl("^p.[[:alpha:]]{3}[[:digit:]]+.*", Exclusion_Variants$Protein) == FALSE &
                                     Exclusion_Variants$Variant_Type != "Fusion" &
                                     Exclusion_Variants$Variant_Type != "CNV"),]
if (nrow(out.DF) > 0) {
  sink(file = out.ouput, append = TRUE, split = FALSE)
  options(max.print=999999)
  
  cat("Exclusion_Variants", "\n")
  cat("Clinical trials that cannot be matching using current pipeline due to HGVS protein nomenclature formatted incorrectly:", "\n")
  out.DF$HGVSGenomic <- paste("g.",out.DF$Position,out.DF$Ref,">",out.DF$Alt, sep="")
  out.DF <- out.DF[, c("Arm_Name","Gene_Name","Protein","HGVS","HGVSGenomic")]
  print(out.DF , quote = FALSE, row.names = FALSE)
  cat("\n")
  
  sink()
}

remove(out.DF)
