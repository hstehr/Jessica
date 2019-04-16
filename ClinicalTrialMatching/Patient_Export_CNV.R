## Output: "patient_id_CNV.tsv"

# Load relevant file from global environment 
STAMP_CNV_structured <- 
  cbind(data.frame(sys.uniqueId = patient.id, stringsAsFactors = FALSE), STAMP_CNV)

# Filter for variants == "AMP" and "DEL"
STAMP_CNV_structured <- STAMP_CNV_structured[which(STAMP_CNV_structured$Status %in% 
                                                     c("DEL","CHECK_DEL","CHECK_AMP","AMP")), ]

# Restructure STAMP dataframe
#----------------------------------------------
if (nrow(STAMP_CNV_structured) == 0) {
  
  # Generate empty dataframe
  STAMP_CNV_structured <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(STAMP_CNV_structured) <- 
    c("PatientID","CNV_Gene","Locus","Tiles","mean.z","mcopies","var.anno","var.type",
      "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus")
  
  # Output indication to file
  cat(paste(patient.id, " does not have any copy number variation entries with statuses of \"AMP\" or \"DEL\".", sep=""),"\n","\n")
  
} else if (nrow(STAMP_CNV_structured) > 0) {
  
  # Restructure STAMP dataframe
  #----------------------------------------------
  STAMP_CNV_structured$var.type <- "CNV"
  STAMP_CNV_structured$Status <- gsub("^CHECK_","",STAMP_CNV_structured$Status)
  
  # Subset columns of interest
  #----------------------------------------------
  colnames_keep <- c("sys.uniqueId","Gene","Locus","Tiles","mean.z","mcopies","Status","var.type")
  STAMP_CNV_structured <- STAMP_CNV_structured[,colnames_keep]
  
  colnames_generic <- c("PatientID","CNV_Gene","Locus","Tiles","mean.z","mcopies","var.anno","var.type")
  colnames(STAMP_CNV_structured) <- colnames_generic
  
  # Missing information
  #----------------------------------------------
  STAMP_CNV_structured$PrimaryTumorSite.Category <- "unknown"
  STAMP_CNV_structured$PrimaryTumorSite <- "unknown"
  STAMP_CNV_structured$VariantPathogenicityStatus <- "NULL"
  STAMP_CNV_structured$HistologicalDx <- "NULL"
  
  remove(colnames_generic,colnames_keep)
}

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_CNV", STAMP_CNV_structured, envir = .GlobalEnv)

## Export entries per patient into .tsv file
#----------------------------------------------
# Extract STAMP entries for individual patient
DF_patient <- STAMP_CNV[which(STAMP_CNV$PatientID == patient.id),]

# Order alphabetically by gene name
DF_patient <- DF_patient[order(DF_patient$CNV_Gene, decreasing = FALSE),]

## Write to local computer
write.table(DF_patient, file = paste(tempdir, patient.id, "_CNV.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

remove(DF_patient,STAMP_CNV_structured)
cat("\n")
