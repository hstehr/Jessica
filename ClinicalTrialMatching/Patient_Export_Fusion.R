## Output: "patient_id_Fusion.tsv"

if (nrow(STAMP_Fusion) == 0) {
  
  # Generate empty dataframe
  STAMP_Fusion_structured <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(STAMP_Fusion_structured) <- 
    c("PatientID","Fusion_Detail","Gene","Break","var.type","var.anno",
      "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus")
  
  # Output indication to file
  cat(paste(patient.id, " does not have any entries that are fusions.", sep=""),"\n","\n")
  
} else if (nrow(STAMP_Fusion) > 0) {
  
  # Load relevant file from global environment 
  STAMP_Fusion_structured <- 
    cbind(data.frame(sys.uniqueId = patient.id, stringsAsFactors = FALSE), STAMP_Fusion)
  
  # Restructure STAMP dataframe
  #----------------------------------------------
  STAMP_Fusion_structured$var.type <- "Fusion"
  STAMP_Fusion_structured$var.anno <- "Fusion"
  
  STAMP_Fusion_structured$sys.label <- paste(STAMP_Fusion_structured$Region1,STAMP_Fusion_structured$Region2,sep="-")
  
  # Subset columns of interest
  #----------------------------------------------
  colnames_keep <- c("sys.uniqueId","sys.label","Region1","Region2","Break1","Break2","var.type","var.anno")
  STAMP_Fusion_structured <- STAMP_Fusion_structured[,colnames_keep]
  
  colnames_generic <- c("PatientID","Fusion_Detail","Gene1","Gene2","Break1","Break2","var.type","var.anno")
  colnames(STAMP_Fusion_structured) <- colnames_generic
  
  # Missing information
  #----------------------------------------------
  STAMP_Fusion_structured$PrimaryTumorSite.Category <- "unknown"
  STAMP_Fusion_structured$PrimaryTumorSite <- "unknown"
  STAMP_Fusion_structured$VariantPathogenicityStatus <- "NULL"
  STAMP_Fusion_structured$HistologicalDx <- "NULL"
  
  # Reformat to long format i.e. single gene per row 
  #----------------------------------------------
  colnames_merged <- c("PatientID","Fusion_Detail","Gene","Break","var.type","var.anno",
                       "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus","HistologicalDx")
  
  Gene1.list <- data.frame(STAMP_Fusion_structured$PatientID,
                           STAMP_Fusion_structured$Fusion_Detail,
                           STAMP_Fusion_structured$Gene1,
                           STAMP_Fusion_structured$Break1,
                           STAMP_Fusion_structured$var.type,
                           STAMP_Fusion_structured$var.anno,
                           STAMP_Fusion_structured$PrimaryTumorSite.Category,
                           STAMP_Fusion_structured$PrimaryTumorSite,
                           STAMP_Fusion_structured$VariantPathogenicityStatus,
                           STAMP_Fusion_structured$HistologicalDx,
                           stringsAsFactors = FALSE)
  colnames(Gene1.list) <- colnames_merged
  
  Gene2.list <- data.frame(STAMP_Fusion_structured$PatientID,
                           STAMP_Fusion_structured$Fusion_Detail,
                           STAMP_Fusion_structured$Gene2,
                           STAMP_Fusion_structured$Break2,
                           STAMP_Fusion_structured$var.type,
                           STAMP_Fusion_structured$var.anno,
                           STAMP_Fusion_structured$PrimaryTumorSite.Category,
                           STAMP_Fusion_structured$PrimaryTumorSite,
                           STAMP_Fusion_structured$VariantPathogenicityStatus,
                           STAMP_Fusion_structured$HistologicalDx,
                           stringsAsFactors = FALSE)
  colnames(Gene2.list) <- colnames_merged
  
  STAMP_Fusion_structured <- rbind(Gene1.list,Gene2.list)
  remove(Gene1.list,Gene2.list,colnames_merged,colnames_keep,colnames_generic)
}

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_Fusion", STAMP_Fusion_structured, envir = .GlobalEnv)

## Export entries per patient into .tsv file
#----------------------------------------------
# Extract STAMP entries for individual patient
DF_patient <- STAMP_Fusion[which(STAMP_Fusion$PatientID == patient.id),]

# Order alphabetically by gene name
DF_patient <- DF_patient[order(DF_patient$Gene, decreasing = FALSE),]

## Write to local computer
write.table(DF_patient, file = paste(tempdir, patient.id, "_Fusion.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

remove(STAMP_Fusion_structured,DF_patient)
cat(paste("Timestamp of patient STAMP processing FINISH: ", Sys.time(), sep=""),"\n","\n")
