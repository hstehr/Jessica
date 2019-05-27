sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP CNV database QC PIPELINE
#################################
## All CNVs marked as "AMP" (i.e. only high-confidence amplifications) for each TRF id from STAMP2. 
## If the TRF is not listed in the file, 
# then either there were no amplifications, or the sample could not be matched for some reason.  
## Existence of sites.addition.CNV not addressed 

# Fusion export for STAMP v2 (130 genes)
#----------------------------------------------
CNV.file="~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-05-15_stamp_cnvs.txt"
STAMP_CNV <- read.csv(file = CNV.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
CNV_Export_timestamp <- format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",sub(".*/", "", CNV.file))), format= "%Y-%m-%d")

# Rename column names of DF
colnames(STAMP_CNV) <- c("PatientID","report_version","CNV_Gene","Variant_Type")

# Integrate with primary tumor site field
#----------------------------------------------
STAMP_CNV <- left_join(STAMP_CNV,
                       unique(STAMP_DF[,c("PatientID","PatientGender","PatientAge","PrimaryTumorSite",
                                          "AssayDateReceived","smpl.specimenType","smpl.percentTumor","HistologicalDx")]),
                       by = "PatientID")

# STAMP_CNV_missing <- STAMP_CNV[is.na(STAMP_CNV$PrimaryTumorSite),]
# STAMP_CNV_missing <- left_join(STAMP_CNV_missing[,1:4],
#                                unique(STAMP_DF[,c("PatientID","PatientGender","PatientAge","PrimaryTumorSite",
#                                                   "AssayDateReceived","smpl.specimenType","smpl.percentTumor","HistologicalDx")]),
#                                by = "PatientID")
# CNV.diff_No <- length(which(is.na(STAMP_CNV_missing$PrimaryTumorSite)))
# cat(paste("Number of CNV entries missing age: ",CNV.diff_No,sep=""))

# Filter for entries with histological dx
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$HistologicalDx), ]
# sort(unique(STAMP_CNV$HistologicalDx))

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$PrimaryTumorSite),]
STAMP_CNV <- STAMP_CNV[!(STAMP_CNV$PrimaryTumorSite %in% c("unknown","none","other primary site")),]

# Collapse similar primary tumor site
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"

# Abbreviate for visualization
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"
# sort(unique(STAMP_CNV$PrimaryTumorSite))

# Filter entries for complete cases = percent tumor
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[which(!is.na(STAMP_CNV$smpl.percentTumor)),]
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$smpl.percentTumor),]
# sort(unique(STAMP_CNV$smpl.percentTumor))

# Filter entries for complete cases = specimen type
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[which(STAMP_CNV$smpl.specimenType != "other"),]
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$smpl.specimenType),]
# sort(unique(STAMP_CNV$smpl.specimenType))

# Collapse similar gene names
#----------------------------------------------
STAMP_CNV$CNV_Gene[which(tolower(STAMP_CNV$CNV_Gene) %in% c("nkx2", "nkx2-1"))] <- "NKX2-1"
# sort(unique(STAMP_CNV$CNV_Gene))

# # Filter for adults
# #----------------------------------------------
# STAMP_CNV <- STAMP_CNV[which(STAMP_CNV$PatientAge >= 18),]
# sort(as.numeric(unique(STAMP_CNV$PatientAge)))

# Examine number of missing fields
#----------------------------------------------
CNV.list <- sort(unique(STAMP_CNV$PatientID))
CNV.diff_No <- length(CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$PatientID)))])
cat(paste("Number of CNV entries without corresponding PatientID in POST-Filter STAMP_DF: ",
          CNV.diff_No,sep=""),"\n","\n")

CNV.diff_No <- nrow(unique(STAMP_CNV[is.na(STAMP_CNV$PrimaryTumorSite),c("PatientID","PrimaryTumorSite")]))
CNV.diff_No <- length(which(is.na(STAMP_CNV$PrimaryTumorSite)))
cat(paste("Number of CNV entries missing PrimaryTumorSite: ",CNV.diff_No,sep=""),"\n","\n")

CNV.list <- sort(unique(STAMP_CNV$PrimaryTumorSite))
CNV.diff_No <- length(CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))])
cat(paste("Number of CNV entries without corresponding PrimaryTumorSite in POST-Filter STAMP_DF: ",
          CNV.diff_No,sep=""),"\n",
    paste(unlist(CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]),collapse=", "),"\n")

# Append elements not found in SNV/Indel DF from CNV DF
CNV.list <- CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]
if (length(CNV.list) > 0) {
  sites.addition.CNV = CNV.list
} else {
  sites.addition.CNV = NULL
}

CNV.list <- sort(unique(STAMP_CNV$CNV_Gene))
CNV.diff_No <- length(CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$VariantGene)))])
cat(paste("Number of CNV entries without corresponding Variant Gene in POST-Filter STAMP_DF: ",
          CNV.diff_No, sep=""),"\n",
    paste(unlist(CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$VariantGene)))]),collapse=", "),"\n","\n")

# Append elements not found in SNV/Indel DF from Fusion DF
CNV.list <- CNV.list[!(CNV.list %in% sort(unique(STAMP_DF$VariantGene)))]
if (length(CNV.list) > 0) {
  genes.addition.CNV = CNV.list
} else {
  genes.addition.CNV = NULL
}

remove(CNV.list,CNV.diff_No,CNV.file)

cat(paste("POST-Filter counts: ",nrow(STAMP_CNV), " total entries and ", length(unique(STAMP_CNV[[1]])), " total test orders", sep=""),"\n","\n")

assign("STAMP_CNV", STAMP_CNV, envir = .GlobalEnv)
write.table(STAMP_CNV, file = paste(tempdir, CNV_Export_timestamp, "_CNV_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
