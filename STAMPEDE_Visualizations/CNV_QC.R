sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP CNV database QC PIPELINE
#################################
## All CNVs marked as "AMP" (i.e. only high-confidence amplifications) for each TRF id from STAMP2. 
## If the TRF is not listed in the file, 
# then either there were no amplifications, or the sample could not be matched for some reason.  
## Existence of sites.addition.CNV not addressed 

# CNV export for STAMP v2 (130 genes)
#----------------------------------------------
source(paste(pipeline.root,"SyapseExport_RetrospectiveAnalysis/CNV_Export_QC.R",sep=""))

STAMP_CNV_missing <- STAMP_CNV[is.na(STAMP_CNV$PrimaryTumorSite),]
STAMP_CNV_missing <- left_join(STAMP_CNV_missing[,1:4],
                               unique(TRF_DF[,c("PatientID","PatientGender","PatientAge","PrimaryTumorSite",
                                                  "AssayDateReceived","smpl.specimenType","smpl.percentTumor","HistologicalDx")]),
                               by = "PatientID")
CNV.diff_No <- length(which(is.na(STAMP_CNV_missing$PrimaryTumorSite)))
cat(paste("Number of CNV entries missing PrimaryTumorSite: ",CNV.diff_No,sep=""))

# Filter for STAMP v2
STAMP_CNV <- STAMP_CNV[which(STAMP_CNV$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_CNV$AssayName))

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# Do not filter for entries from adult patients ie. age >= 18yo

cat(paste("CNV STAMP v2 POST-QC counts: ",nrow(STAMP_CNV), " total entries and ", 
          length(unique(STAMP_CNV[[1]])), " total test orders", sep=""),"\n","\n")

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$PrimaryTumorSite),]
STAMP_CNV <- STAMP_CNV[!(STAMP_CNV$PrimaryTumorSite %in% c("unknown","none","other primary site")),]

# Collapse similar primary tumor site
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"

# Abbreviate for visualization
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
STAMP_CNV$PrimaryTumorSite[which(STAMP_CNV$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"
# sort(unique(STAMP_CNV$PrimaryTumorSite))

# Filter entries for complete cases = percent tumor
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[which(!is.na(STAMP_CNV$smpl.percentTumor)),]
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$smpl.percentTumor),]
STAMP_CNV$smpl.percentTumor <- gsub("^>", "", STAMP_CNV$smpl.percentTumor)
STAMP_CNV$smpl.percentTumor <- gsub("^<", "", STAMP_CNV$smpl.percentTumor)
STAMP_CNV$smpl.percentTumor <- gsub("%$", "", STAMP_CNV$smpl.percentTumor)
STAMP_CNV$smpl.percentTumor <- gsub("(^[[:digit:]]+)([-][[:digit:]]+$)", "\\1", STAMP_CNV$smpl.percentTumor)
STAMP_CNV$smpl.percentTumor <- ceiling(as.numeric(STAMP_CNV$smpl.percentTumor))
# sort(unique(STAMP_CNV$smpl.percentTumor))

# Filter entries for complete cases = specimen type
#----------------------------------------------
STAMP_CNV <- STAMP_CNV[which(STAMP_CNV$smpl.specimenType != "other"),]
STAMP_CNV <- STAMP_CNV[complete.cases(STAMP_CNV$smpl.specimenType),]
# sort(unique(STAMP_CNV$smpl.specimenType))

cat(paste("CNV post-visualization QC-filter: ",nrow(STAMP_CNV), " total entries and ", length(unique(STAMP_CNV[[1]])), " total test orders", sep=""),"\n","\n")

# Examine number of missing fields
#----------------------------------------------
CNV.list <- sort(unique(STAMP_CNV$PatientID))
CNV.diff_No <- length(CNV.list[!(CNV.list %in% sort(unique(TRF_DF$PatientID)))])
cat(paste("Number of CNV entries without corresponding PatientID in POST-Filter TRF_DF: ",
          CNV.diff_No,sep=""),"\n","\n")

CNV.diff_No <- nrow(unique(STAMP_CNV[is.na(STAMP_CNV$PrimaryTumorSite),c("PatientID","PrimaryTumorSite")]))
CNV.diff_No <- length(which(is.na(STAMP_CNV$PrimaryTumorSite)))
cat(paste("Number of CNV entries missing PrimaryTumorSite: ",CNV.diff_No,sep=""),"\n","\n")

# Identify elements not in STAMP v2 annotation file
CNV.list <- sort(unique(STAMP_CNV$CNV_Gene))
CNV.diff_No <- length(CNV.list[!(CNV.list %in% cnv.gene.list.full)])
cat(paste("Number of CNV entries without corresponding Gene in STAMP v2 file '2016-08-23_STAMP2_regions.xlsx': ",
          CNV.diff_No, sep=""),"\n",
    paste(unlist(CNV.list[!(CNV.list %in% cnv.gene.list.full)]),collapse=", "),"\n")

remove(CNV.list,CNV.diff_No)

assign("STAMP_CNV", STAMP_CNV, envir = .GlobalEnv)

write.table(STAMP_CNV, file = paste(tempdir, CNV_Export_timestamp, "_CNV_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
