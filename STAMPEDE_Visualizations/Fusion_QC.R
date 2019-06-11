sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP Fusion database QC PIPELINE
#################################
source(paste(pipeline.root,"SyapseExport_RetrospectiveAnalysis/Fusion_Export_QC.R",sep=""))

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# sort(unique(STAMP_Fusion$HistologicalDx))

# Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_Fusion$AssayName))

# # Filter for adults
# #----------------------------------------------
# STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$PatientAge >= 18),]
# sort(as.numeric(unique(STAMP_Fusion$PatientAge)))

cat(paste("Fusion POST-QC counts: ",nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total test orders", sep=""),"\n","\n")

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$PrimaryTumorSite),]
STAMP_Fusion <- STAMP_Fusion[!(STAMP_Fusion$PrimaryTumorSite %in% c("unknown","none","other primary site")), ]

# Collapse similar primary tumor site
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"

# Abbreviate for STAMPEDE display
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"

STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$PrimaryTumorSite), ]
# sort(unique(STAMP_Fusion$PrimaryTumorSite))

# Filter entries for complete cases = percent tumor
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(!is.na(STAMP_Fusion$smpl.percentTumor)),]
STAMP_Fusion <- STAMP_Fusion[which(grepl("comment", STAMP_Fusion$smpl.percentTumor, ignore.case = TRUE) == FALSE),]
STAMP_Fusion$smpl.percentTumor <- gsub("%$", "", STAMP_Fusion$smpl.percentTumor)
STAMP_Fusion$smpl.percentTumor <- gsub("^>", "", STAMP_Fusion$smpl.percentTumor)
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$smpl.percentTumor != "None"), ]

STAMP_Fusion$smpl.percentTumor <- ceiling(as.numeric(STAMP_Fusion$smpl.percentTumor))
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.percentTumor),]
# sort(unique(STAMP_Fusion$smpl.percentTumor))

# Filter entries for complete cases = specimen type
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$smpl.specimenType != "other"),]
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$smpl.specimenType),]
# sort(unique(STAMP_Fusion$smpl.specimenType))

cat(paste("Fusion post-visualization QC-filter: ",nrow(STAMP_DF), " total entries and ", length(unique(STAMP_DF[[1]])), " total test orders", sep=""),"\n","\n")

# Examine number of missing fields
#----------------------------------------------
Fusion.list <- sort(unique(STAMP_Fusion$PatientID))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PatientID)))])
cat(paste("Number of Fusion entries without corresponding PatientID in POST-Filter STAMP_DF: ",
          Fusion.diff_No,sep=""),"\n","\n")

Fusion.diff_No <- length(which(is.na(STAMP_Fusion$PrimaryTumorSite)))
cat(paste("Number of Fusion entries missing PrimaryTumorSite: ",Fusion.diff_No,sep=""),"\n","\n")

Fusion.list <- sort(unique(STAMP_Fusion$PrimaryTumorSite))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))])
cat(paste("Number of Fusion entries without corresponding PrimaryTumorSite in POST-Filter STAMP_DF: ",
          Fusion.diff_No,sep=""),"\n",
    paste(unlist(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]),collapse=", "),"\n")

# Append elements not found in SNV/Indel DF from Fusion DF
Fusion.list <- Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$PrimaryTumorSite)))]
if (length(Fusion.list) > 0) {
  sites.addition.Fusion = Fusion.list
} else {
  sites.addition.Fusion = NULL
}
assign("sites.addition.Fusion", sites.addition.Fusion, envir = .GlobalEnv)


Fusion.list <- sort(unique(append(STAMP_Fusion$Gene1,STAMP_Fusion$Gene2)))
Fusion.diff_No <- length(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))])
Fusion.list <- Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))]
cat(paste("Number of Fusion entries without corresponding Variant Gene in POST-Filter STAMP_DF: ",
          Fusion.diff_No, sep=""),"\n",
    paste(unlist(Fusion.list[!(Fusion.list %in% sort(unique(STAMP_DF$VariantGene)))]),collapse=", "),"\n","\n")

# Append elements not found in SNV/Indel DF from Fusion DF
if (length(Fusion.list) > 0) {
  genes.addition.Fusion = Fusion.list
} else {
  genes.addition.Fusion = NULL
}

remove(Fusion.list,Fusion.diff_No,Fusion.file)

cat(paste("POST-Filter counts: ",nrow(STAMP_Fusion), " total entries and ", 
          length(unique(STAMP_Fusion[[1]])), " total test orders", sep=""),"\n","\n")

assign("STAMP_Fusion", STAMP_Fusion, envir = .GlobalEnv)

write.table(STAMP_Fusion, file = paste(tempdir, Fusion_Export_timestamp, "_Fusion_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
