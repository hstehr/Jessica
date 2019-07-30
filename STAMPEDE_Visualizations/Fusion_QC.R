sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP Fusion database QC PIPELINE
#################################
source(paste(pipeline.root,"SyapseExport_RetrospectiveAnalysis/Fusion_Export_QC.R",sep=""))

# Remove melting of gene column
col_keep <-  c("PatientID","PatientGender","PatientDOB","PatientAge",
               "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
               "AssayName","sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
               "smpl.amendedString","smpl.amendmentReason","Fusion_Detail")
STAMP_Fusion <- unique(STAMP_Fusion[,col_keep])

# Parse Fusion_Detail
for (row_No in 1:nrow(STAMP_Fusion)) {
  gene_01 <- gsub("(^[[:alnum:]]+)([-].*)","\\1",STAMP_Fusion$Fusion_Detail[row_No])
  gene_02 <- gsub("(^[[:alnum:]]+[-])(.*)","\\2",STAMP_Fusion$Fusion_Detail[row_No])
  
  genes_order <- sort(append(gene_01,gene_02))
  
  STAMP_Fusion$Gene1[row_No] <- genes_order[[1]]
  STAMP_Fusion$Gene2[row_No] <- genes_order[[2]]
}

# Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_Fusion$AssayName))

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# Do not filter for entries from adult patients ie. age >= 18yo

cat(paste("Fusion STAMP v2 POST-QC counts: ",nrow(STAMP_Fusion), " total entries and ", 
          length(unique(STAMP_Fusion[[1]])), " total test orders", sep=""),"\n","\n")

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[complete.cases(STAMP_Fusion$PrimaryTumorSite),]
STAMP_Fusion <- STAMP_Fusion[!(STAMP_Fusion$PrimaryTumorSite %in% c("unknown","none","other primary site")), ]

# Collapse similar primary tumor site
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"
STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"

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

cat(paste("Fusion post-visualization QC-filter: ",nrow(STAMP_Fusion), " total entries and ", 
          length(unique(STAMP_Fusion[[1]])), " total test orders", sep=""),"\n","\n")

# Examine number of missing fields
#----------------------------------------------
Fusion.diff_No <- length(which(is.na(STAMP_Fusion$PrimaryTumorSite)))
cat(paste("Number of Fusion entries missing PrimaryTumorSite: ",Fusion.diff_No,sep=""),"\n","\n")

# Identify elements not in STAMP v2 annotation file
fusion.missing.list <- STAMP_Fusion$Fusion_Detail[!(STAMP_Fusion$Gene1 %in% fusion.gene.list.full |
                                                      STAMP_Fusion$Gene2 %in% fusion.gene.list.full)]
cat(paste("Number of fusion entries without corresponding Gene in STAMP v2 file '2016-08-23_STAMP2_regions.xlsx': ",
          length(fusion.missing.list), sep=""),"\n",
    paste(unlist(fusion.missing.list),collapse=", "),"\n")

remove(Fusion.diff_No,fusion.missing.list)

# Exclude non-listed fusion genes
# fusion.gene.list.full = only genes listed in '2016-08-23_STAMP2_regions.xlsx'
#----------------------------------------------
STAMP_Fusion <- STAMP_Fusion[which(STAMP_Fusion$Gene1 %in% fusion.gene.list.full |
                                     STAMP_Fusion$Gene2 %in% fusion.gene.list.full),]

assign("STAMP_Fusion", STAMP_Fusion, envir = .GlobalEnv)

write.table(STAMP_Fusion, file = paste(tempdir, Fusion_Export_timestamp, "_Fusion_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
