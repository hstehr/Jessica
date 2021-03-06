sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

#################################
## STAMP database QC PIPELINE
#################################

# Syapse export of STAMP v2 (130 genes)
#----------------------------------------------
STAMP.file = paste(data.root, "STAMP/2019-04-30_syapse_export_all_variants_patientNameAndMrnRemoved.csv", sep="")
STAMP_DF <- read.csv(file = STAMP.file, header = TRUE, na.strings = c(""," ","<NA>","NA"), stringsAsFactors = FALSE, sep = ",")
Syapse_Export_timestamp <- 
  format(as.Date(gsub("([[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*$)", "\\1",
                      sub(".*/", "", STAMP.file))), format= "%Y-%m-%d")

# Merge with Syapse export of not STAMP v2 (130 genes)
#----------------------------------------------
# Merge exports from timestamp = c("2018-10-18") and timestamp = c("2019-04-30") - Filtered for "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"
STAMP_DF_old <- read.csv(file ="~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-10-18_syapse_export_all_variants_patientNameAndMrnRemoved.csv", 
                         header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")
STAMP_DF_old <- 
  STAMP_DF_old[which(STAMP_DF_old$smpl.CancerSomaticMutationReport...smpl.assayName != "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"),]

# Incorporate additional columns
DF_TumorSite <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2019-01-31_syapse_export_all.csv",
                         header = TRUE, na.strings = c("NA","None"), stringsAsFactors = FALSE,sep = ",")
# Remove extraneous columns
DF_TumorSite <- DF_TumorSite[,c("UNIQUE_ID","PRIMARY_TUMOR_SITE","HISTOLOGICAL_DIAGNOSIS")]
colnames(DF_TumorSite) <- c("smpl.TestRequest...sys.uniqueId","smpl.Patient...smpl.primaryTumorSite","smpl.Patient...smpl.histologicalDiagnosis")

# Remove duplicate rows   
DF_TumorSite <- DF_TumorSite %>% dplyr::distinct(smpl.TestRequest...sys.uniqueId, .keep_all = TRUE)

# Merge with STAMP entries 
STAMP_DF_old <- left_join(STAMP_DF_old, DF_TumorSite, by = c("smpl.TestRequest...sys.uniqueId"))

# Format dob entries
STAMP_DF$smpl.CancerSomaticMutationReport...base.dob <- 
  format(as.Date(STAMP_DF$smpl.CancerSomaticMutationReport...base.dob, format = "%Y-%m-%d"), "%m/%d/%y")

# Merge datasets 
column_common <- intersect(colnames(STAMP_DF),colnames(STAMP_DF_old))

STAMP_DF <- rbind(STAMP_DF[,column_common],STAMP_DF_old[,column_common])
remove(STAMP_DF_old,column_common,DF_TumorSite)

# Clean up patient data from Syapse
#----------------------------------------------
source(paste(pipeline.root,"SyapseExport_RetrospectiveAnalysis/Syapse_Export_QC.R",sep=""))
source(paste(pipeline.root,"SyapseExport_RetrospectiveAnalysis/Syapse_VariantAnnotate.R",sep=""))

# Filter for STAMP v2
STAMP_DF <- STAMP_DF[which(STAMP_DF$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"), ]
# sort(unique(STAMP_DF$AssayName))

# 2019-05-31 UPDATE: ignore HistologicalDx field for the time being in regards to STAMPEDE
# Do not filter for entries from adult patients ie. age >= 18yo
# Do not filter for entries with proper HGVS nomenclature = keep mutations in genes such as TERT

cat(paste("STAMP v2 POST-QC counts: ",nrow(STAMP_DF), " total entries and ", 
          length(unique(STAMP_DF[[1]])), " total test orders", sep=""),"\n","\n")

# Filter for entries with primary tumor site
#----------------------------------------------
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$PrimaryTumorSite),]
STAMP_DF <- STAMP_DF[!(STAMP_DF$PrimaryTumorSite %in% c("unknown","none","other primary site")), ]

# Collapse similar primary tumor site
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"

# Abbreviate for visualization
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
STAMP_DF$PrimaryTumorSite[which(STAMP_DF$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"
# sort(unique(STAMP_DF$PrimaryTumorSite))

# Filter entries for complete cases = percent tumor
#----------------------------------------------
STAMP_DF <- STAMP_DF[which(!is.na(STAMP_DF$smpl.percentTumor)),]
STAMP_DF <- STAMP_DF[which(grepl("comment", STAMP_DF$smpl.percentTumor, ignore.case = TRUE) == FALSE),]
STAMP_DF$smpl.percentTumor <- gsub("%$", "", STAMP_DF$smpl.percentTumor)
STAMP_DF$smpl.percentTumor <- gsub("^>", "", STAMP_DF$smpl.percentTumor)
STAMP_DF$smpl.percentTumor <- gsub("^<", "", STAMP_DF$smpl.percentTumor)
STAMP_DF <- STAMP_DF[which(!STAMP_DF$smpl.percentTumor %in% c("10-Jun","10-May","15-Oct","5-Mar","20-Oct","2-Jan")),]
STAMP_DF$smpl.percentTumor <- gsub("(^[[:digit:]]+)([-][[:digit:]]+$)", "\\1", STAMP_DF$smpl.percentTumor)
STAMP_DF$smpl.percentTumor <- ceiling(as.numeric(STAMP_DF$smpl.percentTumor))
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$smpl.percentTumor),]
# sort(unique(STAMP_DF$smpl.percentTumor))

# Filter entries for complete cases = specimen type
#----------------------------------------------
STAMP_DF <- STAMP_DF[which(STAMP_DF$smpl.specimenType != "other"),]
STAMP_DF <- STAMP_DF[complete.cases(STAMP_DF$smpl.specimenType),]
# sort(unique(STAMP_DF$smpl.specimenType))

cat(paste("Post-visualization QC-filter: ",nrow(STAMP_DF), " total entries and ", 
          length(unique(STAMP_DF[[1]])), " total test orders", sep=""),"\n","\n")

assign("STAMP_DF", STAMP_DF, envir = .GlobalEnv)

write.table(STAMP_DF, file = paste(tempdir, Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
