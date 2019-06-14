library(gplots)

# STAMP v2 test orders (total)
#----------------------------------------------
TRF_DF_v2 <- TRF_DF[which(TRF_DF$AssayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)"),]
cat(paste("STAMP v2 test orders: ",nrow(TRF_DF_v2),sep=""),"\n","\n")

# Generalized QC filters 
#----------------------------------------------
TRF_DF_QC <- TRF_DF_v2[complete.cases(TRF_DF_v2$PatientGender) & complete.cases(TRF_DF_v2$PatientDOB),]
cat(paste("General post-QC filter test orders: ",nrow(TRF_DF_QC),sep=""),"\n","\n")

complete.list <- sort(unique(append(STAMP_DF_QC$PatientID, 
                                    append(STAMP_Fusion_QC$PatientID, STAMP_CNV_QC$PatientID))))
cat(paste("Outside box n = ",nrow(TRF_DF_QC[!(TRF_DF_QC$PatientID %in% complete.list),]),sep=""),"\n","\n")

venn(list("SNV/Indels" = sort(unique(STAMP_DF_QC$PatientID)),
          "Fusion" = sort(unique(STAMP_Fusion_QC$PatientID)),
          "CNV" = sort(unique(STAMP_CNV_QC$PatientID))))

# Exclude missing primary tumor sites
#----------------------------------------------
TRF_DF_complete <- TRF_DF_QC[complete.cases(TRF_DF_QC$PrimaryTumorSite),]
TRF_DF_complete <- TRF_DF_complete[!(TRF_DF_complete$PrimaryTumorSite %in% c("unknown","none","other primary site")), ]
cat(paste("Exclude primary tumor site test orders: ",nrow(TRF_DF_complete),sep=""),"\n","\n")

complete.list <- sort(unique(append(STAMP_DF_complete$PatientID, 
                                    append(STAMP_Fusion_complete$PatientID, STAMP_CNV_complete$PatientID))))
cat(paste("Outside box n = ",nrow(TRF_DF_complete[!(TRF_DF_complete$PatientID %in% complete.list),]),sep=""),"\n","\n")

venn(list("SNV/Indels" = sort(unique(STAMP_DF_complete$PatientID)),
          "Fusion" = sort(unique(STAMP_Fusion_complete$PatientID)),
          "CNV" = sort(unique(STAMP_CNV_complete$PatientID))))

# AgeFilter ON
#----------------------------------------------
TRF_DF_age <- TRF_DF_complete[which(TRF_DF_complete$PatientAge >= 18),]
cat(paste("AgeFILTER i.e. only adults test orders: ",nrow(TRF_DF_age),sep=""),"\n","\n")

complete.list <- sort(unique(append(STAMP_DF_age$PatientID, 
                                    append(STAMP_Fusion_age$PatientID, STAMP_CNV_age$PatientID))))
cat(paste("Outside box n = ",nrow(TRF_DF_age[!(TRF_DF_age$PatientID %in% complete.list),]),sep=""),"\n","\n")

venn(list("SNV/Indels" = sort(unique(STAMP_DF_age$PatientID)),
          "Fusion" = sort(unique(STAMP_Fusion_age$PatientID)),
          "CNV" = sort(unique(STAMP_CNV_age$PatientID))))

# Exclude Lung cases 
#----------------------------------------------
TRF_DF_noLung <- TRF_DF_age[which(tolower(TRF_DF_age$PrimaryTumorSite) != "lung"),]
cat(paste("Exclude lung cases test orders: ",nrow(TRF_DF_noLung),sep=""),"\n","\n")

complete.list <- sort(unique(append(STAMP_DF_noLung$PatientID, 
                                    append(STAMP_Fusion_noLung$PatientID, STAMP_CNV_noLung$PatientID))))
cat(paste("Outside box n = ",nrow(TRF_DF_noLung[!(TRF_DF_noLung$PatientID %in% complete.list),]),sep=""),"\n","\n")

venn(list("SNV/Indels" = sort(unique(STAMP_DF_noLung$PatientID)),
          "Fusion" = sort(unique(STAMP_Fusion_noLung$PatientID)),
          "CNV" = sort(unique(STAMP_CNV_noLung$PatientID))))

# Lab Replication
#----------------------------------------------
TRF_DF_QC$AssayReportDateReviewed <- as.Date(TRF_DF_QC$AssayReportDateReviewed, format = "%Y-%m-%d")
TRF_DF_LabRep <- TRF_DF_QC[which(TRF_DF_QC$AssayReportDateReviewed >= "2017-07-01" &
                                     TRF_DF_QC$AssayReportDateReviewed <= "2017-12-31"),]

cat(paste("Lab replication test orders: ",nrow(TRF_DF_LabRep),sep=""),"\n","\n")

complete.list <- sort(unique(append(STAMP_DF_LabRep$PatientID, 
                                    append(STAMP_Fusion_LabRep$PatientID, STAMP_CNV_LabRep$PatientID))))
cat(paste("Outside box n = ",nrow(TRF_DF_LabRep[!(TRF_DF_LabRep$PatientID %in% complete.list),]),sep=""),"\n","\n")

venn(list("SNV/Indels" = sort(unique(STAMP_DF_LabRep$PatientID)),
          "Fusion" = sort(unique(STAMP_Fusion_LabRep$PatientID)),
          "CNV" = sort(unique(STAMP_CNV_LabRep$PatientID))))
