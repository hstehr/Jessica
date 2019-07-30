sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

# SNV/Indel info
#----------------------------------------------
pt_mutation_count_fxn (DF = STAMP_DF, cohort=cohort_id, outdir = outdir)

# Fusion info
#----------------------------------------------
fusion_count_fxn (DF_Fusion = STAMP_Fusion, 
                  cohort=cohort_id, 
                  outdir = outdir)

top_fusion_count_fxn (DF_Fusion = STAMP_Fusion, 
                      cohort=cohort_id, 
                      outdir = outdir)

# Gene info
#----------------------------------------------
gene_count_fxn (DF_SNVIndel = STAMP_DF,
                DF_Fusion = STAMP_Fusion,
                DF_CNV = STAMP_CNV,
                cohort=cohort_id, outdir = outdir)

top_gene_count_fxn (DF_SNVIndel = STAMP_DF,
                    DF_Fusion = STAMP_Fusion,
                    DF_CNV = STAMP_CNV,
                    cohort=cohort_id, outdir = outdir)

# Variant info
#----------------------------------------------
top_variant_count_fxn (DF_SNVIndel = STAMP_DF,
                       DF_Fusion = STAMP_Fusion,
                       DF_CNV = STAMP_CNV,
                       cohort=cohort_id, outdir = outdir)

variant_type_distribution_fxn (DF_SNVIndel = STAMP_DF, 
                               DF_Fusion = STAMP_Fusion,
                               DF_CNV = STAMP_CNV,
                               cohort=cohort_id, outdir = outdir)

# Sample info
#----------------------------------------------
col_extract <- c("PatientID","PrimaryTumorSite")
site_count_fxn (DF = rbind(STAMP_DF[,col_extract],
                           STAMP_Fusion[,col_extract],
                           STAMP_CNV[!is.na(STAMP_CNV$PrimaryTumorSite),col_extract]),
                cohort=cohort_id, outdir = outdir)

col_extract <- c("PatientID","PrimaryTumorSite")
top_site_count_fxn (DF = rbind(STAMP_DF[,col_extract],
                               STAMP_Fusion[,col_extract],
                               STAMP_CNV[!is.na(STAMP_CNV$PrimaryTumorSite),col_extract]),
                    cohort=cohort_id, outdir = outdir)

col_extract <- c("PatientID","PatientGender","PatientAge")
gender_age_distribution_fxn (DF = rbind(STAMP_DF[,col_extract],
                                        STAMP_Fusion[,col_extract],
                                        STAMP_CNV[!is.na(STAMP_CNV$PatientAge),col_extract]),
                             cohort=cohort_id, outdir = outdir)

col_extract <- c("PatientID","smpl.specimenType")
specimen_type_stacked_fxn (DF = rbind(STAMP_DF[,col_extract],
                                      STAMP_Fusion[,col_extract],
                                      STAMP_CNV[!is.na(STAMP_CNV$smpl.specimenType),col_extract]),
                           cohort=cohort_id, outdir = outdir)

col_extract <- c("PatientID","smpl.percentTumor")
tumor_purity_count_fxn (DF = rbind(STAMP_DF[,col_extract],
                                   STAMP_Fusion[,col_extract],
                                   STAMP_CNV[!is.na(STAMP_CNV$smpl.percentTumor),col_extract]),
                        cohort=cohort_id, outdir = outdir)

col_extract <- c("PatientID","HistologicalDx")
histologicaldx_distribution_fxn (DF = rbind(STAMP_DF[,col_extract],
                                            STAMP_Fusion[,col_extract],
                                            STAMP_CNV[!is.na(STAMP_CNV$HistologicalDx),col_extract]),
                                 cohort=cohort_id, outdir = outdir) 

col_extract <- c("PatientID","AssayDateReceived")
test_volume_timeline_fxn (DF = rbind(STAMP_DF[,col_extract],
                                     STAMP_Fusion[,col_extract],
                                     STAMP_CNV[!is.na(STAMP_CNV$AssayDateReceived),col_extract]),
                          cohort=cohort_id, outdir = outdir, 
                          PerSite = FALSE)

# Analysis of test volume based on TRF No.
#----------------------------------------------
TRF_volume_timeline_fxn (DF = TRF_DF, cohort=cohort_id, outdir = outdir)
                         
remove(col_extract,cohort_id)
