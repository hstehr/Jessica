sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

# Specify key table
#----------------------------------------------
site.list <- sort(unique(append(unique(STAMP_DF$PrimaryTumorSite),
                                append(unique(STAMP_CNV$PrimaryTumorSite),
                                       unique(STAMP_Fusion$PrimaryTumorSite[which(STAMP_Fusion$Gene1 %in% fusion.gene.list.full | 
                                                                                  STAMP_Fusion$Gene2 %in% fusion.gene.list.full)])))))

site.list <- data.frame(PrimaryTumorSite=site.list, stringsAsFactors = FALSE)
site.list$CohortName <- gsub("[(].*", "", site.list$PrimaryTumorSite)
site.list$CohortName <- gsub("[[:blank:]]and[[:blank:]]", " ", site.list$CohortName)

# Capitalize first letter of word in "CohortName" column
for (row_No in 1:nrow(site.list)) {
  s <- strsplit(site.list$CohortName[row_No], " ")[[1]] 
  site.list$CohortName[row_No] <- 
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
  
  remove(s)
}
# Remove white space
site.list$CohortName <- gsub("[[:blank:]]", "", site.list$CohortName)
site.list <- site.list[order(site.list$PrimaryTumorSite),]

assign("site.list.total", site.list, envir = .GlobalEnv)

# Iterate through each unique primary tumor site
#----------------------------------------------
for (site_num in 1:nrow(site.list.total)) {
  cohort_id = site.list.total$CohortName[site_num]
  site_id = site.list.total$PrimaryTumorSite[site_num]
  cat(paste(site_num,": ", cohort_id, sep=""),"\n")
  
  # Subset DFs
  #----------------------------------------------
  site_DF <- STAMP_DF[which(STAMP_DF$PrimaryTumorSite == site_id),]
  site_Fusion <- STAMP_Fusion[which(STAMP_Fusion$PrimaryTumorSite == tolower(site_id)),]
  site_CNV <- STAMP_CNV[which(STAMP_CNV$PrimaryTumorSite == tolower(site_id)),]
  
  if (nrow(site_DF) > 0) {
    pt_mutation_count_fxn (DF = site_DF, cohort = cohort_id, outdir = outdir)
  }
  
  # Fusion info per site
  #----------------------------------------------
  if (nrow(site_Fusion) > 0) {
    top_fusion_count_fxn (DF_Fusion = site_Fusion, cohort = cohort_id, outdir = outdir)
  }
  
  # Gene info per site 
  #----------------------------------------------
  top_gene_count_fxn (DF_SNVIndel = site_DF,
                      DF_Fusion = site_Fusion,
                      DF_CNV = site_CNV,
                      cohort = cohort_id, outdir = outdir)
  
  # Variant info per site 
  #----------------------------------------------
  top_variant_count_fxn (DF_SNVIndel = site_DF,
                         DF_Fusion = site_Fusion,
                         DF_CNV = site_CNV,
                         cohort = cohort_id, outdir = outdir)
  
  # Sample info per site
  #----------------------------------------------
  col_extract <- c("PatientID","PatientGender","PatientAge")
  gender_age_distribution_fxn (DF = rbind(site_DF[,col_extract],
                                          site_Fusion[,col_extract],
                                          site_CNV[!is.na(site_CNV$PatientGender),col_extract]),
                               cohort = cohort_id, outdir = outdir)
  
  col_extract <- c("PatientID","smpl.specimenType")
  specimen_type_stacked_fxn (DF = rbind(site_DF[,col_extract],
                                        site_Fusion[,col_extract],
                                        site_CNV[!is.na(site_CNV$smpl.specimenType),col_extract]),
                             cohort = cohort_id, outdir = outdir)
  
  col_extract <- c("PatientID","smpl.percentTumor")
  tumor_purity_count_fxn (DF = rbind(site_DF[,col_extract],
                                     site_Fusion[,col_extract],
                                     site_CNV[!is.na(site_CNV$smpl.percentTumor),col_extract]),
                          cohort = cohort_id,outdir = outdir)
  
  col_extract <- c("PatientID","HistologicalDx")
  histologicaldx_distribution_fxn (DF = rbind(site_DF[,col_extract],
                                              site_Fusion[,col_extract],
                                              site_CNV[!is.na(site_CNV$HistologicalDx),col_extract]),
                                   cohort = cohort_id, outdir = outdir) 
  
  col_extract <- c("PatientID","AssayDateReceived")
  test_volume_timeline_fxn (DF = rbind(site_DF[,col_extract],
                                       site_Fusion[,col_extract],
                                       site_CNV[!is.na(site_CNV$AssayDateReceived),col_extract]),
                            cohort = cohort_id, outdir = outdir, 
                            PerSite = TRUE)
  
  remove(cohort_id,site_id,col_extract,site_DF,site_Fusion,site_CNV)
}
remove(row_No,site_num,sites.addition.Fusion,site.list)
