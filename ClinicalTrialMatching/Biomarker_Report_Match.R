# Match by gene > biomarker.condition > biomarker.detail > pathogenicity.status \
# > age.group (i.e. adult) > disease.group.category > disease.site
## Search Output: individual patient file with stop point
## Match Output: "OnCore_Biomarker_Matched.tsv"

cat(paste("Timestamp of Internal clinical trial matching START: ", Sys.time(), sep=""),"\n")

ncol_STAMP <- as.numeric(ncol(STAMP_DF))
ncol_OnCore <- as.numeric(ncol(OnCore_Biomarker_Report))

# Extract gene from OnCore_Biomarker_QC
genes.OnCore_Biomarker <- sort(unique(OnCore_Biomarker_QC$Biomarker_GeneName))

# Output file of matched mutations with clinical trials
DF_Output_OnCore_Biomarker <- data.frame(matrix(NA, ncol = (ncol_STAMP + ncol_OnCore)))
colnames(DF_Output_OnCore_Biomarker) <- append(colnames(STAMP_DF), colnames(OnCore_Biomarker_Report))

## Iterate through each patient_id of patient.list
#----------------------------------------------
for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(outdir_patient, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  # Extract unique genes per patient
  gene.patient <- sort(unique(DF_patient$VariantGene))
  
  # Clinical match result
  DF_patient$OnCore_Report_Status <- NA
  
  ## Iterate through each gene_id of gene.patient
  #----------------------------------------------
  for (gene_num in 1:length(gene.patient)) {
    gene_id <- gene.patient[gene_num]
    
    ## Match Biomarker_GeneName
    #----------------------------------------------
    if (gene_id %in% genes.OnCore_Biomarker) {
      
      # Extract gene_id match from OnCore_Biomarker_Report
      DF_Gene_OnCore_Biomarker <-
        OnCore_Biomarker_QC[which(OnCore_Biomarker_QC$Biomarker_GeneName == gene_id),]
      DF_Gene_OnCore_Biomarker$Disease.Group.category <- as.character(DF_Gene_OnCore_Biomarker$Disease.Group.category)
      DF_Gene_OnCore_Biomarker$Disease.Site <- as.character(DF_Gene_OnCore_Biomarker$Disease.Site)
      
      ## Parse Biomarker_Condition
      #----------------------------------------------
      biomarker.condition.trial <- unique(DF_Gene_OnCore_Biomarker$Biomarker_Condition)
      biomarker.condition.patient <- unique(DF_patient$var.anno[DF_patient$VariantGene == gene_id])
      
      # Remove no match Biomarker_Condition from OnCore_Biomarker_Report
      for (cond_num in 1:length(biomarker.condition.trial)) {
        if (!(biomarker.condition.trial[cond_num] %in% biomarker.condition.patient)) {
          DF_Gene_OnCore_Biomarker <- DF_Gene_OnCore_Biomarker[DF_Gene_OnCore_Biomarker$Biomarker_Condition !=
                                                                 biomarker.condition.trial[cond_num],]
        }
      }
      
      # Update list 
      biomarker.condition.trial <- unique(DF_Gene_OnCore_Biomarker$Biomarker_Condition)
      
      ## Iterate through each bio.cond_id of biomarker.condition.patient
      #----------------------------------------------
      for (cond_num in 1:length(biomarker.condition.patient)) {
        bio.cond_id <- biomarker.condition.patient[cond_num]
        
        ## Match Biomarker_Condition
        #----------------------------------------------
        if (bio.cond_id %in% biomarker.condition.trial) {
          
          ## Parse Biomarker_Detail
          #----------------------------------------------
          biomarker.detail.trial <- 
            unique(DF_Gene_OnCore_Biomarker$Biomarker_Detail[DF_Gene_OnCore_Biomarker$Biomarker_Condition == bio.cond_id])
          biomarker.detail.patient <- 
            unique(DF_patient$VariantHGVSProtein[DF_patient$VariantGene == gene_id & DF_patient$var.anno == bio.cond_id])
          biomarker_detail_name_trial <- NA     # category match
          biomarker_detail_name_patient <- NA   # category match 
          
          ## Iterate through each element of biomarker.detail.trial
          #----------------------------------------------
          for (det_num in 1:length(biomarker.detail.trial)) {
            if (biomarker.detail.trial[det_num] == "All mutations accepted") {
              biomarker_detail_name_trial <- append(biomarker_detail_name_trial, biomarker.detail.trial[det_num])
              biomarker_detail_name_patient <- append(biomarker_detail_name_patient, biomarker.detail.patient)
              
            } else {
              ## Assume Amino acid change
              aa.start_trial <- gsub("(^[[:alpha:]]{3})(.*)","\\1", biomarker.detail.trial[det_num])
              aa.start_patient <- unique(DF_patient$aa.start[DF_patient$VariantGene == gene_id & DF_patient$var.anno == bio.cond_id])
              if (aa.start_trial == "*") {
                aa.end_match <- as.logical("TRUE")
                aa.start_list <- aa.start_patient
              } else {
                aa.start_match <- grepl(aa.start_trial, aa.start_patient)
                aa.start_list <- str_extract(aa.start_patient, aa.start_trial)
              }
              
              var.position_trial <- gsub("(^[[:alpha:]]{3})([[:digit:]]{,4})(.*)","\\2", biomarker.detail.trial[det_num])
              var.position_patient <- unique(DF_patient$var.position[DF_patient$VariantGene == gene_id & DF_patient$var.anno == bio.cond_id])
              if (var.position_trial == "*") {
                var.position_match <- as.logical("TRUE")
                var.position_list <- var.position_patient
              } else {
                var.position_match <- grepl(var.position_trial, var.position_patient)
                var.position_list <- str_extract(var.position_patient, var.position_trial)
              }
              
              aa.end_trial <- gsub("(^[[:alpha:]]{3}[[:digit:]]{,4})(.*)","\\2", biomarker.detail.trial[det_num])
              aa.end_patient <- unique(DF_patient$aa.end[DF_patient$VariantGene == gene_id & DF_patient$var.anno == bio.cond_id])
              if (aa.end_trial == "*") {
                aa.end_match <- as.logical("TRUE")
                aa.end_list <- aa.end_patient
              } else { 
                aa.end_match <- grepl(aa.end_trial, aa.end_patient)
                aa.end_list <- aa.end_patient
              }
              
              if ("TRUE" %in% aa.start_match) { aa.start_match <- TRUE } else { aa.start_match <- FALSE }
              if ("TRUE" %in% var.position_match) { var.position_match <- TRUE } else { var.position_match <- FALSE }
              if ("TRUE" %in% aa.end_match) { aa.end_match <- TRUE } else { aa.end_match <- FALSE }
              
              if (isTRUE(aa.start_match & var.position_match & aa.end_match)) {
                biomarker_detail_name_trial <- append(biomarker_detail_name_trial, biomarker.detail.trial[det_num])
                
                for (start_num in 1:length(aa.start_list)) {
                  for (pos_num in 1:length(var.position_list)) {
                    for (end_num in 1:length(aa.end_list)) {
                      protein_name <- paste("p.",aa.start_list[start_num],var.position_list[pos_num],
                                            aa.end_list[end_num], sep="")
                      biomarker_detail_name_patient <- append(biomarker_detail_name_patient, protein_name)      
                    }
                  }
                }
              }
            }
          }
          
          biomarker_detail_name_trial <- biomarker_detail_name_trial[which(!is.na(biomarker_detail_name_trial))]
          biomarker_detail_name_patient <- unique(biomarker_detail_name_patient[which(!is.na(biomarker_detail_name_patient))])
          
          # Remove no match Biomarker_Detail from OnCore_Biomarker_Report
          for (det_num in 1:length(biomarker.detail.trial)) {
            if (!(biomarker.detail.trial[det_num] %in% biomarker_detail_name_trial)) {
              DF_Gene_OnCore_Biomarker <- DF_Gene_OnCore_Biomarker[DF_Gene_OnCore_Biomarker$Biomarker_Detail !=
                                                                     biomarker.detail.trial[det_num],]
            }
          }
          
          ## Iterate through each bio.detail_id of biomarker.detail.patient
          #----------------------------------------------
          for (det_num in 1:length(biomarker.detail.patient)) {
            bio.detail_id <- biomarker.detail.patient[det_num]
            
            ## Match Biomarker_Detail
            #----------------------------------------------
            if (bio.detail_id %in% biomarker_detail_name_patient) {
              
              ## Assess Pathogenicity Status
              #----------------------------------------------
              pathogenicity_gate <- NA
              
              if (isTRUE(pathogenic_FILTER)) {
                if (DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                DF_patient$var.anno == bio.cond_id &
                                                                DF_patient$VariantHGVSProtein == bio.detail_id)] 
                    %in% pathogenic_accepted) {
                  pathogenicity_gate <- as.logical("TRUE")
                  
                } else {
                  pathogenicity_gate <- as.logical("FALSE")
                  DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id &
                                                          DF_patient$var.anno == bio.cond_id &
                                                          DF_patient$VariantHGVSProtein == bio.detail_id)] <-
                    "Pathogenicity criteria NOT satisfied"
                }
              }
              
              pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                             DF_patient$var.anno == bio.cond_id &
                                                                             DF_patient$VariantHGVSProtein == bio.detail_id)] 
              
              ## Match Pathogenicity Status
              #----------------------------------------------
              if (isTRUE(pathogenicity_gate == TRUE | is.na(pathogenicity_gate))) {
                
                ## Assess Age.Group criteria
                #----------------------------------------------
                age_gate <- NA
                
                if (isTRUE(adult.group_FILTER)) {
                  if (unique(DF_patient$PatientAge) >= 18) {
                    age_gate <- as.logical("TRUE")
                    
                  } else {
                    age_gate <- as.logical("FALSE")
                    DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id &
                                                            DF_patient$var.anno == bio.cond_id &
                                                            DF_patient$VariantHGVSProtein == bio.detail_id &
                                                            DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                      "Age criteria NOT satisfied"
                  }
                }
                
                ## Match Age.Group criteria
                #----------------------------------------------
                if (isTRUE(age_gate == TRUE | is.na(age_gate))) {
                  
                  ## Parse Disease.Group.category
                  #----------------------------------------------
                  # Each element of Disease.Group.category.trial is a single word
                  Disease.category.trial <- unique(DF_Gene_OnCore_Biomarker$Disease.Group.category)
                  # PrimaryTumorSite may match to several categories
                  Disease.category.patient <- unique(DF_patient$PrimaryTumorSite.Category)
                  disease_category_name_trial <- NA     # category match 
                  disease_category_name_patient <- NA   # category match 

                  ## Apply disease.group_FILTER
                  if (isTRUE(disease.group_FILTER)) {
                    
                    ## Iterate through each element of Disease.category.trial
                    #----------------------------------------------
                    for (cat_num in 1:length(Disease.category.trial)) {
                      if (isTRUE(Disease.category.trial[cat_num] == "any site")) {
                        disease_category_name_trial <- append(disease_category_name_trial, Disease.category.trial[cat_num])
                        disease_category_name_patient <- append(disease_category_name_patient, Disease.category.patient)

                      } else if (isTRUE(Disease.category.trial[cat_num] %in% Disease.category.patient)) {
                        disease_category_name_trial <- append(disease_category_name_trial, Disease.category.trial[cat_num])
                        disease_category_name_patient <- append(disease_category_name_patient, 
                                                                str_extract(Disease.category.patient,Disease.category.trial[cat_num]))

                      }
                    }
                    
                    disease_category_name_trial <- disease_category_name_trial[which(!is.na(disease_category_name_trial))]
                    disease_category_name_patient <- unique(disease_category_name_patient[which(!is.na(disease_category_name_patient))])
                    
                  } else {
                    ## Do not apply disease.group_FILTER
                    disease_category_name_trial <- Disease.category.trial
                    disease_category_name_patient <- Disease.category.patient
                  }
                  
                  # Remove no match Disease.Group.category from OnCore_Biomarker_Report
                  for (cat_num in 1:length(Disease.category.trial)) {
                    if (isTRUE(!(Disease.category.trial[cat_num] %in% disease_category_name_trial))) {
                      DF_Gene_OnCore_Biomarker <- DF_Gene_OnCore_Biomarker[DF_Gene_OnCore_Biomarker$Disease.Group.category !=
                                                                             Disease.category.trial[cat_num],]
                    }
                  }
                  
                  ## Iterate through each disease.cat_id of Disease.category.patient
                  #----------------------------------------------
                  for (cat_num in 1:length(Disease.category.patient)) {
                    disease.cat_id <- Disease.category.patient[cat_num]
                    
                    ## Match Disease.Group.category
                    #----------------------------------------------
                    if (isTRUE(disease.cat_id %in% disease_category_name_patient)) {
                      
                      ## Parse Disease.Site
                      #----------------------------------------------
                      # Each element of Disease.Site is a single word
                      Disease.Site.trial <- unique(DF_Gene_OnCore_Biomarker$Disease.Site)
                      # Disease.Site.patient is a single element
                      Disease.Site.patient <- unique(DF_patient$PrimaryTumorSite)
                      disease_site_name_trial <- NA     # category match
                      disease_site_name_patient <- NA   # category match 
                      
                      ## Apply disease.site_FILTER
                      if (isTRUE(disease.group_FILTER & disease.site_FILTER)) {
                        
                        ## Iterate through each element of Disease.Site.trial
                        #----------------------------------------------
                        for (site_num in 1:length(Disease.Site.trial)) {
                          if (isTRUE(Disease.Site.trial[site_num] == "any site")) {
                            disease_site_name_trial <- append(disease_site_name_trial, Disease.Site.trial[site_num])
                            disease_site_name_patient <- append(disease_site_name_patient, Disease.Site.patient)
                            
                          } else if (isTRUE(Disease.Site.trial[site_num] %in% Disease.Site.patient)) {
                            disease_site_name_trial <- append(disease_site_name_trial, Disease.Site.trial[site_num])
                            disease_site_name_patient <- append(disease_site_name_patient, 
                                                                str_extract(Disease.Site.patient, Disease.Site.trial[site_num]))
                          }
                        }
                        
                        disease_site_name_trial <- disease_site_name_trial[which(!is.na(disease_site_name_trial))]
                        disease_site_name_patient <- unique(disease_site_name_patient[which(!is.na(disease_site_name_patient))])
                        
                      } else {
                        ## Do not apply disease.group_FILTER
                        disease_site_name_trial <- Disease.Site.trial
                        disease_site_name_patient <- Disease.Site.patient
                      }
                      
                      # Remove no match Disease.Site from OnCore_Biomarker_Report
                      for (site_num in 1:length(Disease.Site.trial)) {
                        if (isTRUE(!(Disease.Site.trial[site_num] %in% disease_site_name_trial))) {
                          DF_Gene_OnCore_Biomarker <- DF_Gene_OnCore_Biomarker[DF_Gene_OnCore_Biomarker$Disease.Site !=
                                                                                 Disease.Site.trial[site_num],]
                        }
                      }
                      
                      ## Iterate through each disease.site_id of Disease.Site.patient
                      #----------------------------------------------
                      for (site_num in 1:length(Disease.Site.patient)) {
                        disease.site_id <- Disease.Site.patient[site_num]
                        
                        # Corresponding row in patient file
                        pt_rowNo <- which(DF_patient$VariantGene == gene_id &
                                            DF_patient$var.anno == bio.cond_id &
                                            DF_patient$VariantHGVSProtein == bio.detail_id &
                                            DF_patient$VariantPathogenicityStatus == pathogenic_id &
                                            DF_patient$PrimaryTumorSite.Category == disease.cat_id &
                                            DF_patient$PrimaryTumorSite == disease.site_id)
                        
                        ## Match Disease.Site
                        #----------------------------------------------
                        if (isTRUE(disease.site_id %in% disease_site_name_patient)) {
                          DF_patient$OnCore_Report_Status[pt_rowNo] <- "Candidate trial IDENTIFIED"
                      
                          ## Generate output file for candidate MATCH
                          #----------------------------------------------
                          # Corresponding OnCore.No
                          OnCore_No <- unique(DF_Gene_OnCore_Biomarker$OnCore.No)
                          
                          # Corresponding trial INFO in OnCore_Biomarker_Report
                          Trial_INFO <- OnCore_Biomarker_Report[OnCore_Biomarker_Report$OnCore.No %in% OnCore_No,]
                          
                          # Merge patient INFO with trial INFO
                          DF_Output_pre <- data.frame(matrix(NA, ncol = ncol(DF_Output_OnCore_Biomarker)))
                          colnames(DF_Output_pre) <- colnames(DF_Output_OnCore_Biomarker)
                          
                          for (core_num in 1:length(OnCore_No)) {
                            patient <- data.frame(DF_patient[pt_rowNo, c(1:ncol_STAMP)])
                            trial <- data.frame(Trial_INFO[core_num,])
                            DF_Output_pre[core_num,] <- cbind(patient,trial)
                          }
                          
                          # Append to Output file
                          DF_Output_OnCore_Biomarker <- rbind(DF_Output_OnCore_Biomarker, DF_Output_pre)

                        } else {
                          DF_patient$OnCore_Report_Status[pt_rowNo] <- "Disease Site criteria NOT satisfied"
                        }
                      }

                    } else {
                      DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id &
                                                              DF_patient$var.anno == bio.cond_id &
                                                              DF_patient$VariantHGVSProtein == bio.detail_id &
                                                              DF_patient$VariantPathogenicityStatus == pathogenic_id &
                                                              DF_patient$PrimaryTumorSite.Category == disease.cat_id)] <-
                        "Disease Group Category criteria NOT satisfied"
                    }
                  }
                }
              }
              
            } else {
              DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id &
                                                      DF_patient$var.anno == bio.cond_id &
                                                      DF_patient$VariantHGVSProtein == bio.detail_id)] <-
                "Biomarker Detail criteria NOT satisfied"
            }
          }
          
        } else {
          DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id &
                                                  DF_patient$var.anno == bio.cond_id)] <-
            "Biomarker Condition criteria NOT satisfied"
        }
      }
      
    } else {
      DF_patient$OnCore_Report_Status[which(DF_patient$VariantGene == gene_id)] <- 
        "Gene NOT found"
    }
  }
  
  # Write to match results per patient to local computer
  #----------------------------------------------
  write.table(DF_patient, file = paste(outdir_patient, patient_id, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

## Remove rows that are all empty
DF_Output_OnCore_Biomarker <- 
  DF_Output_OnCore_Biomarker[rowSums(is.na(DF_Output_OnCore_Biomarker)) != ncol(DF_Output_OnCore_Biomarker),]

## Write to match results for positive candidacy local computer
#----------------------------------------------
write.table(DF_Output_OnCore_Biomarker, 
            file = paste(outdir_int,"OnCore_Biomarker_Matched_", OnCore_Biomarker_Report_timestamp, "_", filterName, ".tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(DF_Gene_OnCore_Biomarker,DF_Output_pre,DF_patient,patient,trial,Trial_INFO,
       aa.end_list,aa.end_match,aa.end_patient,aa.end_trial,aa.start_list,aa.start_match,
       aa.start_patient,aa.start_trial,age_gate,bio.cond_id,bio.detail_id,
       biomarker_detail_name_patient,biomarker_detail_name_trial,
       biomarker.condition.patient,biomarker.condition.trial,biomarker.detail.patient,
       biomarker.detail.trial,cat_num,cond_num,core_num,det_num,
       disease_category_name_patient,disease_category_name_trial,
       disease_site_name_patient,disease_site_name_trial,disease.cat_id,
       Disease.category.patient,Disease.category.trial,disease.site_id,
       Disease.Site.patient,Disease.Site.trial,end_num,gene_id,gene_num,gene.patient,
       genes.OnCore_Biomarker,ncol_OnCore,ncol_STAMP,OnCore_No,pathogenic_id,
       pathogenicity_gate,patient_id,patient_num,pos_num,protein_name,pt_rowNo,
       site_num,start_num,var.position_list,var.position_match,var.position_patient,
       var.position_trial)

cat(paste("Timestamp of Internal clinical trial matching FINISH: ", Sys.time(), sep=""),"\n","\n")
