# Match by gene > biomarker.condition > biomarker.detail  (if specified, match only for SNVs) \
# > pathogenicity.status > age.group (i.e. adult) > disease.group.category > disease.site
## Search Output: individual patient file with stop point
## Match Output: "OnCore_Matched_SNVIndel.tsv"

# Code SNVIndel dataframe: lines 11-15

if (isTRUE(Internal_match)) {
  cat(paste("Timestamp of Stanford OnCore clinical trial matching START: ", Sys.time(), sep=""),"\n","\n")
  
  ncol_OnCore <- as.numeric(ncol(OnCore_Biomarker_Report))
  ncol_SNVIndel = 13
  SNVIndel.colnames <- c("PatientID","VariantHGVSGenomic","VariantLabel","VariantGene","VariantHGVSCoding",
                         "VariantHGVSProtein","var.type","var.anno","Exon_Number",
                         "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus",
                         "HistologicalDx")
  
  # Extract gene from OnCore_Biomarker_QC
  genes.OnCore_Biomarker <- sort(unique(OnCore_Biomarker_QC$Biomarker_GeneName))
  
  # Output file of matched mutations with clinical trials
  DF_Output_SNVIndel_OnCore <- data.frame(matrix(NA, ncol = (ncol_SNVIndel + ncol_OnCore)))
  colnames(DF_Output_SNVIndel_OnCore) <- append(SNVIndel.colnames, colnames(OnCore_Biomarker_Report))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                           header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient) > 0) {
      # Extract unique genes per patient
      gene.patient <- sort(unique(DF_patient$VariantGene))
      
      # Clinical match result
      DF_patient$OnCore_SNVIndel_Status <- NA
      
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
          if (isTRUE("MUTATION" %in% unique(DF_Gene_OnCore_Biomarker$Biomarker_Condition))) {
            
            # Specify Biomarker_Condition == "MUTATION"
            DF_Gene_OnCore_Biomarker <- 
              DF_Gene_OnCore_Biomarker[which(DF_Gene_OnCore_Biomarker$Biomarker_Condition == "MUTATION"), ]
            
            ## Parse Biomarker_Detail == possible only for SNVs
            #----------------------------------------------
            biomarker.detail.trial <- 
              unique(DF_Gene_OnCore_Biomarker$Biomarker_Detail[DF_Gene_OnCore_Biomarker$Biomarker_Condition == "MUTATION"])
            biomarker.detail.patient <- 
              unique(DF_patient$VariantHGVSProtein[DF_patient$VariantGene == gene_id & DF_patient$var.anno == "MUTATION"])
            biomarker_detail_name_trial <- NA
            biomarker_detail_name_patient <- NA
            
            ## Iterate through each element of biomarker.detail.trial
            #----------------------------------------------
            for (det_num in 1:length(biomarker.detail.trial)) {
              if (biomarker.detail.trial[det_num] == "All mutations accepted") {
                biomarker_detail_name_trial <- append(biomarker_detail_name_trial, biomarker.detail.trial[det_num])
                biomarker_detail_name_patient <- append(biomarker_detail_name_patient, biomarker.detail.patient)
                
              } else {
                ## Assume Amino acid change
                aa.start_trial <- gsub("(^[[:alpha:]])(.*)","\\1", biomarker.detail.trial[det_num])
                aa.start_patient <- gsub("(^p.)([[:alpha:]]{1})(.*)","\\2", 
                                         unique(DF_patient$VariantHGVSProtein[
                                           which(DF_patient$VariantGene == gene_id & DF_patient$var.anno == "MUTATION")]))
                
                if (aa.start_trial == "*") {
                  aa.end_match <- as.logical("TRUE")
                  aa.start_list <- aa.start_patient
                } else {
                  aa.start_match <- grepl(aa.start_trial, aa.start_patient)
                  aa.start_list <- str_extract(aa.start_patient, aa.start_trial)
                }
                
                var.position_trial <- gsub("(^[[:alpha:]]{1})([[:digit:]]{,4})(.*)","\\2", biomarker.detail.trial[det_num])
                var.position_patient <- gsub("(^p.[[:alpha:]]{1})([[:digit:]]{,4})(.*)","\\2", 
                                             unique(DF_patient$VariantHGVSProtein[
                                               which(DF_patient$VariantGene == gene_id & DF_patient$var.anno == "MUTATION")]))
                
                if (var.position_trial == "*") {
                  var.position_match <- as.logical("TRUE")
                  var.position_list <- var.position_patient
                } else {
                  var.position_match <- grepl(var.position_trial, var.position_patient)
                  var.position_list <- str_extract(var.position_patient, var.position_trial)
                }
                
                aa.end_trial <- gsub("(^[[:alpha:]]{3}[[:digit:]]{,4})(.*)","\\2", biomarker.detail.trial[det_num])
                aa.end_patient <- gsub("(^p.[[:alpha:]]{1}[[:digit:]]{,4})(.*)","\\2", 
                                       unique(DF_patient$VariantHGVSProtein[
                                         which(DF_patient$VariantGene == gene_id & DF_patient$var.anno == "MUTATION")]))
                
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
                  remove(end_num,pos_num,protein_name,start_num)
                }
                remove(aa.end_list,aa.end_match,aa.end_patient,aa.end_trial,
                       aa.start_list,aa.start_match,aa.start_patient,aa.start_trial,
                       var.position_list,var.position_match,var.position_patient,var.position_trial)
              }
            }
            
            biomarker_detail_name_trial <- biomarker_detail_name_trial[which(!is.na(biomarker_detail_name_trial))]
            biomarker_detail_name_patient <- unique(biomarker_detail_name_patient[which(!is.na(biomarker_detail_name_patient))])
            
            # Remove no match Biomarker_Detail from OnCore_Biomarker_Report
            for (det_num in 1:length(biomarker.detail.trial)) {
              if (isTRUE(!(biomarker.detail.trial[det_num] %in% biomarker_detail_name_trial))) {
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
              # Account for missing HGVS protein information
              if (bio.detail_id %in% biomarker_detail_name_patient) {
                stop_checkpoint <- as.logical("TRUE")
                stop_checkpoint_na <- as.logical("FALSE")
              } else if (isTRUE(is.na(bio.detail_id) & ("All mutations accepted" %in% biomarker.detail.trial))) {
                stop_checkpoint <- as.logical("TRUE")
                stop_checkpoint_na <- as.logical("TRUE")
              }
              
              if (isTRUE(stop_checkpoint)) {
                
                ## Assess Pathogenicity Status
                #----------------------------------------------
                pathogenicity_gate <- NA
                
                if (isTRUE(pathogenic_FILTER)) {
                  if (isTRUE(DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                         DF_patient$var.anno == "MUTATION" &
                                                                         DF_patient$VariantHGVSProtein == bio.detail_id)] 
                             %in% pathogenic_accepted)) {
                    pathogenicity_gate <- as.logical("TRUE")
                    
                  } else if (isTRUE(DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                                DF_patient$var.anno == "MUTATION" &
                                                                                is.na(DF_patient$VariantHGVSProtein))] 
                                    %in% pathogenic_accepted) & stop_checkpoint_na) {
                    pathogenicity_gate <- as.logical("TRUE")
                    
                  } else {
                    pathogenicity_gate <- as.logical("FALSE")
                    
                    if (isTRUE(stop_checkpoint_na)) {
                      DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                DF_patient$var.anno == "MUTATION" &
                                                                is.na(DF_patient$VariantHGVSProtein))] <-
                        "Pathogenicity criteria NOT satisfied"
                    } else {
                      DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                DF_patient$var.anno == "MUTATION" &
                                                                DF_patient$VariantHGVSProtein == bio.detail_id)] <-
                        "Pathogenicity criteria NOT satisfied"
                    }
                  }
                }
                
                if (isTRUE(stop_checkpoint_na)) {
                  pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                                 DF_patient$var.anno == "MUTATION" &
                                                                                 is.na(DF_patient$VariantHGVSProtein))] 
                } else {
                  pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                                 DF_patient$var.anno == "MUTATION" &
                                                                                 DF_patient$VariantHGVSProtein == bio.detail_id)]   
                }
                
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
                      
                      if (isTRUE(stop_checkpoint_na)) {
                        DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                  DF_patient$var.anno == "MUTATION" &
                                                                  is.na(DF_patient$VariantHGVSProtein) &
                                                                  DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                          "Age criteria NOT satisfied"
                        
                      } else {
                        DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                  DF_patient$var.anno == "MUTATION" &
                                                                  DF_patient$VariantHGVSProtein == bio.detail_id &
                                                                  DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                          "Age criteria NOT satisfied"
                      }
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
                    disease_category_name_trial <- NA
                    disease_category_name_patient <- NA
                    
                    ## Apply disease.group_FILTER
                    # Account for missing primary tumor site information
                    ## PrimaryTumorSite.Category == "unknown" is equivalent to disease.group_FILTER = "FALSE"
                    if (isTRUE(disease.group_FILTER & Disease.category.patient != "unknown")) {
                      
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
                        disease_site_name_trial <- NA
                        disease_site_name_patient <- NA
                        
                        ## Apply disease.site_FILTER
                        # Account for missing primary tumor site information
                        ## PrimaryTumorSite.Category == "unknown" is equivalent to disease.site_FILTER = "FALSE"
                        if (isTRUE(disease.group_FILTER & disease.site_FILTER & Disease.category.patient != "unknown")) {
                          
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
                          ## Do not apply disease.site_FILTER
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
                          # PrimaryTumorSite.Category and PrimaryTumorSite is consistent throughout DF_patient
                          if (isTRUE(stop_checkpoint_na)) {
                            pt_rowNo <- which(DF_patient$VariantGene == gene_id &
                                                DF_patient$var.anno == "MUTATION" &
                                                is.na(DF_patient$VariantHGVSProtein) &
                                                DF_patient$VariantPathogenicityStatus %in% pathogenic_id)
                          } else {
                            pt_rowNo <- which(DF_patient$VariantGene == gene_id &
                                                DF_patient$var.anno == "MUTATION" &
                                                DF_patient$VariantHGVSProtein == bio.detail_id &
                                                DF_patient$VariantPathogenicityStatus %in% pathogenic_id)
                          }
                          
                          ## Match Disease.Site
                          #----------------------------------------------
                          if (isTRUE(disease.site_id %in% disease_site_name_patient)) {
                            
                            ## Generate output file for candidate MATCH
                            #----------------------------------------------
                            # Corresponding OnCore.No
                            OnCore_No <- unique(DF_Gene_OnCore_Biomarker$OnCore.No)
                            
                            DF_patient$OnCore_SNVIndel_Status[pt_rowNo] <- 
                              paste("Candidate trial IDENTIFIED: ", paste(OnCore_No, collapse = ", "), sep="")
                            
                            # Corresponding trial INFO in OnCore_Biomarker_Report
                            Trial_INFO <- OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$OnCore.No %in% OnCore_No),]
                            
                            # Merge patient INFO with trial INFO
                            DF_Output_pre <- data.frame(matrix(NA, ncol = ncol(DF_Output_SNVIndel_OnCore)))
                            colnames(DF_Output_pre) <- colnames(DF_Output_SNVIndel_OnCore)
                            
                            patient <- data.frame(DF_patient[pt_rowNo, c(1:ncol_SNVIndel)])
                            DF_Output_pre <- crossing(patient, Trial_INFO)
                            
                            # Append to Output file
                            DF_Output_SNVIndel_OnCore <- rbind(DF_Output_SNVIndel_OnCore, DF_Output_pre)
                            
                            remove(DF_Output_pre,Trial_INFO,OnCore_No,patient)
                            
                          } else {
                            DF_patient$OnCore_SNVIndel_Status[pt_rowNo] <- "Disease Site criteria NOT satisfied"
                          }
                          remove(disease.site_id,pt_rowNo)
                        }
                        remove(Disease.Site.patient,disease_site_name_patient,disease_site_name_trial,
                               Disease.Site.trial,site_num)
                        
                      } else {
                        if (isTRUE(stop_checkpoint_na)) {
                          DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                    DF_patient$var.anno == "MUTATION" &
                                                                    is.na(DF_patient$VariantHGVSProtein) &
                                                                    DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                            "Disease Group Category criteria NOT satisfied"
                          
                        } else {
                          DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                                    DF_patient$var.anno == "MUTATION" &
                                                                    DF_patient$VariantHGVSProtein == bio.detail_id &
                                                                    DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                            "Disease Group Category criteria NOT satisfied"
                        }
                      }
                      remove(disease.cat_id)
                    }
                    remove(Disease.category.trial,Disease.category.patient,disease_category_name_trial,
                           disease_category_name_patient,cat_num)
                  }
                  remove(age_gate)
                }
                remove(pathogenicity_gate,pathogenic_id)
                
              } else {
                if (isTRUE(stop_checkpoint_na)) {
                  DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                            DF_patient$var.anno == "MUTATION" &
                                                            is.na(DF_patient$VariantHGVSProtein))] <-
                    "Biomarker Detail criteria NOT satisfied"
                  
                } else {
                  DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                            DF_patient$var.anno == "MUTATION" &
                                                            DF_patient$VariantHGVSProtein == bio.detail_id)] <-
                    "Biomarker Detail criteria NOT satisfied"
                }
              }
              remove(bio.detail_id,stop_checkpoint_na,stop_checkpoint)
            }
            remove(biomarker_detail_name_patient,biomarker_detail_name_trial,biomarker.detail.patient,
                   biomarker.detail.trial,det_num)
            
          } else {
            DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id &
                                                      DF_patient$var.anno == "MUTATION")] <-
              "Biomarker Condition criteria NOT satisfied"
          }
          remove(DF_Gene_OnCore_Biomarker)
          
        } else {
          DF_patient$OnCore_SNVIndel_Status[which(DF_patient$VariantGene == gene_id)] <- 
            "Gene NOT found"
        }
        remove(gene_id,gene_num)
      }
      remove(gene.patient)
      
      # Write to match results per patient to local computer
      #----------------------------------------------
      write.table(DF_patient, file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    remove(patient_id,patient_num,DF_patient)
  }
  
  ## Remove rows that are all empty
  DF_Output_SNVIndel_OnCore <- 
    DF_Output_SNVIndel_OnCore[rowSums(is.na(DF_Output_SNVIndel_OnCore)) != ncol(DF_Output_SNVIndel_OnCore),]
  
  ## Write to match results for positive candidacy local computer
  #----------------------------------------------
  write.table(DF_Output_SNVIndel_OnCore, 
              file = paste(tempdir,"OnCore_SNVIndel_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                           groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  remove(genes.OnCore_Biomarker,ncol_OnCore,ncol_SNVIndel,DF_Output_SNVIndel_OnCore,SNVIndel.colnames)
}
