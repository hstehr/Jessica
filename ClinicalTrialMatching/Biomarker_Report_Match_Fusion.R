# Match by gene > biomarker.condition \
# > pathogenicity.status > age.group (i.e. adult) > disease.group.category > disease.site
## Search Output: individual patient file with stop point
## Match Output: "OnCore_Biomarker_Matched.tsv"

# Code Fusion dataframe: lines 11-14

if (isTRUE(Internal_match)) {
  
  ncol_OnCore <- as.numeric(ncol(OnCore_Biomarker_Report))
  ncol_Fusion = 10
  Fusion.colnames <- c("PatientID","Fusion_Detail","Gene","Break","var.type","var.anno",
                       "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus",
                       "HistologicalDx")
  
  # Extract gene from OnCore_Biomarker_QC
  genes.OnCore_Biomarker <- sort(unique(OnCore_Biomarker_QC$Biomarker_GeneName))
  
  # Output file of matched mutations with clinical trials
  DF_Output_Fusion_OnCore <- data.frame(matrix(NA, ncol = (ncol_Fusion + ncol_OnCore)))
  colnames(DF_Output_Fusion_OnCore) <- append(Fusion.colnames, colnames(OnCore_Biomarker_Report))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, "_Fusion.tsv", sep=""),
                           header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient) > 0) {
      # Extract unique genes per patient
      gene.patient <- sort(unique(DF_patient$Gene))
      
      # Clinical match result
      DF_patient$OnCore_Fusion_Status <- NA
      
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
          if (isTRUE("fusion" %in% unique(tolower(DF_Gene_OnCore_Biomarker$Biomarker_Condition)))) {
            
            # Specify Biomarker_Condition
            DF_Gene_OnCore_Biomarker <- 
              DF_Gene_OnCore_Biomarker[which(DF_Gene_OnCore_Biomarker$Biomarker_Condition == "FUSION"), ]
            
            ## Assess Pathogenicity Status
            #----------------------------------------------
            pathogenicity_gate <- NA
            
            if (isTRUE(pathogenic_FILTER)) {
              if (DF_patient$VariantPathogenicityStatus[which(DF_patient$Gene == gene_id &
                                                              DF_patient$var.anno == "Fusion")] 
                  %in% pathogenic_accepted) {
                pathogenicity_gate <- as.logical("TRUE")
                
              } else {
                pathogenicity_gate <- as.logical("FALSE")
                DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id &
                                                        DF_patient$var.anno == "Fusion")] <-
                  "Pathogenicity criteria NOT satisfied"
              }
            }
            
            pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$Gene == gene_id &
                                                                           DF_patient$var.anno == "Fusion")]
            
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
                  DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id &
                                                          DF_patient$var.anno == "Fusion" &
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
                disease_category_name_trial <- NA
                disease_category_name_patient <- NA
                
                ## Apply disease.group_FILTER
                # Account for missing primary tumor site information
                ## PrimaryTumorSite.Category != "unknown" is equivalent to disease.group_FILTER = "FALSE"
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
                    ## PrimaryTumorSite.Category != "unknown" is equivalent to disease.site_FILTER = "FALSE"
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
                      pt_rowNo <- which(DF_patient$Gene == gene_id &
                                          DF_patient$var.anno == "Fusion" &
                                          DF_patient$VariantPathogenicityStatus %in% pathogenic_id)
                      
                      ## Match Disease.Site
                      #----------------------------------------------
                      if (isTRUE(disease.site_id %in% disease_site_name_patient)) {
                        
                        ## Generate output file for candidate MATCH
                        #----------------------------------------------
                        # Corresponding OnCore.No
                        OnCore_No <- unique(DF_Gene_OnCore_Biomarker$OnCore.No)
                        
                        DF_patient$OnCore_Fusion_Status[pt_rowNo] <- 
                          paste("Candidate trial IDENTIFIED: ", paste(OnCore_No, collapse = ", "), sep="")
                        
                        # Corresponding trial INFO in OnCore_Biomarker_Report
                        Trial_INFO <- OnCore_Biomarker_Report[which(OnCore_Biomarker_Report$OnCore.No %in% OnCore_No),]
                        
                        # Merge patient INFO with trial INFO
                        DF_Output_pre <- data.frame(matrix(NA, ncol = ncol(DF_Output_Fusion_OnCore)))
                        colnames(DF_Output_pre) <- colnames(DF_Output_Fusion_OnCore)
                        
                        patient <- data.frame(DF_patient[pt_rowNo, c(1:ncol_Fusion)])
                        DF_Output_pre <- crossing(patient, Trial_INFO)
                        
                        # Append to Output file
                        DF_Output_Fusion_OnCore <- rbind(DF_Output_Fusion_OnCore, DF_Output_pre)
                        
                        remove(DF_Output_pre,Trial_INFO,OnCore_No)
                        
                      } else {
                        DF_patient$OnCore_Fusion_Status[pt_rowNo] <- "Disease Site criteria NOT satisfied"
                      }
                      remove(disease.site_id,pt_rowNo)
                    }
                    remove(Disease.Site.patient,disease_site_name_patient,disease_site_name_trial,
                           Disease.Site.trial,site_num)
                    
                  } else {
                    if (isTRUE(stop_checkpoint_na)) {
                      DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id &
                                                              DF_patient$var.anno == bio.cond_id &
                                                              is.na(DF_patient$VariantHGVSProtein) &
                                                              DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                        "Disease Group Category criteria NOT satisfied"
                      
                    } else {
                      DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id &
                                                              DF_patient$var.anno == bio.cond_id &
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
            DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id &
                                                    DF_patient$var.anno == "Fusion")] <-
              "Biomarker Condition criteria NOT satisfied"
          }
          remove(DF_Gene_OnCore_Biomarker)
          
        } else {
          DF_patient$OnCore_Fusion_Status[which(DF_patient$Gene == gene_id)] <- 
            "Gene NOT found"
        }
        remove(gene_id,gene_num)
      }
      remove(gene.patient)
      
      # Write to match results per patient to local computer
      #----------------------------------------------
      write.table(DF_patient, file = paste(tempdir, patient_id, "_Fusion.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    remove(patient_id,patient_num,DF_patient)
  }
  
  ## Remove rows that are all empty
  DF_Output_Fusion_OnCore <- 
    DF_Output_Fusion_OnCore[rowSums(is.na(DF_Output_Fusion_OnCore)) != ncol(DF_Output_Fusion_OnCore),]
  
  ## Write to match results for positive candidacy local computer
  #----------------------------------------------
  write.table(DF_Output_Fusion_OnCore, 
              file = paste(tempdir,"OnCore_Fusion_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                           groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  remove(genes.OnCore_Biomarker,ncol_OnCore,ncol_Fusion,Fusion.colnames,DF_Output_Fusion_OnCore)
}
