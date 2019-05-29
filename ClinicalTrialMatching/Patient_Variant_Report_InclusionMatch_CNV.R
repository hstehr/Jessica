# Match by gene > variant type > variant annotation (ie. amplification, deletion) \
# > pathogenicity.status > age.group (i.e. adult) > exclusion variant > exclusion nonhotspot
## Search Output: individual patient file with stop point
## Match Output: "NCI_CNV_Variant_Matched.tsv"

# Code CNV dataframe: lines 11-14

if (isTRUE(NCI_match)) {
  
  ncol_InclusionVariants <- as.numeric(ncol(Inclusion_Variants))
  ncol_CNV = 12
  CNV.colnames <- c("PatientID","CNV_Gene","Locus","Tiles","mean.z","mcopies","var.anno","var.type",
                    "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus",
                    "HistologicalDx")
  
  # Extract gene from NCI-MATCH
  genes.Inclusion_Variants <- sort(unique(Inclusion_Variants$Gene_Name))
  
  # Output file of matched mutations with clinical trials
  DF_Output_Patient_CNV_Variant <- data.frame(matrix(NA, ncol = (ncol_CNV + ncol_InclusionVariants)))
  colnames(DF_Output_Patient_CNV_Variant) <- append(CNV.colnames, colnames(Inclusion_Variants))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    DF_Output_Patient_CNV_Variant_int <- data.frame(matrix(NA, ncol = (ncol_CNV + ncol_InclusionVariants)))
    colnames(DF_Output_Patient_CNV_Variant_int) <- append(CNV.colnames, colnames(Inclusion_Variants))
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, "_CNV.tsv", sep=""),
                           header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient) > 0) {
      # Structure annotation for consistency
      DF_patient$var.anno <- gsub("AMP","amplification",DF_patient$var.anno)
      DF_patient$var.anno <- gsub("DEL","deletion",DF_patient$var.anno)
      
      # Extract unique genes per patient
      gene.patient <- sort(unique(DF_patient$CNV_Gene))
      
      # Clinical match result
      DF_patient$NCI_CNV_Variant_Status <- NA
      
      ## Iterate through each gene_id of gene.patient
      #----------------------------------------------
      for (gene_num in 1:length(gene.patient)) {
        gene_id <- gene.patient[gene_num]
        
        ## Match Biomarker_GeneName
        #----------------------------------------------
        if (gene_id %in% genes.Inclusion_Variants) {
          
          # Extract gene_id match from Inclusion_Variants 
          DF_Gene_Patient_Variant <- 
            data.frame(Inclusion_Variants[which(Inclusion_Variants$Gene_Name == gene_id),])
          
          if (isTRUE("CNV" %in% unique(DF_Gene_Patient_Variant$Variant_Type))) {
            
            # Specify var_type == "CNV"
            DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[which(DF_Gene_Patient_Variant$Variant_Type == "CNV"),]
            
            # Modify protein structure for consistency
            DF_Gene_Patient_Variant$Protein <- 
              tolower(gsub("(^[[:upper:]]{2,}.*[[:blank:]]+)([[:alpha:]]+)","\\2",DF_Gene_Patient_Variant$Protein))
            
            ## Parse Variant annotation
            #----------------------------------------------
            var.anno.trial <- unique(DF_Gene_Patient_Variant$Protein)
            var.anno.patient <- unique(DF_patient$var.anno[DF_patient$CNV_Gene == gene_id])
            
            # Remove no match Protein from DF_Gene_Patient_Variant
            for (anno_num in 1:length(var.anno.trial)) {
              if (isTRUE(!(var.anno.trial[anno_num] %in% var.anno.patient))) {
                DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Protein !=
                                                                     var.anno.trial[anno_num],]
              }
            }
            
            # Update list 
            var.anno.trial <- unique(DF_Gene_Patient_Variant$Protein)
            
            ## Iterate through each var.anno_id of var.anno.patient
            #----------------------------------------------
            for (anno_num in 1:length(var.anno.patient)) {
              var.anno_id <- var.anno.patient[anno_num]
              
              ## Match Variant_Anno
              #----------------------------------------------
              if (isTRUE(var.anno_id %in% var.anno.trial & nrow(DF_Gene_Patient_Variant) > 0)) {
                
                ## Assess Pathogenicity Status
                #----------------------------------------------
                pathogenicity_gate <- NA
                
                if (isTRUE(pathogenic_FILTER)) {
                  if (DF_patient$VariantPathogenicityStatus[which(DF_patient$CNV_Gene == gene_id &
                                                                  DF_patient$var.type == "CNV" &
                                                                  DF_patient$var.anno == var.anno_id)] 
                      %in% pathogenic_accepted) {
                    pathogenicity_gate <- as.logical("TRUE")
                    
                  } else {
                    pathogenicity_gate <- as.logical("FALSE")
                    DF_patient$NCI_CNV_Variant_Status[which(DF_patient$CNV_Gene == gene_id &
                                                              DF_patient$var.type == "CNV" &
                                                              DF_patient$var.anno == var.anno_id)] <-
                      "Pathogenicity criteria NOT satisfied"
                  }
                }
                
                pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$CNV_Gene == gene_id &
                                                                               DF_patient$var.type == "CNV" &
                                                                               DF_patient$var.anno == var.anno_id)]
                
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
                      DF_patient$NCI_CNV_Variant_Status[which(DF_patient$CNV_Gene == gene_id &
                                                                DF_patient$var.type == "CNV" &
                                                                DF_patient$var.anno == var.anno_id &
                                                                DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                        "Age criteria NOT satisfied"
                    }
                  }
                  
                  ## Match Age.Group criteria
                  #----------------------------------------------
                  if (isTRUE(age_gate == TRUE | is.na(age_gate))) {
                    
                    # Corresponding row in patient file
                    pt_rowNo <- which(DF_patient$CNV_Gene == gene_id &
                                        DF_patient$var.type == "CNV" &
                                        DF_patient$var.anno == var.anno_id &
                                        DF_patient$VariantPathogenicityStatus == pathogenic_id)
                    
                    arm.match <- unique(DF_Gene_Patient_Variant$Arm_Name
                                        [which(DF_Gene_Patient_Variant$Gene_Name == gene_id &
                                                 DF_Gene_Patient_Variant$Variant_Type == "CNV" &
                                                 DF_Gene_Patient_Variant$Protein == var.anno_id)])
                    
                    DF_patient$NCI_CNV_Variant_Status[pt_rowNo] <- paste("Candidate trial IDENTIFIED: ", 
                                                                         paste(arm.match, collapse = ", "), sep="")
                    
                    ## Assess Exclusions
                    #----------------------------------------------
                    if (length(arm.match) > 0) {
                      
                      for (arm_num in 1:length(arm.match)) {
                        
                        arm_id <- arm.match[arm_num]
                        exclusion_continue <- as.logical("TRUE")
                        
                        # Import SNV/Indel entries per patient
                        #----------------------------------------------
                        DF_patient_SNVIndel <- read.csv(file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                                                        header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                        
                        if (nrow(DF_patient_SNVIndel) > 0) {exclusion_SNV_continue <- as.logical("TRUE")  
                        } else {exclusion_SNV_continue <- as.logical("FALSE")}
                        
                        # Import Fusion entries per patient
                        #----------------------------------------------
                        DF_patient_Fusion <- read.csv(file = paste(tempdir, patient_id, "_Fusion.tsv", sep=""),
                                                      header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                        
                        if (nrow(DF_patient_Fusion) > 0) {exclusion_Fusion_continue <- as.logical("TRUE")  
                        } else {exclusion_Fusion_continue <- as.logical("FALSE")}
                        
                        ## Assess Exclusion Variants
                        #----------------------------------------------
                        # Extract variants for Arm_No in Exclusion_Variants
                        DF_Exclude_Arm <- Exclusion_Variants[which(Exclusion_Variants$Arm_Name == arm_id), ]
                        
                        if (nrow(DF_Exclude_Arm) > 0) {
                          
                          ## Assess Exclusion Variants: CNVs
                          #----------------------------------------------
                          DF_Exclude_Arm_CNV <- DF_Exclude_Arm[which(DF_Exclude_Arm$Variant_Type == "CNV"),]
                          
                          DF_Exclude_Arm_CNV$Protein <- 
                            tolower(gsub("(^[[:upper:]]{2,}.*[[:blank:]]+)([[:alpha:]]+)","\\2",DF_Exclude_Arm_CNV$Protein))
                          
                          if (nrow(DF_Exclude_Arm_CNV) > 0) {
                            
                            # All variant labels for Arm_No in Exclusion_Variants
                            var.exclude <- tolower(sort(unique(paste(DF_Exclude_Arm_CNV$Variant_ID, " ", DF_Exclude_Arm_CNV$Protein, sep=""))))
                            
                            # All variant labels for patient in DF_patient
                            var.patient <- tolower(sort(unique(paste(DF_patient$CNV_Gene, " ", DF_patient$var.anno, sep=""))))
                            
                            # Identify overlapping variants
                            var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                            
                            ## Match Exclusion Variants
                            #----------------------------------------------
                            if (length(var.exclude) > 0) {
                              
                              # Update Match status with exclusion
                              exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                         paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                         " (CNV)", sep="")
                              DF_patient$NCI_CNV_Variant_Status[pt_rowNo] <- exclusion_comment
                              
                              cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                              
                              # Remove excluded ARM_No from DF_Gene_Patient_Variant
                              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                              
                              DF_Output_Patient_CNV_Variant_int <- DF_Output_Patient_CNV_Variant_int[DF_Output_Patient_CNV_Variant_int$Arm_Name != arm_id,]
                              
                              # Checkpoint
                              exclusion_continue <- as.logical("FALSE")
                            }
                          }
                          
                          ## Assess Exclusion Variants: SNVs
                          #----------------------------------------------
                          DF_Exclude_Arm_SNV <- DF_Exclude_Arm[which(DF_Exclude_Arm$Variant_Type == "SNV"),]
                          
                          if (isTRUE(nrow(DF_Exclude_Arm_SNV) > 0 & exclusion_continue & exclusion_SNV_continue)) {
                            
                            # All variant labels for Arm_No in Exclusion_Variants
                            var.exclude <- sort(unique(paste(DF_Exclude_Arm_SNV$Gene_Name, " ", DF_Exclude_Arm_SNV$Genomic, " (", 
                                                             DF_Exclude_Arm_SNV$Variant_Type, ")", sep="")))
                            
                            # All variant labels for patient in DF_patient
                            var.patient <- sort(unique(paste(DF_patient_SNVIndel$VariantGene, " ", DF_patient_SNVIndel$VariantHGVSGenomic, " (",
                                                             DF_patient_SNVIndel$var.type, ")", sep="")))
                            
                            # Identify overlapping variants
                            var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                            
                            ## Match Exclusion Variants
                            #----------------------------------------------
                            if (length(var.exclude) > 0) {
                              
                              # Update Match status with exclusion
                              exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                         paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                         " mutation", sep="")
                              DF_patient$NCI_CNV_Variant_Status[pt_rowNo] <- exclusion_comment
                              
                              cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                              
                              # Remove excluded ARM_No from DF_Gene_Patient_Variant
                              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                              
                              DF_Output_Patient_CNV_Variant_int <- DF_Output_Patient_CNV_Variant_int[DF_Output_Patient_CNV_Variant_int$Arm_Name != arm_id,]
                              
                              # Checkpoint
                              exclusion_continue <- as.logical("FALSE")
                            }
                          }
                          
                          ## Assess Exclusion Variants: Fusions
                          #----------------------------------------------
                          DF_Exclude_Arm_Fusion <- DF_Exclude_Arm[which(DF_Exclude_Arm$Variant_Type == "Fusion"),]
                          
                          if (isTRUE(nrow(DF_Exclude_Arm_Fusion) > 0 & exclusion_continue & exclusion_Fusion_continue)) {
                            
                            # All variant labels for Arm_No in Exclusion_Variants
                            var.exclude <- sort(unique(paste(DF_Exclude_Arm_Fusion$Gene_Name, " (", 
                                                             DF_Exclude_Arm_Fusion$Variant_Type, ")", sep="")))
                            
                            # All variant labels for patient in DF_patient
                            var.patient <- sort(unique(paste(DF_patient_Fusion$Gene, " (",
                                                             DF_patient_Fusion$var.type, ")", sep="")))
                            
                            # Identify overlapping variants
                            var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                            
                            ## Match Exclusion Variants
                            #----------------------------------------------
                            if (length(var.exclude) > 0) {
                              
                              # Update Match status with exclusion
                              exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                         paste(sub("/.*", "", var.exclude), collapse = " & "), sep="")
                              DF_patient$NCI_CNV_Variant_Status[pt_rowNo] <- exclusion_comment
                              
                              cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                              
                              # Remove excluded ARM_No from DF_Gene_Patient_Variant
                              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                              
                              DF_Output_Patient_CNV_Variant_int <- DF_Output_Patient_CNV_Variant_int[DF_Output_Patient_CNV_Variant_int$Arm_Name != arm_id,]
                              
                              # Checkpoint
                              exclusion_continue <- as.logical("FALSE")
                            }
                          }
                          remove(DF_Exclude_Arm_SNV, DF_Exclude_Arm_CNV, DF_Exclude_Arm_Fusion)
                        }
                        remove(DF_Exclude_Arm)
                        
                        ## Assess Disease Exclusion Codes
                        #----------------------------------------------
                        # dx.patient == "NULL" indicates that histological dx is not available 
                        dx.patient = unique(DF_patient$HistologicalDx)
                        tumorsite.patient = unique(DF_patient$PrimaryTumorSite)
                        
                        ## HistologicalDx == "NULL" is equivalent to disease.site_FILTER = "FALSE"
                        if (isTRUE(!is.null(dx.patient) & exclusion_continue & disease.code_FILTER)) {
                          
                          # Extract variants for Arm_No in Disease Exclusion Codes
                          DF_Disease_Exclude_Arm <- Disease_Exclusion_Codes[Disease_Exclusion_Codes$Arm_Name == arm_id,]
                          DF_Disease_Exclude_patient <-
                            HistologicalDxCategory[which(HistologicalDxCategory$histologicalDiagnosis == dx.patient &
                                                           HistologicalDxCategory$primaryTumorSite == tumorsite.patient),]
                          
                          # Identify overlapping classifications 
                          CTEP.CATEGORY <- 
                            DF_Disease_Exclude_patient$CTEP.CATEGORY[which(DF_Disease_Exclude_patient$CTEP.CATEGORY %in% tolower(DF_Disease_Exclude_Arm$CTEP.CATEGORY))]
                          CTEP.CATEGORY <- CTEP.CATEGORY[which(!is.na(CTEP.CATEGORY))]
                          
                          CTEP.SUBCATEGORY <- 
                            DF_Disease_Exclude_patient$CTEP.SUBCATEGORY[which(DF_Disease_Exclude_patient$CTEP.SUBCATEGORY %in% tolower(DF_Disease_Exclude_Arm$CTEP.SUBCATEGORY))]
                          CTEP.SUBCATEGORY <- CTEP.SUBCATEGORY[which(!is.na(CTEP.SUBCATEGORY))]
                          
                          CTEP.TERM <- 
                            DF_Disease_Exclude_patient$CTEP.TERM[which(DF_Disease_Exclude_patient$CTEP.TERM %in% tolower(DF_Disease_Exclude_Arm$CTEP.TERM))]
                          CTEP.TERM <- CTEP.TERM[which(!is.na(CTEP.TERM))]
                          
                          SHORT.NAME <- 
                            DF_Disease_Exclude_patient$SHORT.NAME[which(DF_Disease_Exclude_patient$SHORT.NAME %in% tolower(DF_Disease_Exclude_Arm$SHORT.NAME))]
                          SHORT.NAME <- SHORT.NAME[which(!is.na(SHORT.NAME))]
                          
                          histologicaldx.match <- 
                            gsub("(, )+", ", ", paste(CTEP.CATEGORY, CTEP.SUBCATEGORY, CTEP.TERM, SHORT.NAME,sep=", "))
                          histologicaldx.match <- 
                            gsub("^, ", "", histologicaldx.match)
                          histologicaldx.match <- 
                            gsub(",[[:blank:]]*$", "", histologicaldx.match)
                          
                          ## Match Disease Exclusion Codes
                          #----------------------------------------------
                          if (isTRUE((length(CTEP.CATEGORY) + length(CTEP.SUBCATEGORY) + 
                                      length(CTEP.TERM) + length(SHORT.NAME)) > 0)) {
                            
                            # Update Match status with exclusion
                            disease_comment <- paste("DISQUALIFIED from ", arm_id, " due to disease exclusion criteria", sep="")
                            DF_patient$NCI_CNV_Variant_Status[pt_rowNo] <- disease_comment
                            
                            cat(paste(patient_id, ": ", DF_patient$NCI_CNV_Variant_Status[pt_rowNo], sep=""))
                            cat("\n", paste("Patient Dx (primary tumor site - histological dx) of ", tumorsite.patient, " ", dx.patient, 
                                            " matched with following CTEP categories: ", 
                                            histologicaldx.match, sep=""),"\n","\n")
                            
                            # Remove excluded ARM_No from DF_Gene_Patient_Variant
                            DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                            
                            DF_Output_Patient_CNV_Variant_int <- 
                              DF_Output_Patient_CNV_Variant_int[DF_Output_Patient_CNV_Variant_int$Arm_Name != arm_id,]
                            
                            # Checkpoint
                            exclusion_continue <- as.logical("FALSE")
                            
                            remove(disease_comment)
                          }
                          remove(CTEP.CATEGORY,CTEP.SUBCATEGORY,CTEP.TERM,SHORT.NAME,histologicaldx.match)
                        }
                        remove(tumorsite.patient,dx.patient)
                        
                        ## Generate output file for candidate MATCH
                        #----------------------------------------------
                        if (isTRUE(exclusion_continue)) {
                          
                          # Corresponding patient INFO - assume single row
                          patient_INFO <- DF_patient[which(DF_patient$CNV_Gene == gene_id &
                                                             DF_patient$var.type == "CNV" &
                                                             DF_patient$var.anno == var.anno_id), c(1:ncol_CNV)]
                          
                          if (nrow(patient_INFO) > 0) {
                            
                            # Matched trial INFO
                            trial_INFO_matched <- 
                              DF_Gene_Patient_Variant[which(DF_Gene_Patient_Variant$Arm_Name == arm_id &
                                                              DF_Gene_Patient_Variant$Gene_Name == gene_id &
                                                              DF_Gene_Patient_Variant$Variant_Type == "CNV" &
                                                              DF_Gene_Patient_Variant$Protein == var.anno_id), ]
                            
                            # Append merged patient INFO with trial INFO to Output file
                            DF_Output_Patient_CNV_Variant_int <- rbind.fill(DF_Output_Patient_CNV_Variant_int,
                                                                            crossing(patient_INFO, trial_INFO_matched))
                            remove(trial_INFO_matched)
                          }
                          remove(patient_INFO)
                        }
                        remove(exclusion_SNV_continue,exclusion_continue,exclusion_Fusion_continue,
                               DF_patient_Fusion,DF_patient_SNVIndel)
                      }
                      remove(arm_id,arm_num)
                    }
                    remove(arm.match)
                    if (isTRUE(exists("var.patient.list"))) {remove(var.patient.list)}
                    if (isTRUE(exists("var.patient"))) {remove(var.patient)}
                    if (isTRUE(exists("var.exclude"))) {remove(var.exclude)}
                  }
                  remove(age_gate)
                }
                remove(pathogenic_id,pathogenicity_gate)
                
              } else {
                DF_patient$NCI_CNV_Variant_Status[which(DF_patient$CNV_Gene == gene_id &
                                                          DF_patient$var.type == "CNV" & 
                                                          DF_patient$var.anno == var.anno_id)] <-
                  "Variant annotation (ie. protein) criteria NOT satisfied"
              }
            }
            remove(var.anno.trial,var.anno.patient,var.anno_id,anno_num)
            
          } else {
            DF_patient$NCI_CNV_Variant_Status[which(DF_patient$CNV_Gene == gene_id &
                                                      DF_patient$var.type == "CNV")] <-
              "Variant Type criteria NOT satisfied"
          }
          remove(DF_Gene_Patient_Variant)
          
        } else {
          DF_patient$NCI_CNV_Variant_Status[which(DF_patient$CNV_Gene == gene_id)] <- 
            "Gene NOT found"
        }
      }
      
      DF_Output_Patient_CNV_Variant <- rbind.fill(DF_Output_Patient_CNV_Variant, DF_Output_Patient_CNV_Variant_int) 
      remove(gene.patient,gene_num,gene_id)
      
      # Write to match results per patient to local computer
      #----------------------------------------------
      write.table(DF_patient, file = paste(tempdir, patient_id, "_CNV.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE)
    }
    remove(DF_Output_Patient_CNV_Variant_int)
  }
  
  ## Remove rows that are all empty or duplicates
  DF_Output_Patient_CNV_Variant <- 
    unique(DF_Output_Patient_CNV_Variant[rowSums(is.na(DF_Output_Patient_CNV_Variant)) != ncol(DF_Output_Patient_CNV_Variant),])
  
  ## Write to match results for positive candidacy local computer
  #----------------------------------------------
  write.table(DF_Output_Patient_CNV_Variant, 
              file = paste(tempdir,"NCIMatch_CNV_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                           dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  remove(DF_patient,genes.Inclusion_Variants,ncol_InclusionVariants,ncol_CNV,patient_id,patient_num,
         DF_Output_Patient_CNV_Variant,CNV.colnames)
  
  if (isTRUE(exists("pt_rowNo"))){remove(pt_rowNo)}
  if (isTRUE(exists("exclusion_comment"))){remove(exclusion_comment)}
}
