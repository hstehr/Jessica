# Match by gene > variant type > Genomic_Region \
# > pathogenicity.status > age.group (i.e. adult) > exclusion variant > exclusion nonhotspot
## Search Output: individual patient file with stop point
## Match Output: "NCI_SNVIndel_Variant_Matched.tsv"

# Code SNVIndel dataframe: lines 12-16

if (isTRUE(NCI_match)) {
  cat(paste("Timestamp of NCI-MATCH trial matching START: ", Sys.time(), sep=""),"\n","\n")
  
  ncol_InclusionVariants <- as.numeric(ncol(Inclusion_Variants))
  ncol_SNVIndel = 13
  SNVIndel.colnames <- c("PatientID","VariantHGVSGenomic","VariantLabel","VariantGene","VariantHGVSCoding",
                         "VariantHGVSProtein","var.type","var.anno","Exon_Number",
                         "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus",
                         "HistologicalDx")
  
  # Extract gene from NCI-MATCH
  genes.Inclusion_Variants <- sort(unique(Inclusion_Variants$Gene_Name))
  
  # Output file of matched mutations with clinical trials
  DF_Output_Patient_Variant <- data.frame(matrix(NA, ncol = (ncol_SNVIndel + ncol_InclusionVariants)))
  colnames(DF_Output_Patient_Variant) <- append(SNVIndel.colnames, colnames(Inclusion_Variants))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) { 
    patient_id <- patient.list[patient_num]
    
    DF_Output_Patient_Variant_int <- data.frame(matrix(NA, ncol = (ncol_SNVIndel + ncol_InclusionVariants)))
    colnames(DF_Output_Patient_Variant_int) <- append(SNVIndel.colnames, colnames(Inclusion_Variants))
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                           header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient) > 0) {
      # Extract unique genes per patient
      gene.patient <- sort(unique(DF_patient$VariantGene))
      
      # Clinical match result
      DF_patient$NCI_SNVIndel_Variant_Status <- NA
      
      ## Iterate through each gene_id of gene.patient
      #----------------------------------------------
      for (gene_num in 1:length(gene.patient)) {
        gene_id <- gene.patient[gene_num]
        
        ## Match Biomarker_GeneName
        #----------------------------------------------
        if (gene_id %in% genes.Inclusion_Variants) {
          
          # Extract gene_id match from Inclusion_Variants
          DF_Gene_Patient_Variant <- data.frame(Inclusion_Variants[which(Inclusion_Variants$Gene_Name == gene_id),])
          
          # Remove trailing whitespace from HGVS nomenclature
          DF_Gene_Patient_Variant$Protein <- gsub("[[:space:]]*$", "", DF_Gene_Patient_Variant$Protein)
          
          ## Parse Variant_Type
          #----------------------------------------------
          var.type.trial <- unique(DF_Gene_Patient_Variant$Variant_Type)
          var.type.patient <- unique(DF_patient$var.type[DF_patient$VariantGene == gene_id])
          
          # Remove no match Variant_Type from DF_Gene_Patient_Variant
          for (type_num in 1:length(var.type.trial)) {
            if (isTRUE(!(var.type.trial[type_num] %in% var.type.patient))) {
              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Variant_Type !=
                                                                   var.type.trial[type_num],]
            }
          }
          
          # Update list 
          var.type.trial <- unique(DF_Gene_Patient_Variant$Variant_Type)
          
          ## Iterate through each var.type_id of var.type.patient
          #----------------------------------------------
          for (type_num in 1:length(var.type.patient)) {
            var.type_id <- var.type.patient[type_num]
            
            ## Match Variant_Type
            #----------------------------------------------
            if (isTRUE(var.type_id %in% var.type.trial & nrow(DF_Gene_Patient_Variant) > 0)) {
              
              ## Parse GenomicRegion
              #----------------------------------------------
              genomic.trial <- 
                unique(DF_Gene_Patient_Variant$Genomic[DF_Gene_Patient_Variant$Variant_Type == var.type_id])
              genomic.patient <- 
                unique(DF_patient$VariantHGVSGenomic[DF_patient$VariantGene == gene_id & DF_patient$var.type == var.type_id])
              
              # Remove no match Genomic Region from DF_Gene_Patient_Variant
              if (isTRUE(length(genomic.trial) > 0)) {
                for (pro_num in 1:length(genomic.trial)) {
                  if (!(genomic.trial[pro_num] %in% genomic.patient)) {
                    DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Genomic !=
                                                                         genomic.trial[pro_num],]
                  }
                }
                
                # Update list 
                genomic.trial <- unique(DF_Gene_Patient_Variant$Genomic)
              }
              
              ## Iterate through each genomic_id of genomic.patient
              #----------------------------------------------
              for (gen_num in 1:length(genomic.patient)) {
                genomic_id <- genomic.patient[gen_num]
                
                ## Match Genomic Region
                #----------------------------------------------
                if (isTRUE(genomic_id %in% genomic.trial & nrow(DF_Gene_Patient_Variant) > 0)) {
                  
                  ## Assess Pathogenicity Status
                  #----------------------------------------------
                  pathogenicity_gate <- NA
                  
                  if (isTRUE(pathogenic_FILTER)) {
                    if (DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                    DF_patient$var.type == var.type_id &
                                                                    DF_patient$VariantHGVSGenomic == genomic_id)] 
                        %in% pathogenic_accepted) {
                      pathogenicity_gate <- as.logical("TRUE")
                      
                    } else {
                      pathogenicity_gate <- as.logical("FALSE")
                      DF_patient$NCI_SNVIndel_Variant_Status[which(DF_patient$VariantGene == gene_id &
                                                                     DF_patient$var.type == var.type_id &
                                                                     DF_patient$VariantHGVSGenomic == genomic_id)] <-
                        "Pathogenicity criteria NOT satisfied"
                    }
                  }
                  
                  pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                                 DF_patient$var.type == var.type_id &
                                                                                 DF_patient$VariantHGVSGenomic == genomic_id)]
                  
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
                        DF_patient$NCI_SNVIndel_Variant_Status[which(DF_patient$VariantGene == gene_id &
                                                                       DF_patient$var.type == var.type_id &
                                                                       DF_patient$VariantHGVSGenomic == genomic_id &
                                                                       DF_patient$VariantPathogenicityStatus == pathogenic_id)] <-
                          "Age criteria NOT satisfied"
                      }
                    }
                    
                    ## Match Age.Group criteria
                    #----------------------------------------------
                    if (isTRUE(age_gate == TRUE | is.na(age_gate))) {
                      
                      # Corresponding row in patient file
                      pt_rowNo <- which(DF_patient$VariantGene == gene_id &
                                          DF_patient$var.type == var.type_id &
                                          DF_patient$VariantHGVSGenomic == genomic_id &
                                          DF_patient$VariantPathogenicityStatus == pathogenic_id)
                      
                      arm.match <- unique(DF_Gene_Patient_Variant$Arm_Name
                                          [which(DF_Gene_Patient_Variant$Gene_Name == gene_id &
                                                   DF_Gene_Patient_Variant$Variant_Type == var.type_id &
                                                   DF_Gene_Patient_Variant$Genomic == genomic_id)])
                      
                      DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo] <- paste("Candidate trial IDENTIFIED: ", 
                                                                                paste(arm.match, collapse = ", "), sep="")
                      
                      ## Assess Exclusions
                      #----------------------------------------------
                      if (length(arm.match) > 0) {
                        
                        for (arm_num in 1:length(arm.match)) {
                          
                          arm_id <- arm.match[arm_num]
                          exclusion_continue <- as.logical("TRUE")
                          
                          # Import CNV entries per patient
                          #----------------------------------------------
                          file_name = paste(tempdir, patient_id, "_CNV.tsv", sep="")
                          
                          if (file.exists(file_name)) {
                            DF_patient_CNV <- read.csv(file = file_name, header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                            
                            if (nrow(DF_patient_CNV) > 0) { exclusion_CNV_continue <- as.logical("TRUE")  
                            } else { exclusion_CNV_continue <- as.logical("FALSE")}
                            
                          } else { exclusion_CNV_continue <- as.logical("FALSE")}
                          
                          # Import Fusion entries per patient
                          #----------------------------------------------
                          file_name = paste(tempdir, patient_id, "_Fusion.tsv", sep="")
                          
                          if (file.exists(file_name)) {
                            DF_patient_Fusion <- read.csv(file = file_name, header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                            
                            if (nrow(DF_patient_Fusion) > 0) {exclusion_Fusion_continue <- as.logical("TRUE")  
                            } else {exclusion_Fusion_continue <- as.logical("FALSE")}
                            
                          } else {exclusion_Fusion_continue <- as.logical("FALSE")}
                          
                          ## Assess Exclusion Variants
                          #----------------------------------------------
                          # Extract variants for Arm_No in Exclusion_Variants
                          DF_Exclude_Arm <- Exclusion_Variants[Exclusion_Variants$Arm_Name == arm_id,]
                          
                          if (nrow(DF_Exclude_Arm) > 0) {
                            
                            ## Assess Exclusion Variants: SNVs
                            #----------------------------------------------
                            DF_Exclude_Arm_SNV <- DF_Exclude_Arm[which(DF_Exclude_Arm$Variant_Type == "SNV"),]
                            
                            if (nrow(DF_Exclude_Arm_SNV) > 0) {
                              
                              # All variant labels for Arm_No in Exclusion_Variants
                              var.exclude <- sort(unique(paste(DF_Exclude_Arm_SNV$Gene_Name, " ", DF_Exclude_Arm_SNV$Genomic, " (", 
                                                               DF_Exclude_Arm_SNV$Variant_Type, ")", sep="")))
                              
                              # All variant labels for patient in DF_patient
                              var.patient <- sort(unique(paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSGenomic, " (",
                                                               DF_patient$var.type, ")", sep="")))
                              
                              # Identify overlapping variants
                              var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                              
                              ## Match Exclusion Variants
                              #----------------------------------------------
                              if (length(var.exclude) > 0) {
                                
                                # Update Match status with exclusion
                                exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                           paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                           " mutation", sep="")
                                DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo] <- exclusion_comment
                                
                                cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                
                                # Remove excluded ARM_No from DF_Gene_Patient_Variant
                                DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                                
                                DF_Output_Patient_Variant_int <- DF_Output_Patient_Variant_int[DF_Output_Patient_Variant_int$Arm_Name != arm_id,]
                                
                                # Checkpoint
                                exclusion_continue <- as.logical("FALSE")
                                
                                remove(exclusion_comment)
                              }
                            }
                            
                            ## Assess Exclusion Variants: CNVs
                            #----------------------------------------------
                            DF_Exclude_Arm_CNV <- DF_Exclude_Arm[which(DF_Exclude_Arm$Variant_Type == "CNV"), ]
                            DF_Exclude_Arm_CNV$Protein <- 
                              tolower(gsub("(^[[:upper:]]{2,}.*[[:blank:]]+)([[:alpha:]]+)","\\2",DF_Exclude_Arm_CNV$Protein))
                            
                            if (isTRUE(nrow(DF_Exclude_Arm_CNV) > 0 & exclusion_continue & exclusion_CNV_continue)) {
                              
                              # Structure DF
                              DF_patient_CNV$var.anno <- gsub("AMP","amplification",DF_patient_CNV$var.anno)
                              DF_patient_CNV$var.anno <- gsub("DEL","deletion",DF_patient_CNV$var.anno)
                              
                              # All variant labels for Arm_No in Exclusion_Variants
                              var.exclude <- sort(unique(paste(tolower(DF_Exclude_Arm_CNV$Gene_Name), " (", 
                                                               DF_Exclude_Arm_CNV$Protein, ")", sep="")))
                              
                              # All variant labels for patient in DF_patient
                              var.patient <- sort(unique(paste(tolower(DF_patient_CNV$CNV_Gene), " (",
                                                               DF_patient_CNV$var.anno, ")", sep="")))
                              
                              # Identify overlapping variants
                              var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                              
                              ## Match Exclusion Variants
                              #----------------------------------------------
                              if (length(var.exclude) > 0) {
                                
                                # Update Match status with exclusion
                                exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                           paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                           " CNV", sep="")
                                DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo] <- exclusion_comment
                                
                                cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                
                                # Remove excluded ARM_No from DF_Gene_Patient_Variant
                                DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                                
                                DF_Output_Patient_Variant_int <- DF_Output_Patient_Variant_int[DF_Output_Patient_Variant_int$Arm_Name != arm_id,]
                                
                                # Checkpoint
                                exclusion_continue <- as.logical("FALSE")
                                
                                remove(exclusion_comment)
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
                                DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo] <- exclusion_comment
                                
                                cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                
                                # Remove excluded ARM_No from DF_Gene_Patient_Variant
                                DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                                
                                DF_Output_Patient_Variant_int <- DF_Output_Patient_Variant_int[DF_Output_Patient_Variant_int$Arm_Name != arm_id,]
                                
                                # Checkpoint
                                exclusion_continue <- as.logical("FALSE")
                                
                                remove(exclusion_comment)
                              }
                            }
                            remove(DF_Exclude_Arm_SNV, DF_Exclude_Arm_CNV, DF_Exclude_Arm_Fusion,
                                   var.patient,var.exclude)
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
                              DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo] <- disease_comment
                              
                              cat(paste(patient_id, ": ", DF_patient$NCI_SNVIndel_Variant_Status[pt_rowNo], sep=""))
                              cat("\n", paste("Patient Dx (primary tumor site - histological dx) of ", tumorsite.patient, " ", dx.patient, 
                                              " matched with following CTEP categories: ", 
                                              histologicaldx.match, sep=""),"\n","\n")
                              
                              # Remove excluded ARM_No from DF_Gene_Patient_Variant
                              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]
                              
                              DF_Output_Patient_Variant_int <- DF_Output_Patient_Variant_int[DF_Output_Patient_Variant_int$Arm_Name != arm_id,]
                              
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
                            patient_INFO <- 
                              DF_patient[which(DF_patient$VariantGene == gene_id &
                                                 DF_patient$var.type == var.type_id &
                                                 DF_patient$VariantHGVSGenomic == genomic_id), SNVIndel.colnames]
                            
                            if (nrow(patient_INFO) > 0) {
                              
                              # Matched trial INFO
                              trial_INFO_matched <- 
                                DF_Gene_Patient_Variant[which(DF_Gene_Patient_Variant$Arm_Name == arm_id &
                                                                DF_Gene_Patient_Variant$Gene_Name == gene_id &
                                                                DF_Gene_Patient_Variant$Variant_Type == var.type_id &
                                                                DF_Gene_Patient_Variant$Genomic == genomic_id),]
                              
                              # Append merged patient INFO with trial INFO to Output file
                              DF_Output_Patient_Variant_int <- rbind.fill(DF_Output_Patient_Variant_int,
                                                                          crossing(patient_INFO, trial_INFO_matched))
                              remove(trial_INFO_matched)
                            }
                            remove(patient_INFO)
                          }
                          if (exists("DF_patient_CNV")) {remove(DF_patient_CNV)}
                          if (exists("DF_patient_Fusion")) {remove(DF_patient_Fusion)}
                          remove(exclusion_CNV_continue,exclusion_Fusion_continue,exclusion_continue,file_name)
                        }
                        remove(arm_id,arm_num)
                      }
                      remove(arm.match,pt_rowNo)
                    }
                    remove(age_gate)
                  }
                  remove(pathogenic_id,pathogenicity_gate)
                  
                } else {
                  DF_patient$NCI_SNVIndel_Variant_Status[which(DF_patient$VariantGene == gene_id &
                                                                 DF_patient$var.type == var.type_id &
                                                                 DF_patient$VariantHGVSGenomic == genomic_id)] <-
                    "Genomic region criteria NOT satisfied"
                }
                remove(gen_num)
              }
              remove(pro_num,genomic.patient,genomic.trial,genomic_id)
              
            } else {
              DF_patient$NCI_SNVIndel_Variant_Status[which(DF_patient$VariantGene == gene_id &
                                                             DF_patient$var.type == var.type_id)] <-
                "Variant Type criteria NOT satisfied"
            }
          }
          remove(var.type.trial,var.type.patient,var.type_id,type_num,DF_Gene_Patient_Variant)
          
        } else {
          DF_patient$NCI_SNVIndel_Variant_Status[which(DF_patient$VariantGene == gene_id)] <- 
            "Gene NOT found"
        }
      }
      
      DF_Output_Patient_Variant <- rbind.fill(DF_Output_Patient_Variant, DF_Output_Patient_Variant_int) 
      remove(DF_Output_Patient_Variant_int)
      
      # Write to match results per patient to local computer
      #----------------------------------------------
      write.table(DF_patient, file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE)
    }
  }
  
  ## Remove rows that are all empty or duplicates
  DF_Output_Patient_Variant <- 
    unique(DF_Output_Patient_Variant[rowSums(is.na(DF_Output_Patient_Variant)) != ncol(DF_Output_Patient_Variant),])
  
  ## Write to match results for positive candidacy local computer
  #----------------------------------------------
  write.table(DF_Output_Patient_Variant, 
              file = paste(tempdir,"NCI_SNVIndel_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                           dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)
  
  remove(DF_patient,gene.patient,gene_id,gene_num,genes.Inclusion_Variants,
         ncol_InclusionVariants,ncol_SNVIndel,patient_id,patient_num,
         DF_Output_Patient_Variant,SNVIndel.colnames)
}
