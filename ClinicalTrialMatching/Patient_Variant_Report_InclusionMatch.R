# Match by gene > variant type/function > pathogenicity.status \
# > age.group (i.e. adult) > exclusion variant <Gene Protein Variant_Type>
## Search Output: individual patient file with stop point
## Match Output: "NCIMatch_NonHotspot_Matched.tsv"

# HGVSGenomic STAMP entry does not always correspond to CHR:POS_REF_ALT of NCI-MATCH criteria
# Reason: positions are +/- 1 base pair off and/or missing information
# NCI-MATCH Protein entries without HGVS nomenclature will be missed 

cat(paste("Timestamp of NCI-MATCH trial Inclusion Variant matching START: ", Sys.time(), sep=""),"\n","\n")

ncol_STAMP <- as.numeric(ncol(STAMP_DF))
ncol_InclusionVariants <- as.numeric(ncol(Inclusion_Variants))

# Extract gene from OnCore_Biomarker_QC
genes.Inclusion_Variants <- sort(unique(Inclusion_Variants$Gene_Name))

# Output file of matched mutations with clinical trials
DF_Output_Patient_Variant <- data.frame(matrix(NA, ncol = (ncol_STAMP + ncol_InclusionVariants)))
colnames(DF_Output_Patient_Variant) <- append(colnames(STAMP_DF), colnames(Inclusion_Variants))

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
  DF_patient$Patient_Variant_Inclusion_Status <- NA
  
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
        if (isTRUE(var.type_id %in% var.type.trial)) {
          
          ## Parse Protein
          #----------------------------------------------
          protein.trial <- 
            unique(DF_Gene_Patient_Variant$Protein[DF_Gene_Patient_Variant$Variant_Type == var.type_id])
          protein.patient <- 
            unique(DF_patient$VariantHGVSProtein[DF_patient$VariantGene == gene_id & DF_patient$var.type == var.type_id])
          
          # Remove no match Protein from DF_Gene_Patient_Variant
          for (pro_num in 1:length(protein.trial)) {
            if (!(protein.trial[pro_num] %in% protein.patient)) {
              DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Protein !=
                                                                   protein.trial[pro_num],]
            }
          }
          
          # Update list 
          protein.trial <- unique(DF_Gene_Patient_Variant$Protein)
          
          ## Iterate through each protein_id of protein.patient
          #----------------------------------------------
          for (pro_num in 1:length(protein.trial)) {
            protein_id <- protein.patient[pro_num]
            
            ## Match Protein
            #----------------------------------------------
            if (isTRUE(protein_id %in% protein.trial)) {
              
              ## Assess Pathogenicity Status
              #----------------------------------------------
              pathogenicity_gate <- NA
              
              if (isTRUE(pathogenic_FILTER)) {
                if (DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                DF_patient$var.type == var.type_id &
                                                                DF_patient$VariantHGVSProtein == protein_id)] 
                    %in% pathogenic_accepted) {
                  pathogenicity_gate <- as.logical("TRUE")
                  
                } else {
                  pathogenicity_gate <- as.logical("FALSE")
                  DF_patient$Patient_Variant_Inclusion_Status[which(DF_patient$VariantGene == gene_id &
                                                                      DF_patient$var.type == var.type_id &
                                                                      DF_patient$VariantHGVSProtein == protein_id)] <-
                    "Pathogenicity criteria NOT satisfied"
                }
              }
              
              pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                             DF_patient$var.type == var.type_id &
                                                                             DF_patient$VariantHGVSProtein == protein_id)]
              
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
                    DF_patient$Patient_Variant_Inclusion_Status[which(DF_patient$VariantGene == gene_id &
                                                                        DF_patient$var.type == var.type_id &
                                                                        DF_patient$VariantHGVSProtein == protein_id &
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
                                      DF_patient$VariantHGVSProtein == protein_id &
                                      DF_patient$VariantPathogenicityStatus == pathogenic_id)
                  
                  DF_patient$Patient_Variant_Inclusion_Status[pt_rowNo] <- "Candidate trial IDENTIFIED"
                  
                  ## Assess Exclusion Variants
                  #----------------------------------------------
                  arm.match <- unique(DF_Gene_Patient_Variant$Arm_Name
                                      [which(DF_Gene_Patient_Variant$Gene_Name == gene_id &
                                               DF_Gene_Patient_Variant$Variant_Type == var.type_id &
                                               DF_Gene_Patient_Variant$Protein == protein_id)])
                  
                  if (length(arm.match) > 0) {
                    for (arm_num in 1:length(arm.match)) {
                      arm_id <- arm.match[arm_num]
                      
                      # Extract variants for Arm_No in Exclusion_Variants
                      DF_Exclude_Arm <- Exclusion_Variants[Exclusion_Variants$Arm_Name == arm_id,]
                      
                      if (nrow(DF_Exclude_Arm) > 0) {
                        
                        # All variant labels for Arm_No in Exclusion_Variants
                        var.exclude <- sort(unique(paste(DF_Exclude_Arm$Gene_Name, " ", DF_Exclude_Arm$Protein, " (", 
                                                          DF_Exclude_Arm$Variant_Type, ")", sep="")))

                        # All variant labels for patient in DF_patient
                        var.patient <- paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSProtein, " (",
                                             DF_patient$var.type, ")", sep="")
                        var.patient <- sort(unique(var.patient))
                        
                        # Identify overlapping variants
                        var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                        
                        ## Match Exclusion Variants
                        #----------------------------------------------
                        if (length(var.exclude) > 0) {
                          
                          # Update Match status with exclusion
                          exclusion_comment <- paste(DF_patient$Patient_Variant_Inclusion_Status[pt_rowNo],
                                                     "; DISQUALIFIED from ", arm_id, " due to ", 
                                                     paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                     " mutation", sep="")
                          DF_patient$Patient_Variant_Inclusion_Status[pt_rowNo] <- exclusion_comment
                          
                          print(paste(patient_id, ": ", exclusion_comment, sep=""))
                          
                          # Remove excluded ARM_No from DF_Gene_Patient_Variant
                          DF_Gene_Patient_Variant <- DF_Gene_Patient_Variant[DF_Gene_Patient_Variant$Arm_Name != arm_id,]                      
                        }
                      }
                    }  
                  } else {
                    
                    # Candidate match removed due to exclusion at earlier iteration
                    arm.match <- unique(Inclusion_Variants$Arm_Name
                                        [which(Inclusion_Variants$Gene_Name == gene_id &
                                                 Inclusion_Variants$Variant_Type == var.type_id &
                                                 Inclusion_Variants$Protein == protein_id)])
                    
                    for (arm_num in 1:length(arm.match)) {
                      arm_id <- arm.match[arm_num]
                      
                      # Extract variants for Arm_No in Exclusion_Variants
                      DF_Exclude_Arm <- Exclusion_Variants[Exclusion_Variants$Arm_Name == arm_id,]
                      
                      if (nrow(DF_Exclude_Arm) > 0) {
                        
                        # All variant labels for Arm_No in Exclusion_Variants
                        var.exclude <- sort(unique(paste(DF_Exclude_Arm$Gene_Name, " ", DF_Exclude_Arm$Protein, " (", 
                                                          DF_Exclude_Arm$Variant_Type, ")", sep="")))

                        # All variant labels for patient in DF_patient
                        var.patient <- sort(unique(paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSProtein, " (",
                                             DF_patient$var.type, ")", sep="")))
                        
                        # Identify overlapping variants
                        var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                        
                        ## Match Exclusion Variants
                        #----------------------------------------------
                        if (length(var.exclude) > 0) {
                          
                          # Update Match status with exclusion
                          exclusion_comment <- paste(DF_patient$Patient_Variant_Inclusion_Status[pt_rowNo],
                                                     "; DISQUALIFIED from ", arm_id, " due to ", 
                                                     paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                     " mutation", sep="")
                          DF_patient$Patient_Variant_Inclusion_Status[pt_rowNo] <- exclusion_comment
                          
                          print(paste(patient_id, ": ", exclusion_comment, sep=""))
                        }
                      }
                    }
                  }
                }  
              }
              
            } else {
              DF_patient$Patient_Variant_Inclusion_Status[which(DF_patient$VariantGene == gene_id &
                                                                  DF_patient$var.type == var.type_id &
                                                                  DF_patient$VariantHGVSProtein == protein_id)] <-
                "Protein criteria NOT satisfied"
            }
          }
          
        } else {
          DF_patient$Patient_Variant_Inclusion_Status[which(DF_patient$VariantGene == gene_id &
                                                              DF_patient$var.type == var.type_id)] <-
            "Variant Type criteria NOT satisfied"
        }
      }
      
    } else {
      DF_patient$Patient_Variant_Inclusion_Status[which(DF_patient$VariantGene == gene_id)] <- 
        "Gene NOT found"
    }
  }
  
  # Write to match results per patient to local computer
  #----------------------------------------------
  write.table(DF_patient, file = paste(outdir_patient, patient_id, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ## Generate output file for candidate MATCH
  #----------------------------------------------
  if (exists("DF_Gene_Patient_Variant")) {
    if (nrow(DF_Gene_Patient_Variant) > 0) {
      
      for (match_num in 1:nrow(DF_Gene_Patient_Variant)) {
        
        match.gene <- DF_Gene_Patient_Variant$Gene_Name[match_num]
        match.type <- DF_Gene_Patient_Variant$Variant_Type[match_num]
        match.protein <- DF_Gene_Patient_Variant$Protein[match_num]
        
        # Corresponding patient INFO - assume single row
        patient_INFO <- DF_patient[which(DF_patient$VariantGene == match.gene &
                                           DF_patient$var.type == match.type &
                                           DF_patient$VariantHGVSProtein == match.protein), c(1:ncol_STAMP)]
        
        if (nrow(patient_INFO) > 0) {
          # Merge patient INFO with trial INFO
          DF_Output_pre <- cbind(patient_INFO,DF_Gene_Patient_Variant[match_num,])
          
          # Append to Output file
          DF_Output_Patient_Variant <- rbind.fill(DF_Output_Patient_Variant, DF_Output_pre) 
        }
      }
    }
    remove(DF_Gene_Patient_Variant)
  }
}

## Remove rows that are all empty
DF_Output_Patient_Variant <- 
  DF_Output_Patient_Variant[rowSums(is.na(DF_Output_Patient_Variant)) != ncol(DF_Output_Patient_Variant),]

## Write to match results for positive candidacy local computer
#----------------------------------------------
write.table(DF_Output_Patient_Variant, 
            file = paste(outdir_int,"NCIMatch_Variant_Matched_", Patient_Variant_Report_timestamp, "_", filterName, ".tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(DF_Output_pre,DF_patient,patient_INFO,age_gate,arm_id,
       arm_num,arm.match,exclusion_comment,gene_id,gene_num,gene.patient,
       genes.Inclusion_Variants,match_num,match.gene,match.protein,match.type,
       ncol_InclusionVariants,ncol_STAMP,patient_id,patient_num,pro_num,protein_id,
       protein.patient,protein.trial,pt_rowNo,type_num,var.exclude,var.patient,var.type_id,
       var.type.patient,var.type.trial,pathogenicity_gate,pathogenic_id,DF_Exclude_Arm)

cat(paste("Timestamp of NCI-MATCH trial Inclusion Variant matching FINISH: ", Sys.time(), sep=""),"\n","\n")
