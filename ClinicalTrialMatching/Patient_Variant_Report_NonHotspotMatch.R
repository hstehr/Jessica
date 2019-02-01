# Match by gene > variant type/function > pathogenicity.status \
# > age.group (i.e. adult) > exclusion variant <Gene Protein Variant_Type>
## Search Output: individual patient file with stop point
## Match Output: "NCIMatch_NonHotspot_Matched.tsv"

cat(paste("Timestamp of NCI-MATCH trial Inclusion NonHotspot matching START: ", Sys.time(), sep=""),"\n","\n")

ncol_STAMP <- as.numeric(ncol(STAMP_DF))
ncol_InclusionNonHotspot <- as.numeric(ncol(Inclusion_NonHotspot_Rules))

# Extract gene from OnCore_Biomarker_QC
genes.Inclusion_NonHotspot <- sort(unique(Inclusion_NonHotspot_Rules$Gene_Name))

# Output file of matched mutations with clinical trials
DF_Output_Patient_NonHotspot <- data.frame(matrix(NA, ncol = (ncol_STAMP + ncol_InclusionNonHotspot)))
colnames(DF_Output_Patient_NonHotspot) <- append(colnames(STAMP_DF), colnames(Inclusion_NonHotspot_Rules))

## Iterate through each patient_id of patient.list
#----------------------------------------------
# patient_num = which(patient.list == "TRF-3932")
for (patient_num in 1:length(patient.list)) {
  patient_id <- patient.list[patient_num]
  
  # Import STAMP entries per patient
  DF_patient <- read.csv(file = paste(outdir_patient, patient_id, ".tsv", sep=""),
                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
  
  # Extract unique genes per patient
  gene.patient <- sort(unique(DF_patient$VariantGene))
  
  # Clinical match result
  DF_patient$Patient_Variant_NonHotspot_Status <- NA
  
  # Generate output per unique gene
  for (gene_num in 1:length(gene.patient)) {
    gene_id <- gene.patient[gene_num]
    
    ## Extract matches from DF_Inclusion_Variants = DF_Gene_Patient_NonHotspot
    #----------------------------------------------
    if (gene_id %in% genes.Inclusion_NonHotspot) {
      
      # Extract gene_id match from Inclusion_Variants
      DF_Gene_Patient_NonHotspot <- 
        data.frame(Inclusion_NonHotspot_Rules[which(Inclusion_NonHotspot_Rules$Gene_Name == gene_id),])
      
      ## Parse Variant_Type
      #----------------------------------------------
      var.type.trial <- unique(DF_Gene_Patient_NonHotspot$Function[!is.na(DF_Gene_Patient_NonHotspot$Function)])
      var.type.patient <- unique(DF_patient$var.type[DF_patient$VariantGene == gene_id])
      
      ## Iterate through each bio.cond_id of var.type.patient
      #----------------------------------------------
      for (type_num in 1:length(var.type.patient)) {
        var.type_id <- var.type.patient[type_num]
        
        # Translate terminology to common syntax
        var.type_id_syn <- var.type_id
        var.type_id_syn <- gsub("Deletion", "nonframeshiftDeletion", var.type_id_syn)
        var.type_id_syn <- gsub("Insertion", "nonframeshiftInsertion", var.type_id_syn)
        var.type_id_syn <- gsub("Delins", "nonframeshiftDeletion,nonframeshiftInsertion", var.type_id_syn)
        var.type_id_syn <- unlist(strsplit(paste(sub("/.*", "", var.type_id_syn), collapse = " ,"), ","))
        var.type_id_syn <- var.type_id_syn[which(var.type_id_syn %in% var.type.trial)]
        
        ## Match Variant_Type
        #----------------------------------------------
        if (isTRUE(var.type_id_syn %in% var.type.trial)) {
          
          ## Assess Pathogenicity Status
          #----------------------------------------------
          pathogenicity_gate <- NA
          
          if (isTRUE(pathogenic_FILTER)) {
            if (DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                            DF_patient$var.type == var.type_id)] 
                %in% pathogenic_accepted) {
              pathogenicity_gate <- as.logical("TRUE")
              
            } else {
              pathogenicity_gate <- as.logical("FALSE")
              DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                   DF_patient$var.type == var.type_id)] <-
                "Pathogenicity criteria NOT satisfied"
            }
          }
          
          pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                         DF_patient$var.type == var.type_id)]
          
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
                DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                     DF_patient$var.type == var.type_id &
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
                                  DF_patient$VariantPathogenicityStatus == pathogenic_id)
              
              DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- "Candidate trial IDENTIFIED"
              
              ## Assess Exclusion Variants
              #----------------------------------------------
              arm.match <- unique(DF_Gene_Patient_NonHotspot$Arm_Name
                                  [which(DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                           DF_Gene_Patient_NonHotspot$Function == var.type_id_syn)])
              
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
                    var.patient <- sort(unique(paste(DF_patient$VariantGene, " ", DF_patient$VariantHGVSProtein, " (",
                                                     DF_patient$var.type, ")", sep="")))
                    
                    # Identify overlapping variants
                    var.exclude <- var.exclude[which(var.exclude %in% var.patient)]
                    
                    ## Match Exclusion Variants
                    #----------------------------------------------
                    if (length(var.exclude) > 0) {
                      
                      # Update Match status with exclusion
                      exclusion_comment <- paste(DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo],
                                                 "; DISQUALIFIED from ", arm_id, " due to ", 
                                                 paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                 " mutation", sep="")
                      DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                      
                      print(paste(patient_id, ": ", exclusion_comment, sep=""))
                      
                      # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                      DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]                      
                    }
                  }
                }  
              } else {
                
                # Candidate match removed due to exclusion at earlier iteration
                arm.match <- unique(Inclusion_NonHotspot_Rules$Arm_Name
                                    [which(Inclusion_NonHotspot_Rules$Gene_Name == gene_id &
                                             Inclusion_NonHotspot_Rules$Function == var.type_id_syn)])
                
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
                      exclusion_comment <- paste(DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo],
                                                 "; DISQUALIFIED from ", arm_id, " due to ", 
                                                 paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                 " mutation", sep="")
                      DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                      
                      print(paste(patient_id, ": ", exclusion_comment, sep=""))
                    }
                  }
                }
              }
            }  
          }
          
        }  else {
          DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                               DF_patient$var.type == var.type_id)] <-
            "Variant Type criteria NOT satisfied"
        }
      }
      
    } else {
      DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id)] <- 
        "Gene NOT found"
    }
  }
  
  # Write to match results per patient to local computer
  #----------------------------------------------
  write.table(DF_patient, file = paste(outdir_patient, patient_id, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ## Generate output file for candidate MATCH
  #----------------------------------------------
  if (exists("DF_Gene_Patient_NonHotspot")) {
    if (nrow(DF_Gene_Patient_NonHotspot) > 0) {
      
      for (match_num in 1:nrow(DF_Gene_Patient_NonHotspot)) {
        
        match.gene <- DF_Gene_Patient_NonHotspot$Gene_Name[match_num]
        match.type_syn <- DF_Gene_Patient_NonHotspot$Function[match_num]
        
        # Translate terminology to common syntax
        match.type <- match.type_syn
        match.type <- gsub("nonframeshiftDeletion","Deletion,Delins", match.type)
        match.type <- gsub("nonframeshiftInsertion", "Insertion,Delins",match.type)
        match.type <- unlist(strsplit(paste(sub("/.*", "", match.type), collapse = " ,"), ","))
        
        # Corresponding patient INFO - assume single row
        patient_INFO <- DF_patient[which(DF_patient$VariantGene == match.gene &
                                           DF_patient$var.type %in% match.type), c(1:ncol_STAMP)]
        
        if (nrow(patient_INFO) > 0) {
          # Merge patient INFO with trial INFO
          DF_Output_pre <- cbind(patient_INFO,DF_Gene_Patient_NonHotspot[match_num,])
          
          # Append to Output file
          DF_Output_Patient_NonHotspot <- rbind.fill(DF_Output_Patient_NonHotspot, DF_Output_pre)      
        }
      }
    }
    remove(DF_Gene_Patient_NonHotspot)
  }
}

## Remove rows that are all empty
DF_Output_Patient_NonHotspot <- 
  DF_Output_Patient_NonHotspot[rowSums(is.na(DF_Output_Patient_NonHotspot)) != ncol(DF_Output_Patient_NonHotspot),]

## Write to local computer
#----------------------------------------------
write.table(DF_Output_Patient_NonHotspot, 
            file = paste(outdir_int,"NCIMatch_NonHotspot_Matched_", Patient_Variant_Report_timestamp, "_", filterName, ".tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(DF_Exclude_Arm,DF_Output_pre,DF_patient,patient_INFO,age_gate,arm_id,arm_num,arm.match,
       exclusion_comment,gene_id,gene_num,gene.patient,genes.Inclusion_NonHotspot,match_num,
       match.gene,match.type,match.type_syn,ncol_InclusionNonHotspot,ncol_STAMP,pathogenic_id,
       pathogenicity_gate,patient_id,patient_num,pt_rowNo,type_num,var.exclude,var.patient,
       var.type_id,var.type_id_syn,var.type.patient, var.type.trial)

cat(paste("Timestamp of NCI-MATCH trial Inclusion NonHotspot matching FINISH: ", Sys.time(), sep=""),"\n","\n")
