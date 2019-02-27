# Match by gene > variant type/function > pathogenicity.status > age.group (i.e. adult) \
# > exclusion variant <Gene Protein Variant_Type> \
# > exclusion nonhotspot <Gene Exon_Number Variant_Type>
## Search Output: individual patient file with stop point
## Match Output: "NCIMatch_NonHotspot_Matched.tsv"

if (isTRUE(NCI_match)) {
  
  ncol_STAMP <- as.numeric(ncol(STAMP_DF))
  ncol_InclusionNonHotspot <- as.numeric(ncol(Inclusion_NonHotspot_Rules))
  
  # Extract gene from OnCore_Biomarker_QC
  genes.Inclusion_NonHotspot <- sort(unique(Inclusion_NonHotspot_Rules$Gene_Name))
  
  # Output file of matched mutations with clinical trials
  DF_Output_Patient_NonHotspot <- data.frame(matrix(NA, ncol = (ncol_STAMP + ncol_InclusionNonHotspot)))
  colnames(DF_Output_Patient_NonHotspot) <- append(colnames(STAMP_DF), colnames(Inclusion_NonHotspot_Rules))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    DF_Output_Patient_NonHotspot_int <- data.frame(matrix(NA, ncol = (ncol_STAMP + ncol_InclusionNonHotspot)))
    colnames(DF_Output_Patient_NonHotspot_int) <- append(colnames(STAMP_DF), colnames(Inclusion_NonHotspot_Rules))
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, ".tsv", sep=""),
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
            
            ## Parse Exon Number
            #----------------------------------------------
            exon.trial <- 
              unique(DF_Gene_Patient_NonHotspot$Exon[DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                       DF_Gene_Patient_NonHotspot$Function == var.type_id_syn])
            exon.patient <- 
              unique(DF_patient$Exon_Number[DF_patient$VariantGene == gene_id & DF_patient$var.type == var.type_id])
            
            # Remove no match Exon Number from DF_Gene_Patient_NonHotspot
            if (isTRUE(length(exon.trial) > 0)) {
              for (exon_num in 1:length(exon.trial)) {
                if (!(exon.trial[exon_num] %in% exon.patient)) {
                  DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Exon !=
                                                                             exon.trial[exon_num],]
                }
              }
              
              # Update list 
              exon.trial <- unique(DF_Gene_Patient_NonHotspot$Exon)
            }
            
            ## Iterate through each exon_id of exon.patient
            #----------------------------------------------
            for (exon_num in 1:length(exon.patient)) {
              exon_id <- exon.patient[exon_num]
              
              ## Match Protein
              #----------------------------------------------
              if (isTRUE(exon_id %in% exon.trial & nrow(DF_Gene_Patient_NonHotspot) > 0)) {
                
                ## Assess Pathogenicity Status
                #----------------------------------------------
                pathogenicity_gate <- NA
                
                if (isTRUE(pathogenic_FILTER)) {
                  if (DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                  DF_patient$var.type == var.type_id &
                                                                  DF_patient$Exon_Number == exon_id)] 
                      %in% pathogenic_accepted) {
                    pathogenicity_gate <- as.logical("TRUE")
                    
                  } else {
                    pathogenicity_gate <- as.logical("FALSE")
                    DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                         DF_patient$var.type == var.type_id &
                                                                         DF_patient$Exon_Number == exon_id)] <-
                      "Pathogenicity criteria NOT satisfied"
                  }
                }
                
                pathogenic_id <- DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                               DF_patient$var.type == var.type_id &
                                                                               DF_patient$Exon_Number == exon_id)]
                
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
                                                                           DF_patient$Exon_Number == exon_id &
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
                                        DF_patient$Exon_Number == exon_id &
                                        DF_patient$VariantPathogenicityStatus == pathogenic_id)
                    
                    arm.match <- unique(DF_Gene_Patient_NonHotspot$Arm_Name
                                        [which(DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                 DF_Gene_Patient_NonHotspot$Function == var.type_id_syn &
                                                 DF_Gene_Patient_NonHotspot$Exon == exon_id)])
                    
                    DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- paste("Candidate trial IDENTIFIED: ", 
                                                                                    paste(arm.match, collapse = ", "), sep="")
                    
                    ## Assess Exclusion Variants
                    #----------------------------------------------
                    if (length(arm.match) > 0) {
                      for (arm_num in 1:length(arm.match)) {
                        arm_id <- arm.match[arm_num]
                        exclusion_continue <- as.logical("TRUE")
                        
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
                            exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                       paste(sub("/.*", "", var.exclude), collapse = " & "), 
                                                       " mutation", sep="")
                            DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                            
                            cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                            
                            # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                            DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                            
                            DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                            
                            # Checkpoint
                            exclusion_continue <- as.logical("FALSE")
                          }
                        }
                        
                        ## Assess Exclusion NonHotspot
                        #----------------------------------------------
                        if (isTRUE(exclusion_continue)) {
                          
                          # Extract variants for Arm_No in Exclusion_NonHotspot_Rules
                          DF_Exclude_Arm <- 
                            Exclusion_NonHotspot_Rules[Exclusion_NonHotspot_Rules$Arm_Name == arm_id,]
                          
                          if (nrow(DF_Exclude_Arm) > 0) {
                            
                            # All variant gene-exon for Arm_No in Exclusion_NonHotspot_Rules
                            var.exclude <- unique(data.frame(Gene_Name=DF_Exclude_Arm$Gene_Name, 
                                                             Exon=DF_Exclude_Arm$Exon,
                                                             Function=DF_Exclude_Arm$Function, stringsAsFactors = FALSE))
                            # Translate terminology to common syntax
                            for (row_No in nrow(var.exclude)) {
                              var.exclude$var.type[row_No] = gsub("nonframeshiftDeletion", "Deletion,Delins", var.exclude$Function[row_No])
                              var.exclude$var.type[row_No] = gsub("nonframeshiftInsertion", "Insertion,Delins", var.exclude$Function[row_No])
                            }
                            # Strip var.type into components
                            var.exclude <- cSplit(var.exclude, "var.type", ",", 
                                                  stripWhite = TRUE, type.convert="as.character",
                                                  drop = FALSE)
                            ## Convert to long format 
                            var.exclude <- melt(var.exclude, id.vars=c("Gene_Name","Exon","Function"))
                            var.exclude <- var.exclude[which(var.exclude$variable != "var.type" & !is.na(var.exclude$value)), 
                                                       c("Gene_Name","Exon","Function","value")]
                            
                            var.exclude.list <- sort(unique(paste(var.exclude$Gene_Name, " ", var.exclude$Exon, " (", 
                                                                  var.exclude$value, ")", sep="")))
                            
                            # All variant gene-exon for patient in DF_patient
                            var.patient <- unique(data.frame(VariantGene=DF_patient$VariantGene,
                                                             Exon_Number=DF_patient$Exon_Number,
                                                             var.type=DF_patient$var.type, stringsAsFactors = FALSE))
                            
                            var.patient.list <- sort(unique(paste(var.patient$VariantGene, " ", var.patient$Exon_Number, " (",
                                                                  var.patient$var.type, ")", sep="")))
                            
                            # Identify overlapping variants
                            var.exclude.list <- var.exclude.list[which(var.exclude.list %in% var.patient.list)]
                            
                            ## Match Exclusion Variants
                            #----------------------------------------------
                            if (length(var.exclude.list) > 0) {
                              
                              # Update Match status with exclusion
                              exclusion_comment <- paste("DISQUALIFIED from ", arm_id, " due to ", 
                                                         paste(sub("/.*", "", var.exclude.list), collapse = " & "), 
                                                         " mutation", sep="")
                              DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                              
                              cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                              
                              # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                              DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                              
                              DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                              
                              # Checkpoint
                              exclusion_continue <- as.logical("FALSE")
                            }
                          }
                        }
                        
                        ## Assess Disease Exclusion Codes
                        #----------------------------------------------
                        dx.patient = unique(DF_patient$HistologicalDx)
                        tumorsite.patient = unique(DF_patient$PrimaryTumorSite)
                        
                        if (isTRUE(!is.na(dx.patient) & exclusion_continue & disease.code_FILTER)) {
                          
                          # Extract variants for Arm_No in Disease Exclusion Codes
                          DF_Disease_Exclude_Arm <- Disease_Exclusion_Codes[Disease_Exclusion_Codes$Arm_Name == arm_id,]
                          DF_Disease_Exclude_patient <-
                            HistologicalDxCategory[which(HistologicalDxCategory$histologicalDiagnosis == dx.patient &
                                                           HistologicalDxCategory$primaryTumorSite == tumorsite.patient),]
                          
                          # Identify overlapping classifications 
                          CTEP.CATEGORY <- DF_Disease_Exclude_patient$CTEP.CATEGORY[which(DF_Disease_Exclude_patient$CTEP.CATEGORY %in% tolower(DF_Disease_Exclude_Arm$CTEP.CATEGORY))]
                          CTEP.CATEGORY <- CTEP.CATEGORY[which(!is.na(CTEP.CATEGORY))]
                          
                          CTEP.SUBCATEGORY <- DF_Disease_Exclude_patient$CTEP.SUBCATEGORY[which(DF_Disease_Exclude_patient$CTEP.SUBCATEGORY %in% tolower(DF_Disease_Exclude_Arm$CTEP.SUBCATEGORY))]
                          CTEP.SUBCATEGORY <- CTEP.SUBCATEGORY[which(!is.na(CTEP.SUBCATEGORY))]
                          
                          CTEP.TERM <- DF_Disease_Exclude_patient$CTEP.TERM[which(DF_Disease_Exclude_patient$CTEP.TERM %in% tolower(DF_Disease_Exclude_Arm$CTEP.TERM))]
                          CTEP.TERM <- CTEP.TERM[which(!is.na(CTEP.TERM))]
                          
                          SHORT.NAME <- DF_Disease_Exclude_patient$SHORT.NAME[which(DF_Disease_Exclude_patient$SHORT.NAME %in% tolower(DF_Disease_Exclude_Arm$SHORT.NAME))]
                          SHORT.NAME <- SHORT.NAME[which(!is.na(SHORT.NAME))]
                          
                          histologicaldx.match <- gsub("(, )+", ", ", paste(CTEP.CATEGORY, CTEP.SUBCATEGORY, CTEP.TERM, SHORT.NAME,sep=", "))
                          histologicaldx.match <- gsub("^, ", "", histologicaldx.match)
                          histologicaldx.match <- gsub(",[[:blank:]]*$", "", histologicaldx.match)
                          
                          ## Match Disease Exclusion Codes
                          #----------------------------------------------
                          if (isTRUE((length(CTEP.CATEGORY) + length(CTEP.SUBCATEGORY) + 
                                      length(CTEP.TERM) + length(SHORT.NAME)) > 0)) {
                            
                            # Update Match status with exclusion
                            disease_comment <- paste("DISQUALIFIED from ", arm_id, " due to disease exclusion criteria", sep="")
                            DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo] <- disease_comment
                            
                            cat(paste(patient_id, ": ", DF_patient$Patient_Variant_NonHotspot_Status[pt_rowNo], sep=""))
                            cat("\n", paste("Patient Dx (primary tumor site - histological dx) of ", tumorsite.patient, " ", dx.patient, 
                                            " matched with following CTEP categories: ", 
                                            histologicaldx.match, sep=""),"\n","\n")
                            
                            # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                            DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                            
                            DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                            
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
                          patient_INFO <- DF_patient[which(DF_patient$VariantGene == gene_id &
                                                             DF_patient$var.type == var.type_id &
                                                             DF_patient$Exon_Number == exon_id), c(1:ncol_STAMP)]
                          
                          if (nrow(patient_INFO) > 0) {
                            
                            # Matched trial INFO
                            trial_INFO_matched <- DF_Gene_Patient_NonHotspot[which(DF_Gene_Patient_NonHotspot$Arm_Name == arm_id &
                                                                                     DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                                                     DF_Gene_Patient_NonHotspot$Function == var.type_id_syn &
                                                                                     DF_Gene_Patient_NonHotspot$Exon == exon_id),]
                            
                            # Append merged patient INFO with trial INFO to Output file
                            for (row_No in 1:nrow(trial_INFO_matched)) {
                              DF_Output_Patient_NonHotspot_int <- rbind.fill(DF_Output_Patient_NonHotspot_int, 
                                                                             cbind(patient_INFO,trial_INFO_matched[row_No,]))
                            }
                            remove(trial_INFO_matched)
                          }
                          remove(patient_INFO,row_No)
                        }
                      }
                      remove(arm_id,arm_num)
                    }  
                    remove(arm.match)
                    if (isTRUE(exists("var.patient.list"))) {remove(var.patient.list)}
                    if (isTRUE(exists("var.exclude.list"))) {remove(var.exclude.list)}
                    if (isTRUE(exists("var.patient"))) {remove(var.patient)}
                    if (isTRUE(exists("var.exclude"))) {remove(var.exclude)}
                  }
                  remove(age_gate)
                }
                remove(pathogenic_id,pathogenicity_gate)
                
              }  else {
                DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                     DF_patient$var.type == var.type_id &
                                                                     DF_patient$Exon_Number == exon_id)] <-
                  "Variant Exon criteria NOT satisfied"
              }
              remove(exon_id,exon_num,exon.patient,exon.trial)
            }
            
          }  else {
            DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                 DF_patient$var.type == var.type_id)] <-
              "Variant Type criteria NOT satisfied"
          }
        }
        remove(var.type.trial,var.type.patient,var.type_id_syn,var.type_id,
               type_num,DF_Gene_Patient_NonHotspot)
        
      } else {
        DF_patient$Patient_Variant_NonHotspot_Status[which(DF_patient$VariantGene == gene_id)] <- 
          "Gene NOT found"
      }
    }
    
    DF_Output_Patient_NonHotspot <- rbind.fill(DF_Output_Patient_NonHotspot, DF_Output_Patient_NonHotspot_int) 
    remove(DF_Output_Patient_NonHotspot_int)
    
    # Write to match results per patient to local computer
    #----------------------------------------------
    write.table(DF_patient, file = paste(tempdir, patient_id, ".tsv", sep=""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  ## Remove rows that are all empty
  DF_Output_Patient_NonHotspot <- 
    DF_Output_Patient_NonHotspot[rowSums(is.na(DF_Output_Patient_NonHotspot)) != ncol(DF_Output_Patient_NonHotspot),]
  
  ## Write to local computer
  #----------------------------------------------
  write.table(DF_Output_Patient_NonHotspot, 
              file = paste(tempdir,"NCIMatch_NonHotspot_Matched_", Patient_Variant_Report_timestamp, "_", 
                           dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  remove(DF_patient,gene_id,gene_num,gene.patient,genes.Inclusion_NonHotspot,
         ncol_InclusionNonHotspot,ncol_STAMP,patient_id,patient_num)
  
  if (isTRUE(exists("pt_rowNo"))){remove(pt_rowNo)}
  if (isTRUE(exists("DF_Exclude_Arm"))){remove(DF_Exclude_Arm)}
  if (isTRUE(exists("exclusion_comment"))){remove(exclusion_comment)}
}
