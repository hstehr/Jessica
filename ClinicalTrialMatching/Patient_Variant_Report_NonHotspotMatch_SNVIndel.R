# Match by gene > variant function/oncominevariantclass > Exon No \
# > pathogenicity.status > age.group (i.e. adult) > exclusion variant > exclusion nonhotspot
## Search Output: individual patient file with stop point
## Match Output: "NCI_SNVIndel_NonHotspot_Matched.tsv"

# Code SNVIndel dataframe: lines 11-15

if (isTRUE(NCI_match)) {
  
  ncol_InclusionNonHotspot <- as.numeric(ncol(Inclusion_NonHotspot_Rules))
  ncol_SNVIndel = 13
  SNVIndel.colnames <- c("PatientID","VariantHGVSGenomic","VariantLabel","VariantGene","VariantHGVSCoding",
                         "VariantHGVSProtein","var.type","var.anno","Exon_Number",
                         "PrimaryTumorSite.Category","PrimaryTumorSite","VariantPathogenicityStatus",
                         "HistologicalDx")
  
  # Extract gene from NCI-MATCH
  genes.Inclusion_NonHotspot <- sort(unique(Inclusion_NonHotspot_Rules$Gene_Name))
  
  # Output file of matched mutations with clinical trials
  DF_Output_Patient_NonHotspot <- data.frame(matrix(NA, ncol = (ncol_SNVIndel + ncol_InclusionNonHotspot)))
  colnames(DF_Output_Patient_NonHotspot) <- append(SNVIndel.colnames, colnames(Inclusion_NonHotspot_Rules))
  
  ## Iterate through each patient_id of patient.list
  #----------------------------------------------
  for (patient_num in 1:length(patient.list)) {
    patient_id <- patient.list[patient_num]
    
    DF_Output_Patient_NonHotspot_int <- data.frame(matrix(NA, ncol = (ncol_SNVIndel + ncol_InclusionNonHotspot)))
    colnames(DF_Output_Patient_NonHotspot_int) <- append(SNVIndel.colnames, colnames(Inclusion_NonHotspot_Rules))
    
    # Import STAMP entries per patient
    DF_patient <- read.csv(file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                           header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(DF_patient) > 0) {
      # Extract unique genes per patient
      gene.patient <- sort(unique(DF_patient$VariantGene))
      
      # Clinical match result
      DF_patient$NCI_SNVIndel_NonHotspot_Status <- NA
      
      ## Iterate through each gene_id of gene.patient
      #----------------------------------------------
      for (gene_num in 1:length(gene.patient)) {
        gene_id <- gene.patient[gene_num]
        
        ## Match Biomarker_GeneName
        #----------------------------------------------
        if (gene_id %in% genes.Inclusion_NonHotspot) {
          
          # Extract gene_id match from Inclusion_NonHotspot_Rules
          DF_Gene_Patient_NonHotspot <- 
            data.frame(Inclusion_NonHotspot_Rules[which(Inclusion_NonHotspot_Rules$Gene_Name == gene_id),])
          
          ## Parse Variant_Type
          #----------------------------------------------
          var.type.trial <- unique(DF_Gene_Patient_NonHotspot$Function[!is.na(DF_Gene_Patient_NonHotspot$Function)])
          var.group.trial <- unique(tolower(DF_Gene_Patient_NonHotspot$oncominevariantclass[!is.na(DF_Gene_Patient_NonHotspot$oncominevariantclass)]))
          
          var.type.patient <- unique(DF_patient$var.type[DF_patient$VariantGene == gene_id])
          
          ## Iterate through each bio.cond_id of var.type.patient
          #----------------------------------------------
          for (type_num in 1:length(var.type.patient)) {
            var.type_id <- var.type.patient[type_num]
            
            # Translate terminology to common syntax
            if (isTRUE(!is.na(var.type_id))) {
              var.type_id_syn.patient <- c("deleterious")  # Any mutation may be deleterious >> cast wide net
              
              if (isTRUE("Deletion" %in% var.type_id)) {
                var.type_id_syn.patient <- append(var.type_id_syn.patient, "nonframeshiftDeletion")
              }
              
              if (isTRUE("Insertion" %in% var.type_id)) {
                var.type_id_syn.patient <- append(var.type_id_syn.patient, "nonframeshiftInsertion")
              }
              
              if (isTRUE("Delins" %in% var.type_id)) {
                var.type_id_syn.patient <- append(var.type_id_syn.patient, c("nonframeshiftDeletion","nonframeshiftInsertion"))
              }
              
              var.type_id_syn.patient <- unique(var.type_id_syn.patient)
            }
            
            var.type_id_syn.patient <- 
              var.type_id_syn.patient[which(var.type_id_syn.patient %in% append(var.type.trial, var.group.trial))]
            
            if (isTRUE(length(var.type_id_syn.patient) > 0)) {next_continue <- as.logical("TRUE")
            } else {next_continue <- as.logical("FALSE")}
            
            ## Match Variant_Type
            #----------------------------------------------
            if (isTRUE(next_continue)) {
              for (type_syn_num in 1:length(var.type_id_syn.patient)) {
                var.type_id_syn <- var.type_id_syn.patient[type_syn_num]
                
                ## Parse Exon Number
                #----------------------------------------------
                exon.trial <- 
                  unique(DF_Gene_Patient_NonHotspot$Exon[which(DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                                 (DF_Gene_Patient_NonHotspot$Function == var.type_id_syn |
                                                                    tolower(DF_Gene_Patient_NonHotspot$oncominevariantclass) == var.type_id_syn))])
                
                exon.patient <- unique(DF_patient$Exon_Number[
                  DF_patient$VariantGene == gene_id & DF_patient$var.type == var.type_id])
                
                # Remove no match Exon Number from DF_Gene_Patient_NonHotspot
                if (isTRUE(length(exon.trial) > 0)) {
                  
                  for (exon_num in 1:length(exon.trial)) {
                    exon_id <- exon.trial[exon_num]
                    
                    if (isTRUE(!is.na(exon_id))) {
                      if (!(exon_id %in% exon.patient)) {
                        
                        DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[which(
                          (DF_Gene_Patient_NonHotspot$Exon != exon_id | is.na(DF_Gene_Patient_NonHotspot$Exon)) &
                            (DF_Gene_Patient_NonHotspot$Function == var.type_id_syn | 
                               tolower(DF_Gene_Patient_NonHotspot$oncominevariantclass) == var.type_id_syn)),]
                      }
                    }
                  }
                }
                
                ## Iterate through each exon_id of exon.patient
                #----------------------------------------------
                for (exon_num in 1:length(exon.patient)) {
                  exon_id <- exon.patient[exon_num]
                  
                  ## Match Exon Number
                  #----------------------------------------------
                  if (isTRUE((exon_id %in% exon.trial | isTRUE("TRUE" %in% is.na(exon.trial))) & 
                             nrow(DF_Gene_Patient_NonHotspot) > 0)) {
                    
                    ## Iterate through each genomic_id of genomic.list
                    #----------------------------------------------
                    genomic.list <- sort(unique(DF_patient$VariantHGVSGenomic[which(DF_patient$VariantGene == gene_id &
                                                                                      DF_patient$var.type == var.type_id &
                                                                                      DF_patient$Exon_Number == exon_id)]))
                    for (genomic_num in 1:length(genomic.list)) {
                      genomic_id <- genomic.list[genomic_num]
                      
                      ## Assess Pathogenicity Status
                      #----------------------------------------------
                      pathogenicity_gate <- NA
                      
                      if (isTRUE(pathogenic_FILTER)) {
                        if (isTRUE(DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                               DF_patient$var.type == var.type_id &
                                                                               DF_patient$Exon_Number == exon_id & 
                                                                               DF_patient$VariantHGVSGenomic == genomic_id)] 
                                   %in% pathogenic_accepted)) {
                          pathogenicity_gate <- as.logical("TRUE")
                          
                        } else {
                          pathogenicity_gate <- as.logical("FALSE")
                          DF_patient$NCI_SNVIndel_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                            DF_patient$var.type == var.type_id &
                                                                            DF_patient$Exon_Number == exon_id & 
                                                                            DF_patient$VariantHGVSGenomic == genomic_id)] <-
                            "Pathogenicity criteria NOT satisfied"
                        }
                      }
                      
                      pathogenic_id <- unique(DF_patient$VariantPathogenicityStatus[which(DF_patient$VariantGene == gene_id &
                                                                                            DF_patient$var.type == var.type_id &
                                                                                            DF_patient$Exon_Number == exon_id & 
                                                                                            DF_patient$VariantHGVSGenomic == genomic_id)])
                      
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
                            DF_patient$NCI_SNVIndel_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                              DF_patient$var.type == var.type_id &
                                                                              DF_patient$Exon_Number == exon_id & 
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
                                              DF_patient$Exon_Number == exon_id & 
                                              DF_patient$VariantHGVSGenomic == genomic_id &
                                              DF_patient$VariantPathogenicityStatus == pathogenic_id)
                          
                          arm.match <- unique(append(
                            DF_Gene_Patient_NonHotspot$Arm_Name[which(DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                                        DF_Gene_Patient_NonHotspot$Function == var.type_id_syn &
                                                                        DF_Gene_Patient_NonHotspot$Exon == exon_id)],
                            DF_Gene_Patient_NonHotspot$Arm_Name[which(DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                                        tolower(DF_Gene_Patient_NonHotspot$oncominevariantclass) == var.type_id_syn)]
                          ))
                          
                          DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo] <- paste("Candidate trial IDENTIFIED: ", 
                                                                                       paste(arm.match, collapse = ", "), sep="")
                          
                          ## Assess Exclusions
                          #----------------------------------------------
                          if (length(arm.match) > 0) {
                            
                            for (arm_num in 1:length(arm.match)) {
                              
                              arm_id <- arm.match[arm_num]
                              exclusion_continue <- as.logical("TRUE")
                              
                              # Import CNV entries per patient
                              #----------------------------------------------
                              DF_patient_CNV <- read.csv(file = paste(tempdir, patient_id, "_CNV.tsv", sep=""),
                                                         header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                              
                              if (nrow(DF_patient_CNV) > 0) {exclusion_CNV_continue <- as.logical("TRUE")  
                              } else {exclusion_CNV_continue <- as.logical("FALSE")}
                              
                              # Import Fusion entries per patient
                              #----------------------------------------------
                              DF_patient_Fusion <- read.csv(file = paste(tempdir, patient_id, "_Fusion.tsv", sep=""),
                                                            header = TRUE, na.strings = c("NA"), stringsAsFactors = FALSE, sep = "\t")
                              
                              if (nrow(DF_patient_Fusion) > 0) {exclusion_Fusion_continue <- as.logical("TRUE")  
                              } else {exclusion_Fusion_continue <- as.logical("FALSE")}
                              
                              ## Assess Exclusion Variants
                              #----------------------------------------------
                              # Extract variants for Arm_No in Exclusion_Variants
                              DF_Exclude_Arm <- Exclusion_Variants[Exclusion_Variants$Arm_Name == arm_id,]
                              
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
                                  DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                                  
                                  cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                  
                                  # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                                  DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                                  
                                  DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                                  
                                  # Checkpoint
                                  exclusion_continue <- as.logical("FALSE")
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
                                  DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                                  
                                  cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                  
                                  # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                                  DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                                  
                                  DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                                  
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
                                  DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo] <- exclusion_comment
                                  
                                  cat(paste(patient_id, ": ", exclusion_comment, sep=""),"\n","\n")
                                  
                                  # Remove excluded ARM_No from DF_Gene_Patient_NonHotspot
                                  DF_Gene_Patient_NonHotspot <- DF_Gene_Patient_NonHotspot[DF_Gene_Patient_NonHotspot$Arm_Name != arm_id,]
                                  
                                  DF_Output_Patient_NonHotspot_int <- DF_Output_Patient_NonHotspot_int[DF_Output_Patient_NonHotspot_int$Arm_Name != arm_id,]
                                  
                                  # Checkpoint
                                  exclusion_continue <- as.logical("FALSE")
                                }
                              }
                              remove(DF_Exclude_Arm_SNV, DF_Exclude_Arm_CNV, DF_Exclude_Arm, DF_Exclude_Arm_Fusion)
                              
                              ## Assess Disease Exclusion Codes
                              #----------------------------------------------
                              # dx.patient == "NULL" indicates that histological dx is not available 
                              dx.patient = unique(DF_patient$HistologicalDx)
                              tumorsite.patient = unique(DF_patient$PrimaryTumorSite)
                              
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
                                  DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo] <- disease_comment
                                  
                                  cat(paste(patient_id, ": ", DF_patient$NCI_SNVIndel_NonHotspot_Status[pt_rowNo], sep=""))
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
                                                                   DF_patient$Exon_Number == exon_id & 
                                                                   DF_patient$VariantHGVSGenomic == genomic_id), c(1:ncol_SNVIndel)]
                                
                                if (nrow(patient_INFO) > 0) {
                                  
                                  # Matched trial INFO
                                  trial_INFO_matched <- 
                                    DF_Gene_Patient_NonHotspot[which(DF_Gene_Patient_NonHotspot$Arm_Name == arm_id &
                                                                       DF_Gene_Patient_NonHotspot$Gene_Name == gene_id &
                                                                       (DF_Gene_Patient_NonHotspot$Function == var.type_id_syn |
                                                                          tolower(DF_Gene_Patient_NonHotspot$oncominevariantclass) == var.type_id_syn) &
                                                                       (DF_Gene_Patient_NonHotspot$Exon == exon_id | is.na(DF_Gene_Patient_NonHotspot$Exon))),]
                                  
                                  # Append merged patient INFO with trial INFO to Output file
                                  DF_Output_Patient_NonHotspot_int <- rbind.fill(DF_Output_Patient_NonHotspot_int, 
                                                                                 crossing(patient_INFO, trial_INFO_matched))
                                  remove(trial_INFO_matched)
                                }
                                remove(patient_INFO)
                              }
                              remove(exclusion_CNV_continue,exclusion_continue,exclusion_Fusion_continue,
                                     DF_patient_CNV,DF_patient_Fusion)
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
                      remove(pathogenic_id,pathogenicity_gate,genomic_id,genomic_num)
                    }
                    remove(genomic.list)
                    
                  }  else {
                    DF_patient$NCI_SNVIndel_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                      DF_patient$var.type == var.type_id &
                                                                      DF_patient$Exon_Number == exon_id)] <-
                      "Variant Exon criteria NOT satisfied"
                  }
                  remove(exon_id,exon_num)
                }
                remove(var.type_id_syn,exon.trial,exon.patient)
              }
              remove(type_syn_num)
              
            } else {
              DF_patient$NCI_SNVIndel_NonHotspot_Status[which(DF_patient$VariantGene == gene_id &
                                                                DF_patient$var.type == var.type_id)] <-
                "Variant Type criteria NOT satisfied"
            }
            remove(next_continue)
          }
          remove(type_num,DF_Gene_Patient_NonHotspot,var.type.patient,var.group.trial,
                 var.type.trial,var.type_id)
          
        } else {
          DF_patient$NCI_SNVIndel_NonHotspot_Status[which(DF_patient$VariantGene == gene_id)] <- 
            "Gene NOT found"
        }
      }
      
      DF_Output_Patient_NonHotspot <- rbind.fill(DF_Output_Patient_NonHotspot, DF_Output_Patient_NonHotspot_int) 
      remove(DF_Output_Patient_NonHotspot_int)
      
      # Write to match results per patient to local computer
      #----------------------------------------------
      write.table(DF_patient, file = paste(tempdir, patient_id, "_SNVIndel.tsv", sep=""),
                  append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
                  quote = FALSE)
    }
  }
  
  ## Remove rows that are all empty
  DF_Output_Patient_NonHotspot <- 
    DF_Output_Patient_NonHotspot[rowSums(is.na(DF_Output_Patient_NonHotspot)) != ncol(DF_Output_Patient_NonHotspot),]
  
  ## Write to local computer
  #----------------------------------------------
  write.table(DF_Output_Patient_NonHotspot, 
              file = paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                           dxName, "_", pathoName, "_",  ageName, ".tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  remove(DF_patient,gene_id,gene_num,gene.patient,genes.Inclusion_NonHotspot,
         ncol_InclusionNonHotspot,ncol_SNVIndel,patient_id,patient_num,
         DF_Output_Patient_NonHotspot,SNVIndel.colnames)
  
  if (isTRUE(exists("pt_rowNo"))){remove(pt_rowNo)}
  if (isTRUE(exists("exclusion_comment"))){remove(exclusion_comment)}
  if (isTRUE(exists("var.type_id_syn.patient"))){remove(var.type_id_syn.patient)}
  
}
