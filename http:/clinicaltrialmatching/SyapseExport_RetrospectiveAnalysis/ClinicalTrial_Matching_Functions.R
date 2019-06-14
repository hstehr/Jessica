#################################
## Customized color palattes
#################################
## http://tools.medialab.sciences-po.fr/iwanthue/
## Parameters: n=XX soft (k-means), colorblind-friendly, hard (force vector repulsion algorithm)
## http://tools.medialab.sciences-po.fr/iwanthue/theory.php
custom.hues.60 = c("#135300","#ff8af4","#51b949","#621488","#d1dd59","#00238a","#b6b520","#3a6ae1","#5d990d","#6f7ffa",
                   "#e0991e","#0060c6","#75ec95","#81007f","#00741f","#c155c4","#008b48","#de51b5","#009f69","#d0237c",
                   "#01e1e8","#c90b55","#01c4bc","#ad0032","#019368","#ff65bc","#094600","#a28eff","#9d8400","#1996ff",
                   "#b86200","#0159b2","#e5d56f","#30165a","#f4d07c","#00347a","#ff7b4c","#0191d2","#863800","#6cb3ff",
                   "#626f00","#be9aff","#8a984f","#930070","#7ea6ff","#8c0039","#0161a5","#ff6f96","#b6abff","#ac005f",
                   "#968cd3","#ff64aa","#654a8a","#ff97b5","#580054","#ffa7f9","#601041","#ff85cd","#78004d","#af535e")

custom.hues.30 = c("#47b041","#4e40b2","#9be876","#830074","#529100","#867cf7","#b7a606","#013c9f","#ded86d","#007ee5",
                   "#ed852d","#006baf","#af3109","#5ceacc","#dd3283","#00b973","#f95fbb","#004d0b","#fe93ff","#738700",
                   "#ce9fff","#685f00","#ff83da","#a96500","#ae75b2","#e9d387","#8a0038","#ff798b","#612500","#c71e3b")

custom.hues.10 = c("#d52373","#41e185","#310b69","#5b9503","#c460d1","#004002","#ff8179","#5d003f","#f07138","#7d2200")

custom.hues.9 = c("#90001c","#2b8b13","#692697","#f9cf6c","#01478c","#b34200","#c42d93","#ddae6e","#ff85a9")

custom.hues.8 = c("#9d006c","#006938","#8e007d","#fbcd7e","#6b1f4e","#cb511c","#ff8fd5","#ff5d74")

custom.hues.7 = c("#ca372d","#0294ea","#5abb48","#9350c4","#b1e28f","#ff88f7","#ffbf63")

custom.hues.6 = c("#ff7eaf","#6d8b00","#470064","#f5d248","#a30028","#ffa16a")

custom.hues.5 = c("#2d1783","#a84c00","#a092ff","#b20049","#00326f")

custom.hues.4 = c("#bb7438","#7f64b9","#72ac5c","#b94b75")

custom.hues.3 = c("#f3a632","#60e9d9","#ba005b")

custom.hues.2 = c("#a2001e","#9acd46")

#################################
## FUNCTIONS
#################################
shorthand_visual_fxn <- function (DF) {
  # Collapse similar primary tumor site
  DF$PrimaryTumorSite[which(DF$PrimaryTumorSite %in% c("colon","colon and rectum"))] <- "colon and rectum"
  DF$PrimaryTumorSite[which(DF$PrimaryTumorSite %in% c("liver","hepatocellular (liver)"))] <- "liver"
  DF$PrimaryTumorSite[which(DF$PrimaryTumorSite %in% c("testes","testis"))] <- "testes"
  
  # Abbreviate for visualization
  DF$PrimaryTumorSite[which(DF$PrimaryTumorSite == "central nervous system (brain/spinal cord)")] <- "cns (brain/spinal cord)"
  DF$PrimaryTumorSite[which(DF$PrimaryTumorSite == "hematologic and lymphatic neoplasm")] <- "hematologic and lymphoid"
  # sort(unique(DF$PrimaryTumorSite))
  
  assign("DF", DF, envir = .GlobalEnv)
}

parameter_filter_fxn <- function() {
  ## Filters APPLIED
  if (isTRUE(adult.group_FILTER)) { ageName <- "AdultGroupON" } else { ageName <- "AdultGroupOFF" }
  if (isTRUE(pathogenic_FILTER)) { pathoName <- "pathogenicON" } else { pathoName <- "pathogenicOFF" }
  if (isTRUE(disease.group_FILTER)) { groupName <- "diseaseFILTER_groupON" } else { groupName <- "diseaseFILTER_groupOFF" }
  if (isTRUE(disease.group_FILTER & disease.site_FILTER)) { siteName <- "siteON" } else { siteName <- "siteOFF" }
  if (isTRUE(disease.code_FILTER)) { dxName <- "histologicaldxON" } else { dxName <- "histologicaldxOFF" }
  filterName <- paste(groupName, siteName, "_", dxName, "_", pathoName, "_",  ageName, sep="")
  
  assign("groupName", groupName, envir = .GlobalEnv)
  assign("siteName", siteName, envir = .GlobalEnv)
  assign("dxName", dxName, envir = .GlobalEnv)
  assign("pathoName", pathoName, envir = .GlobalEnv)
  assign("ageName", ageName, envir = .GlobalEnv)
  assign("filterName", filterName, envir = .GlobalEnv)
  
  ## Print parameters to output
  if (isTRUE(Internal_match & NCI_match)) {
    cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
        "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
        "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
        "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
        "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
        "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
        "Outdirectory: ", data.root, "\n",
        "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
        "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
    
  } else if (isTRUE(Internal_match)) {
    cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
        "Stanford Internal Trial Timestamp: ", OnCore_Biomarker_Report_timestamp, "\n",
        "\t", "FILTERs: disease group matched: ", disease.group_FILTER, "; disease site matched: ", disease.site_FILTER, "\n",
        "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
        "Outdirectory: ", data.root, "\n",
        "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
        "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
    
  } else if (isTRUE(NCI_match)) {
    cat("Syapse Timestamp: ", Syapse_Export_timestamp, "\n", "\n",
        "NCI-MATCH Trial Timestamp: ", Patient_Variant_Report_timestamp, "\n",
        "\t", "FILTERs: disease code matched: ", disease.code_FILTER, "\n",
        "FILTERs: adult patients (age >= 18yo): ",adult.group_FILTER, "; pathogenic variants: ", pathogenic_FILTER, "\n", "\n",
        "Outdirectory: ", data.root, "\n",
        "Temporary Files within outdir: ", gsub(data.root, "", tempdir), "\n",
        "Patient Match Results within outdir: ", gsub(data.root, "", outdir), "\n","\n")
    
  } else {
    sink(file = err.output, append = TRUE, split = FALSE)
    options(max.print=999999)
    
    cat("NO TRIAL has been specified to be matched. Check parameter specification.", "\n","\n")
    
    sink()
    
    sink(file = out.output, append = TRUE, split = FALSE)
    options(max.print=999999)
    
  }
}

Iterate_Fxn <- function(STAMP_DF_iterate,STAMP_Fusion_iterate,STAMP_CNV_iterate) {
  
  SNVIndel.list <- sort(unique(STAMP_DF_iterate$PatientID))
  Fusion.list <- sort(unique(STAMP_Fusion_iterate$PatientID))
  CNV.list <- sort(unique(STAMP_CNV_iterate$PatientID))
  
  patient.list.full <- sort(unique(append(SNVIndel.list,append(Fusion.list,CNV.list))))
  assign("patient.list.full", patient.list.full, envir = .GlobalEnv)
  
  if (length(patient.list.full) > 0) {
    
    if (isTRUE(Internal_match)) {
      source(paste(script.root,"Biomarker_Report_QC.R",sep=""))
      
      patient.list=SNVIndel.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Biomarker_Report_Match_SNVIndel.R",sep="")) 
      
      patient.list=CNV.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Biomarker_Report_Match_CNV.R",sep=""))
      
      patient.list=Fusion.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Biomarker_Report_Match_Fusion.R",sep=""))
    }
    
    if (isTRUE(NCI_match)) {
      source(paste(script.root,"Patient_Variant_Report_QC.R",sep=""))
      
      patient.list=SNVIndel.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Patient_Variant_Report_InclusionMatch_SNVIndel.R",sep=""))
      source(paste(script.root,"Patient_Variant_Report_NonHotspotMatch_SNVIndel.R",sep=""))
      
      patient.list=CNV.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Patient_Variant_Report_InclusionMatch_CNV.R",sep=""))
      source(paste(script.root,"Patient_Variant_Report_NonHotspotMatch_CNV.R",sep=""))
      
      patient.list=Fusion.list
      assign("patient.list", patient.list, envir = .GlobalEnv)
      source(paste(script.root,"Patient_Variant_Report_InclusionMatch_Fusion.R",sep=""))
    }
    
    patient.list=patient.list.full
    # Generate OUTPUT files
    source(paste(script.root,"ClinicalTrial_Output_Details.R",sep=""))
    source(paste(script.root,"ClinicalTrial_Output_tsv.R",sep=""))
    source(paste(script.root,"ClinicalTrial_Output_Candidates.R",sep=""))
  }
}

plot_Trial_Count <- function(DF, trialColumn, fileName_pre, plotTitle, ymax, outdir) {
  
  # Y-axis parameters
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # X-axis parameters
  x_max = max(as.numeric(DF[[trialColumn]]))
  if (x_max > 2) { image_width = x_max +2
  } else { image_width = x_max +3
  }
  
  # Specify palette 
  custom.palette = get(paste("custom.hues.",nrow(DF),sep=""))
  
  # Specify match number
  total_matched = sum(DF$No.Orders[which(DF[[trialColumn]] != 0)])
  
  plot <- ggplot(DF, aes(x=get(trialColumn), y=No.Orders, fill=as.factor(get(trialColumn)))) +
    geom_bar(stat="identity") +
    geom_text(aes(label=No.Orders), vjust=-0.75, size=4.0) +
    
    labs(title = plotTitle,
         subtitle = paste("Match rate: ",
                          format(round(as.numeric(100*total_matched/sum(DF$No.Orders)),1),nsmall = 1),
                          "% (n = ", total_matched, "/", sum(DF$No.Orders), ")", sep="")) +
    
    xlab("No. Clinical Trials Matched") +
    scale_y_continuous(name="No. Test Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) +
    
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
          plot.subtitle = element_text(hjust=1, face="bold",size=12),
          legend.text=element_text(size=10),
          axis.text=element_text(size=12,face="bold"), 
          axis.title=element_text(size=12,face="bold"),
          
          legend.position="none") +
    
    scale_fill_manual(values = custom.palette)
  
  tiff(filename = paste(outdir, fileName_pre, "_", Sys.Date(), ".tiff", sep=""),
       width = image_width, height = 7, units = "in", res = 350)
  print(plot)
  dev.off()
}

OnCore_Tally_fxn <- function() {
  
  colnames_extract <- c("PatientID","OnCore.No","Match_Detail","Gene")
  colnames_extract.snv <- c("PatientID","OnCore.No","Match_Detail","VariantGene")
  colnames_extract.cnv <- c("PatientID","OnCore.No","Match_Detail","CNV_Gene")

  # Import candidate matches and specify empty dataframes
  #---------------------------------------------- 
  # SNV Indels
  file_name = paste(tempdir,"OnCore_SNVIndel_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_OnCore) > 0) {continue_SNVIndel <- as.logical("TRUE")
    } else {continue_SNVIndel <- as.logical("FALSE")}
    
  } else {continue_SNVIndel <- as.logical("FALSE")}
  
  # CNVs
  file_name = paste(tempdir,"OnCore_CNV_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_OnCore) > 0) {continue_CNV <- as.logical("TRUE")
    } else {continue_CNV <- as.logical("FALSE")}
    
  } else {continue_CNV <- as.logical("FALSE")}
  
  # Fusions
  file_name = paste(tempdir,"OnCore_Fusion_Matched_", OnCore_Biomarker_Report_timestamp, "_", 
                    groupName,siteName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_Fusion_OnCore <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_Fusion_OnCore) > 0) {continue_Fusion <- as.logical("TRUE")
    } else {continue_Fusion <- as.logical("FALSE")}
    
  } else {continue_Fusion <- as.logical("FALSE")}
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_OnCore_FINAL <- data.frame(matrix(ncol = (length(colnames_extract) +1), nrow = 0))
  colnames(Output_OnCore_FINAL) <- c(colnames_extract,"Gene")
  
  if (isTRUE(continue_SNVIndel)) {
    Output_SNVIndel_OnCore$Match_Detail <- "SNVIndel_OnCore"
    
    Output_SNVIndel_OnCore <- Output_SNVIndel_OnCore[,colnames_extract.snv]
    colnames(Output_SNVIndel_OnCore) <- colnames_extract
    
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_SNVIndel_OnCore)
  }
  
  if (isTRUE(continue_CNV)) {
    Output_CNV_OnCore$Match_Detail <- "CNV_OnCore"
    
    Output_CNV_OnCore <- Output_CNV_OnCore[,colnames_extract.cnv]
    colnames(Output_CNV_OnCore) <- colnames_extract
    
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_CNV_OnCore)
  }
  
  if (isTRUE(continue_Fusion)) {
    Output_Fusion_OnCore$Match_Detail <- "Fusion_OnCore"
    Output_OnCore_FINAL <- rbind(Output_OnCore_FINAL, Output_Fusion_OnCore[,colnames_extract])
  }
  
  Output_OnCore_FINAL <- unique(Output_OnCore_FINAL[,])

  # Breakdown by trial & gene
  OnCoreMatch_tally <- data.frame(Output_OnCore_FINAL %>% group_by(OnCore.No,Gene) %>% tally())
  colnames(OnCoreMatch_tally) <- c("OnCore.No","Gene","No.Orders")
  OnCoreMatch_tally <- OnCoreMatch_tally[order(OnCoreMatch_tally$No.Orders, decreasing = TRUE),]
  
  assign("OnCore_SNVIndel_Matched_TrialGene", OnCoreMatch_tally, envir = .GlobalEnv)
  
  write.table(OnCoreMatch_tally, file = paste(outdir, "OnCore_SNVIndel_Matched_TrialGene.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)
  
  # Breakdown by trial 
  OnCoreMatch_tally <- data.frame(Output_OnCore_FINAL %>% group_by(OnCore.No) %>% tally())
  colnames(OnCoreMatch_tally) <- c("OnCore.No","No.Orders")
  OnCoreMatch_tally <- OnCoreMatch_tally[order(OnCoreMatch_tally$No.Orders, decreasing = TRUE),]
  
  write.table(OnCoreMatch_tally, file = paste(outdir, "OnCore_SNVIndel_Matched_Trial.tsv", sep=""),
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  assign("OnCore_SNVIndel_Matched_Trial", OnCoreMatch_tally, envir = .GlobalEnv)
  
  remove(Output_OnCore_FINAL,OnCoreMatch_tally)
}

LabRep_Summary_fxn <- function() {
  
  colnames_extract <- c("PatientID","Arm_Name","Variant_Type","VariantGene","HistologicalDx","PrimaryTumorSite",
                        "VariantPathogenicityStatus","Variant_Detail",
                        "VariantHGVSProtein","VariantHGVSCoding","VariantHGVSGenomic")
  
  colnames_extract.snv <- c("PatientID","Arm_Name","var.type","VariantGene","HistologicalDx","PrimaryTumorSite",
                            "VariantPathogenicityStatus","Variant_Detail",
                            "VariantHGVSProtein","VariantHGVSCoding","VariantHGVSGenomic")
  
  colnames_extract.cnv <- c("PatientID","Arm_Name","var.type","CNV_Gene","HistologicalDx","PrimaryTumorSite",
                            "VariantPathogenicityStatus","var.anno",
                            "VariantHGVSProtein","VariantHGVSCoding","VariantHGVSGenomic")
  
  colnames_extract.fusion <- c("PatientID","Arm_Name","var.type","Gene_Name","HistologicalDx","PrimaryTumorSite",
                               "VariantPathogenicityStatus","Fusion_Detail",
                               "VariantHGVSProtein","VariantHGVSCoding","VariantHGVSGenomic")
  
  colnames_keep <- c("PatientID","Arm_Name","VariantLabel","Variant_Type",
                     "HistologicalDx","PrimaryTumorSite","VariantPathogenicityStatus")
  
  # Import candidate matches and specify empty dataframes
  #---------------------------------------------- 
  # SNVIndel_Variant
  file_name = paste(tempdir,"NCI_SNVIndel_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_Variant) > 0) {
      continue_SNVIndel_Variant <- as.logical("TRUE")
      
      Output_SNVIndel_Variant$Variant_Detail <- NA
      Output_SNVIndel_Variant <- unique(Output_SNVIndel_Variant[,colnames_extract.snv])
      colnames(Output_SNVIndel_Variant) <- colnames_extract
      
    } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
    
  } else {continue_SNVIndel_Variant <- as.logical("FALSE")}
  
  # CNV_Variant
  file_name = paste(tempdir,"NCIMatch_CNV_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_Variant) > 0) {
      continue_CNV_Variant <- as.logical("TRUE")
      
      Output_CNV_Variant$VariantHGVSProtein <- NA
      Output_CNV_Variant$VariantHGVSCoding <- NA
      Output_CNV_Variant$VariantHGVSGenomic <- NA
      Output_CNV_Variant <- unique(Output_CNV_Variant[,colnames_extract.cnv])
      colnames(Output_CNV_Variant) <- colnames_extract
      
    } else {continue_CNV_Variant <- as.logical("FALSE")}
    
  } else {continue_CNV_Variant <- as.logical("FALSE")}
  
  # Fusion_Variant
  file_name = paste(tempdir,"NCI_Fusion_Variant_Matched_", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_Fusion_Variant <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_Fusion_Variant) > 0) {
      continue_Fusion_Variant <- as.logical("TRUE")
      
      Output_Fusion_Variant$VariantHGVSProtein <- NA
      Output_Fusion_Variant$VariantHGVSCoding <- NA
      Output_Fusion_Variant$VariantHGVSGenomic <- NA
      Output_Fusion_Variant <- unique(Output_Fusion_Variant[,colnames_extract.fusion])
      colnames(Output_Fusion_Variant) <- colnames_extract
      
    } else {continue_Fusion_Variant <- as.logical("FALSE")}
    
  } else {continue_Fusion_Variant <- as.logical("FALSE")}
  
  # SNVIndel_NonHotspot
  file_name = paste(tempdir,"NCI_SNVIndel_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_SNVIndel_NonHotspot <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_SNVIndel_NonHotspot) > 0) {
      continue_SNVIndel_NonHotspot <- as.logical("TRUE")
      
      Output_SNVIndel_NonHotspot$Variant_Detail <- NA
      Output_SNVIndel_NonHotspot <- unique(Output_SNVIndel_NonHotspot[,colnames_extract.snv])
      colnames(Output_SNVIndel_NonHotspot) <- colnames_extract
      
    } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
    
  } else {continue_SNVIndel_NonHotspot <- as.logical("FALSE")}
  
  # CNV_NonHotspot
  file_name = paste(tempdir,"NCI_CNV_NonHotspot_Matched", Patient_Variant_Report_timestamp, "_", 
                    dxName, "_", pathoName, "_",  ageName, ".tsv", sep="")
  
  if (file.exists(file_name)) {
    Output_CNV_NonHotspot <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    if (nrow(Output_CNV_NonHotspot) > 0) {
      continue_CNV_NonHotspot <- as.logical("TRUE")
      
      Output_CNV_NonHotspot$VariantHGVSProtein <- NA
      Output_CNV_NonHotspot$VariantHGVSCoding <- NA
      Output_CNV_NonHotspot$VariantHGVSGenomic <- NA
      Output_CNV_NonHotspot <- unique(Output_CNV_NonHotspot[,colnames_extract.cnv])
      colnames(Output_CNV_NonHotspot) <- colnames_extract
      
    } else {continue_CNV_NonHotspot <- as.logical("FALSE")}
    
  } else {continue_CNV_NonHotspot <- as.logical("FALSE")}
  
  # Merge into single dataframe
  #---------------------------------------------- 
  Output_NCI_FINAL <- data.frame(matrix(ncol = length(colnames_extract), nrow = 0))
  colnames(Output_NCI_FINAL) <- colnames_extract
  
  if (isTRUE(continue_SNVIndel_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_Variant)
  }
  
  if (isTRUE(continue_CNV_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_Variant)
  }
  
  if (isTRUE(continue_Fusion_Variant)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_Fusion_Variant)
  }
  
  if (isTRUE(continue_SNVIndel_NonHotspot)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_SNVIndel_NonHotspot)
  }
  
  if (isTRUE(continue_CNV_NonHotspot)) {
    Output_NCI_FINAL <- rbind(Output_NCI_FINAL, Output_CNV_NonHotspot)
  }
  
  # Remove duplicates
  Output_NCI_FINAL <- unique(Output_NCI_FINAL[order(Output_NCI_FINAL$Arm_Name, decreasing = FALSE),])
  
  Output_NCI_FINAL$Variant_Detail <- gsub("amplification","AMP",Output_NCI_FINAL$Variant_Detail)
  Output_NCI_FINAL$Variant_Detail <- gsub("deletion","DEL",Output_NCI_FINAL$Variant_Detail)
  Output_NCI_FINAL$VariantPathogenicityStatus <- gsub("NULL",NA,Output_NCI_FINAL$VariantPathogenicityStatus)
  
  Output_NCI_FINAL$VariantLabel <- NA
  for (row_No in 1:nrow(Output_NCI_FINAL)) {
    
    if (Output_NCI_FINAL$Variant_Type[row_No] == "CNV") {
      Output_NCI_FINAL$VariantLabel[row_No] <- paste(Output_NCI_FINAL$VariantGene[row_No], 
                                                     " ",Output_NCI_FINAL$Variant_Detail[row_No],sep="")
      
    } else if (Output_NCI_FINAL$Variant_Type[row_No] == "Fusion") {
      Output_NCI_FINAL$VariantLabel[row_No] <- paste(Output_NCI_FINAL$VariantGene[row_No], 
                                                     " (",Output_NCI_FINAL$Variant_Detail[row_No]," Fusion)",sep="")
      
    } else {
      Output_NCI_FINAL$VariantLabel[row_No] <- paste(Output_NCI_FINAL$VariantGene[row_No]," ",
                                                     Output_NCI_FINAL$VariantHGVSCoding[row_No]," (",
                                                     Output_NCI_FINAL$VariantHGVSProtein[row_No],")",sep="")
    }
  }
  
  Output_NCI_FINAL <- unique(Output_NCI_FINAL[,colnames_keep])
  
  # Import patient DOB
  Output_NCI_FINAL <- left_join(Output_NCI_FINAL,TRF_DF[,c("PatientID","PatientDOB")],
                                by = "PatientID")
  
  assign("Output_NCI_FINAL", Output_NCI_FINAL, envir = .GlobalEnv)
  
  write.table(Output_NCI_FINAL, file = paste(outdir,"Output_Patient_List.tsv",sep=""), 
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)
}
