# Generate single level pie graph
# Illustrate age and variant type distribution
# DF variable assigned as global variable in pipeline

cat(paste("Timestamp of data visualization generation START: ", Sys.time(), sep=""),"\n")

library("gridExtra")
library("gtable")
library("grid")
library("ggrepel")
library("RColorBrewer")

## STAMP database: age distribution by gender
#----------------------------------------------
if (isTRUE(static.plots_FILTER)) {
  # Subset relevant columns
  DF_subset <- data.frame(patient.id = DF$PatientID, age = DF$PatientAge,
                          gender = DF$PatientGender, stringsAsFactors = FALSE)
  # Remove duplicate entries 
  DF_subset <- DF_subset %>% dplyr::distinct(patient.id, .keep_all = TRUE)
  
  # Generate histogram
  #----------------------------------------------
  plot <- ggplot(DF_subset, aes(DF_subset$age, fill=DF_subset$gender)) +
    geom_histogram(bins = (max(DF_subset$age) - min(DF_subset$age) +1), col="gray") +
    
    labs(title = "STAMP Database: Age Distribution",
         subtitle = paste("N = ", nrow(DF_subset), sep="")) +
    
    scale_x_continuous(name="Current Age", breaks = seq(0, 100, 5)) +
    scale_y_continuous(name="Number of Samples", breaks = seq(0,150,10)) + 
    scale_fill_discrete(name = "Gender") +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=15),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.955, 0.88),
          legend.text=element_text(size=10),
          
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  # Generate output table
  #----------------------------------------------
  # Specify Age.Cohort
  DF_subset$Age.Cohort <- NA
  for (row_No in 1:nrow(DF_subset)) {
    if (isTRUE(DF_subset$age[row_No] < 18)) {
      DF_subset$Age.Cohort[row_No] <- "Child (< 18yo)"
    } else if (isTRUE(DF_subset$age[row_No] < 65)) {
      DF_subset$Age.Cohort[row_No] <- "Adult (18-64yo)"
    } else {
      DF_subset$Age.Cohort[row_No] <- "Older Adult (>= 65yo)"
    }
  }
  
  # Tabulate frequency
  DF_subset_tabulate <- data.frame(DF_subset %>% group_by(Age.Cohort,gender) %>% tally())
  colnames(DF_subset_tabulate) <- c("Age.Cohort","Gender","No.Samples")
  
  # Convert table to wide format
  DF_subset_tabulate <- data.frame(cast(DF_subset_tabulate, Age.Cohort ~ Gender),
                                   stringsAsFactors = FALSE)
  DF_subset_tabulate$Total <- DF_subset_tabulate$Female + DF_subset_tabulate$Male
  
  DF_subset_tabulate <- DF_subset_tabulate[order(DF_subset_tabulate$Total),]
  
  Output.table <- tableGrob(DF_subset_tabulate, rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # Save to local computer
  #----------------------------------------------
  tiff(filename = paste(outdir, Syapse_Export_timestamp,"_syapse_GenderAgeDistribution_",
                        cohort,".tiff", sep=""),
       width = 15, height = 7.5, units = "in", res = 200)
  grid.arrange(plot, Output.table, 
               heights = c(2, 0.4), ncol = 1, nrow = 2)
  dev.off()
  
  remove(DF_subset,DF_subset_tabulate,Output.table,plot,row_No)
} 

## STAMP database: pathogenicity status + variant type
#----------------------------------------------
if (isTRUE(static.plots_FILTER)) {
  # Subset relevant columns
  DF_subset <- DF[,c("PatientID","VariantLabel","VariantGene","VariantHGVSCoding",
                     "VariantPathogenicityStatus","var.type")]
  
  # Specify VariantPathogenicityStatus
  #----------------------------------------------
  pathogenic <- c("Likely Pathogenic","Pathogenic")
  VUS <- c("Unknown","Unknown significance")
  benign <- c("Likely Benign")
  
  for (row_No in 1:nrow(DF_subset)) {
    if (DF_subset$VariantPathogenicityStatus[row_No] %in% pathogenic) {
      DF_subset$VariantPathogenicityStatus[row_No] <- "pathogenic"
    } else if (DF_subset$VariantPathogenicityStatus[row_No] %in% VUS) {
      DF_subset$VariantPathogenicityStatus[row_No] <- "unknown_significance"
    } else if (DF_subset$VariantPathogenicityStatus[row_No] %in% benign) {
      DF_subset$VariantPathogenicityStatus[row_No] <- "benign"
    }
  }
  
  # Specify var.type.general
  #----------------------------------------------
  InFrame <- c("Deletion","Delins","Duplication","Insertion")
  for (row_No in 1:nrow(DF_subset)) {
    if (isTRUE(DF_subset$var.type[row_No] %in% InFrame)) {
      DF_subset$var.type[row_No] <- paste("InFrame_", DF_subset$var.type[row_No], sep="")
    }
  }
  
  for (row_No in 1:nrow(DF_subset)) {
    if (isTRUE(DF_subset$var.type[row_No] == "Synonymous")) {
      DF_subset$var.type.general[row_No] <- "Synonymous"
      
    } else if (isTRUE(grepl("^c.-[[:digit:]]+[ATCG]>[ATCG]", DF_subset$VariantHGVSCoding[row_No]))) {
      DF_subset$var.type.general[row_No] <- "Upstream Mutation"
      
    } else if (isTRUE(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+.*", DF_subset$VariantHGVSCoding[row_No]))) {
      DF_subset$var.type.general[row_No] <- "Intronic Mutation"
      
    } else if (isTRUE(grepl("InFrame", DF_subset$var.type[row_No]))) {
      DF_subset$var.type.general[row_No] <- "In-Frame Mutation"
      
    } else if (isTRUE(grepl("Frameshift", DF_subset$var.type[row_No]))) {
      DF_subset$var.type.general[row_No] <- "Frameshift Mutation"
      
    } else if (isTRUE(DF_subset$var.type[row_No] == "SNV")) {
      DF_subset$var.type.general[row_No] <- "SNV"
    }
  }
  
  # Tabulate frequency
  #----------------------------------------------
  DF_subset_tabulate <- 
    data.frame(DF_subset %>% group_by(VariantPathogenicityStatus,var.type.general,var.type) %>% tally())
  colnames(DF_subset_tabulate) <- c("Pathogenicity.Status","Mutation.Type","Variant.Type","No.Entries")
  
  for (row_No in 1:nrow(DF_subset_tabulate)) {
    DF_subset_tabulate$No.Genes[row_No] <- 
      length(unique(
        DF_subset$VariantGene[
          DF_subset$VariantPathogenicityStatus == DF_subset_tabulate$Pathogenicity.Status[row_No] &
            DF_subset$var.type.general == DF_subset_tabulate$Mutation.Type[row_No] &
            DF_subset$var.type == DF_subset_tabulate$Variant.Type[row_No]]))
  }
  
  for (row_No in 1:nrow(DF_subset_tabulate)) {
    DF_subset_tabulate$No.Samples[row_No] <- 
      length(unique(
        DF_subset$PatientID[
          DF_subset$VariantPathogenicityStatus == DF_subset_tabulate$Pathogenicity.Status[row_No] &
            DF_subset$var.type.general == DF_subset_tabulate$Mutation.Type[row_No] &
            DF_subset$var.type == DF_subset_tabulate$Variant.Type[row_No]]))
  }
  
  # Calculate percentages 
  #----------------------------------------------
  DF_subset_tabulate$Percentage <-
    round((100 * DF_subset_tabulate$No.Entries / sum(DF_subset_tabulate$No.Entries)), digits = 2)
  
  # Change order of rows
  for (row_No in 1:nrow(DF_subset_tabulate)) {
    if (isTRUE(DF_subset_tabulate$Mutation.Type[row_No] %in% 
               c("Intronic Mutation","Synonymous","Upstream Mutation"))) {
      DF_subset_tabulate$Lollipop[row_No] <- "Not_Mapped"
    } else {
      DF_subset_tabulate$Lollipop[row_No] <- "Mapped"
    }
    
    if (isTRUE(DF_subset_tabulate$Pathogenicity.Status[row_No] == "pathogenic")) {
      DF_subset_tabulate$Candidate[row_No] <- "Match-able"
    } else {
      DF_subset_tabulate$Candidate[row_No] <- "Not_Match-able"
    }
  }
  DF_subset_tabulate <- DF_subset_tabulate[order(DF_subset_tabulate$Lollipop, decreasing = FALSE),]
  DF_subset_tabulate <- DF_subset_tabulate[order(DF_subset_tabulate$Candidate, decreasing = FALSE),]
  
  for (row_No in 1:nrow(DF_subset_tabulate)) {
    if (row_No == 1) {
      DF_subset_tabulate$ymax[row_No] = DF_subset_tabulate$Percentage[row_No] 
      DF_subset_tabulate$ymin[row_No] = 0
    } else {
      DF_subset_tabulate$ymax[row_No] = DF_subset_tabulate$Percentage[row_No] + DF_subset_tabulate$ymax[row_No -1] 
      DF_subset_tabulate$ymin[row_No] = DF_subset_tabulate$ymax[row_No -1] 
    }
  }
  
  # Calculate position
  #----------------------------------------------
  for (row_No in 1:nrow(DF_subset_tabulate)) {
    if (row_No == 1) {
      DF_subset_tabulate$pos[row_No] <- 0.5*DF_subset_tabulate$Percentage[row_No]
    } else {
      DF_subset_tabulate$pos[row_No] <- sum(DF_subset_tabulate$Percentage[c(1:(row_No -1))]) + 
        0.5*DF_subset_tabulate$Percentage[row_No]
    }
  }
  
  # Annotation edits
  #----------------------------------------------
  # Modify text
  DF_subset_tabulate$Variant.Type <- gsub("^Frameshift_", "", DF_subset_tabulate$Variant.Type)
  DF_subset_tabulate$Variant.Type <- gsub("^InFrame_", "", DF_subset_tabulate$Variant.Type)
  
  DF_subset_tabulate$Mutation.Type <- gsub("[[:blank:]]Mutation$", "", DF_subset_tabulate$Mutation.Type)
  
  # Specify label
  DF_subset_tabulate$plotLabel <- 
    paste(DF_subset_tabulate$Mutation.Type, DF_subset_tabulate$Variant.Type, sep=" ")
  
  # Modify text
  DF_subset_tabulate$plotLabel <- gsub("SNV SNV", "SNV", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("Synonymous Synonymous", "Synonymous", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("Frameshift", "Fs", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("In-Frame ", "", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("Deletion", "Del", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("Insertion", "Ins", DF_subset_tabulate$plotLabel)
  DF_subset_tabulate$plotLabel <- gsub("Duplication", "Dup", DF_subset_tabulate$plotLabel)
  
  # Generate pie plot
  #----------------------------------------------
  colourCount = length(unique(append(DF_subset_tabulate$Pathogenicity.Status, 
                                     append(DF_subset_tabulate$Mutation.Type, DF_subset_tabulate$Variant.Type))))
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  plot <- ggplot(DF_subset_tabulate) + 
    geom_rect(aes(fill=Variant.Type, ymax=ymax, ymin=ymin, xmax=9, xmin=6)) +
    geom_rect(aes(fill=Mutation.Type, ymax=ymax, ymin=ymin, xmax=6, xmin=3)) +
    geom_rect(aes(fill=Pathogenicity.Status, ymax=ymax, ymin=ymin, xmax=3, xmin=0)) +
    # coord_polar(theta="y",start=55, direction = 1) +
    
    labs(title = "STAMP Database: Variant Distribution",
         subtitle = paste("N = ", sum(DF_subset_tabulate$No.Entries), " entries (N = ", 
                          length(unique(DF_subset$PatientID)), " samples; N = ",
                          length(unique(DF_subset$VariantGene)), " genes)", sep="")) +
    
    # geom_text_repel(aes(x =9.85, y = pos, label = plotLabel),
    #                 direction='y', nudge_x = 6.5,
    #                 segment.size = 0.25, show.legend = TRUE, color = 'black', fontface = "bold") +
    
    scale_fill_manual(values = getPalette(colourCount)) +
    
    theme_minimal() +
    theme(aspect.ratio=1,
          
          legend.position="bottom",
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          
          plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=12),
          
          axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank())
  
  # Generate output table
  #----------------------------------------------
  Output.table <- tableGrob(DF_subset_tabulate[,c("Pathogenicity.Status","Mutation.Type","Variant.Type",
                                                  "No.Entries","No.Genes","No.Samples")], rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),base_size = 8))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # Save to local computer
  #----------------------------------------------
  tiff(filename = paste(outdir, Syapse_Export_timestamp, "_syapse_PathoVarTypeDistribution_",
                        cohort,".tiff", sep=""),
       width = 14, height = 10, units = "in", res = 150)
  grid.arrange(plot, Output.table, ncol=2, widths=c(1.5, 1))
  dev.off()
  
  remove(DF_subset,DF_subset_tabulate,Output.table,plot,benign,colourCount,
         InFrame,pathogenic,row_No,VUS,getPalette)
}

## Order volume plots
#----------------------------------------------
month.list <- c("2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08","2015-09","2015-10","2015-11","2015-12",
                "2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08","2016-09","2016-10","2016-11","2016-12",
                "2017-01","2017-02","2017-03","2017-04","2017-05","2017-06","2017-07","2017-08","2017-09","2017-10","2017-11","2017-12",
                "2018-01","2018-02","2018-03","2018-04","2018-05","2018-06","2018-07","2018-08","2018-09","2018-10","2018-11","2018-12",
                "2019-01","2019-02","2019-03","2019-04")

month.label <- c("Jan 2015","Feb 2015","Mar 2015","Apr 2015","May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
                 "Jan 2016","Feb 2016","Mar 2016","Apr 2016","May 2016","Jun 2016","Jul 2016","Aug 2016","Sep 2016","Oct 2016","Nov 2016","Dec 2016",
                 "Jan 2017","Feb 2017","Mar 2017","Apr 2017","May 2017","Jun 2017","Jul 2017","Aug 2017","Sep 2017","Oct 2017","Nov 2017","Dec 2017",
                 "Jan 2018","Feb 2018","Mar 2018","Apr 2018","May 2018","Jun 2018","Jul 2018","Aug 2018","Sep 2018","Oct 2018","Nov 2018","Dec 2018",
                 "Jan 2019","Feb 2019","Mar 2019","Apr 2019")

## STAMP database: order volume distribution by assay
if (isTRUE(static.plots_FILTER)) {
  
  # Subset relevant columns
  DF_subset <- data.frame(patient.id = DF$PatientID, 
                          assay = DF$AssayName,
                          DateOrdered = DF$AssayReportDateReviewed,
                          stringsAsFactors = FALSE)
  # Remove duplicate entries 
  DF_subset <- DF_subset %>% dplyr::distinct(patient.id, .keep_all = TRUE)
  
  # Date format
  DF_subset$DateOrdered <- as.Date(DF_subset$DateOrdered, "%m/%d/%y")
  
  ## Subsitute for missing info 
  missing.id <- DF_subset$patient.id[is.na(DF_subset$DateOrdered)]
  for (elem_No in 1:length(missing.id)) {
    patient_id = missing.id[elem_No]
    DateOrdered.alt <- unique(DF$AssayDateReceived[which(DF$PatientID == patient_id)])
    
    DF_subset$DateOrdered[which(DF_subset$patient.id == patient_id)] <- DateOrdered.alt
  }
  
  # Convert to month/year
  DF_subset$DateOrdered <- gsub("([[:digit:]]{,4}[-][[:digit:]]{,2})(.*)","\\1", DF_subset$DateOrdered)
  
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_subset %>% group_by(DateOrdered,assay) %>% tally())
  colnames(DF_tabulate) <- c("AssayReportDateReviewed","AssayName","No.Samples")
  
  months.missing <- setdiff(month.list, unique(DF_tabulate$AssayReportDateReviewed))
  if (length(months.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(AssayReportDateReviewed = months.missing,
                                    AssayName = "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)",
                                    No.Samples = 0, stringsAsFactors = FALSE))
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$AssayReportDateReviewed, decreasing = FALSE),]
  
  # HISTOGRAM: Monthly
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Samples)/10)*10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 50
    } else {y_increment = 100
    }
  } else {y_increment = 250
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Samples) == 1)) {comment = "order"
  } else {comment = "orders"
  }
  
  # Reformat "AssayReportDateReviewed" 
  Month_Key <- data.frame(month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                          digit=seq(1,12,1), 
                          stringsAsFactors = FALSE)
  for (row_No in 1:nrow(DF_tabulate)) {
    row_Month <- Month_Key$month[which(Month_Key$digit == as.numeric(gsub("(^[[:digit:]]{4}[-])([[:digit:]]{2})","\\2",DF_tabulate$AssayReportDateReviewed[row_No])))]
    row_Year <- gsub("(^[[:digit:]]{4})(.*)","\\1",DF_tabulate$AssayReportDateReviewed[row_No])
    DF_tabulate$DateReviewed[row_No] <- paste(row_Month,row_Year,sep=" ")
  }
  DF_tabulate$DateReviewed <- factor(DF_tabulate$DateReviewed, levels = month.label)
  
  plot <- ggplot(DF_tabulate, aes(x=DateReviewed, y=No.Samples, fill=AssayName)) +
    geom_bar(stat="identity") +
    
    labs(title = "Monthly Order Volume",
         subtitle = paste("N = ", sum(DF_tabulate$No.Samples), " ", comment, sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Frequency of Orders", breaks = seq(0,round(ymax/100)*100,y_increment), 
                       limits=c(0,round(ymax/100)*100)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="bottom")
  
  # Save to local computer
  #----------------------------------------------
  tiff(filename = paste(outdir, Syapse_Export_timestamp,"_syapse_MonthOrderVolume_",
                        cohort,".tiff", sep=""),
       width = 15, height = 7.5, units = "in", res = 200)
  grid.arrange(plot)
  dev.off()
}

## Cumulative order volume plot
if (isTRUE(static.plots_FILTER)) {
  
  # Subset relevant columns
  DF_subset <- data.frame(patient.id = DF$PatientID, 
                          DateOrdered = DF$AssayReportDateReviewed,
                          stringsAsFactors = FALSE)
  # Remove duplicate entries 
  DF_subset <- DF_subset %>% dplyr::distinct(patient.id, .keep_all = TRUE)
  
  # Date format
  DF_subset$DateOrdered <- as.Date(DF_subset$DateOrdered, "%m/%d/%y")
  
  ## Subsitute for missing info 
  missing.id <- DF_subset$patient.id[is.na(DF_subset$DateOrdered)]
  for (elem_No in 1:length(missing.id)) {
    patient_id = missing.id[elem_No]
    DateOrdered.alt <- unique(DF$AssayDateReceived[which(DF$PatientID == patient_id)])
    
    DF_subset$DateOrdered[which(DF_subset$patient.id == patient_id)] <- DateOrdered.alt
  }
  
  # Convert to month/year
  DF_subset$DateOrdered <- gsub("([[:digit:]]{,4}[-][[:digit:]]{,2})(.*)","\\1", DF_subset$DateOrdered)
  
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_subset %>% group_by(DateOrdered) %>% tally())
  colnames(DF_tabulate) <- c("AssayReportDateReviewed","No.Samples")
  
  months.missing <- setdiff(month.list, unique(DF_tabulate$AssayReportDateReviewed))
  if (length(months.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(AssayReportDateReviewed = months.missing,
                                    No.Samples = 0, stringsAsFactors = FALSE))
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$AssayReportDateReviewed, decreasing = FALSE),]
  
  for (row_No in 1:nrow(DF_tabulate)) {
    if (isTRUE(row_No == 1)) {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Samples[row_No]
    } else {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Samples[row_No] + DF_tabulate$CumulativeCount[row_No -1]
    }
  }
  
  # HISTOGRAM: Monthly
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$CumulativeCount)/10)*10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 50
    } else {y_increment = 100
    }
  } else {y_increment = 250
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Samples) == 1)) {comment = "order"
  } else {comment = "orders"
  }
  
  # Reformat "AssayReportDateReviewed" 
  Month_Key <- data.frame(month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                          digit=seq(1,12,1), 
                          stringsAsFactors = FALSE)
  for (row_No in 1:nrow(DF_tabulate)) {
    row_Month <- Month_Key$month[which(Month_Key$digit == as.numeric(gsub("(^[[:digit:]]{4}[-])([[:digit:]]{2})","\\2",DF_tabulate$AssayReportDateReviewed[row_No])))]
    row_Year <- gsub("(^[[:digit:]]{4})(.*)","\\1",DF_tabulate$AssayReportDateReviewed[row_No])
    DF_tabulate$DateReviewed[row_No] <- paste(row_Month,row_Year,sep=" ")
  }
  DF_tabulate$DateReviewed <- factor(DF_tabulate$DateReviewed,levels = month.label)
  
  plot <- ggplot(DF_tabulate, aes(x=DateReviewed, y=CumulativeCount)) +
    geom_point(aes(y=CumulativeCount, color="red")) +
    geom_line(aes(y=CumulativeCount, group=1), color="red") +
    
    labs(title = "Cumulative Order Volume",
         subtitle = paste("N = ", sum(DF_tabulate$No.Samples), " ", comment, sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Frequency of Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="none")
  
  # Save to local computer
  #----------------------------------------------
  tiff(filename = paste(outdir, Syapse_Export_timestamp,"_syapse_CumOrderVolume_",
                        cohort,".tiff", sep=""),
       width = 15, height = 4, units = "in", res = 200)
  grid.arrange(plot)
  dev.off()
}

remove(DF,cohort)
cat(paste("Timestamp of data visualization generation FINISH: ", Sys.time(), sep=""),"\n")
