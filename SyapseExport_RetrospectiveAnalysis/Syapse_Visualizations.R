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
  
  DF_subset <- unique(DF_subset[,c("patient.id","age","gender")])
  colnames(DF_subset) <- c("PatientID","PatientAge","PatientGender")
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency = age of patient
  DF_tabulate_full <- data.frame(DF_subset %>% group_by(PatientAge,PatientGender) %>% tally())
  colnames(DF_tabulate_full) <- c("PatientAge","Gender","No.Orders")
  row_append <- data.frame(PatientAge=setdiff(seq(1,100), unique(DF_tabulate_full$PatientAge)),
                           Gender=unique(DF_tabulate_full$Gender[[1]]),
                           No.Orders=0, stringsAsFactors = FALSE)
  DF_tabulate_full <- rbind(DF_tabulate_full,row_append)
  DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$PatientAge, decreasing = FALSE),]
  
  gender.missing <- setdiff(c("Female","Male"), unique(DF_tabulate_full$Gender))
  if (length(gender.missing) > 0) {
    DF_tabulate_full <- rbind(DF_tabulate_full,
                              data.frame(PatientAge=50, Gender=gender.missing,No.Orders=0, stringsAsFactors = FALSE))
  }
  
  # Specify Age.Cohort
  DF_subset$Age.Cohort <- NA
  for (row_No in 1:nrow(DF_subset)) {
    if (isTRUE(DF_subset$PatientAge[row_No] < 18)) {DF_subset$Age.Cohort[row_No] <- "Child (< 18yo)"
    } else if (isTRUE(DF_subset$PatientAge[row_No] < 65)) {DF_subset$Age.Cohort[row_No] <- "Adult (18-64yo)"
    } else {DF_subset$Age.Cohort[row_No] <- "Older Adult (>= 65yo)"
    }
  }
  
  # Tabulate frequency = age cohorts
  DF_tabulate <- data.frame(DF_subset %>% group_by(Age.Cohort,PatientGender) %>% tally())
  colnames(DF_tabulate) <- c("Age.Cohort","Gender","No.Orders")
  # Convert table to wide format
  DF_tabulate <- data.frame(cast(DF_tabulate, Age.Cohort ~ Gender), stringsAsFactors = FALSE)
  
  if (isTRUE(!("Male" %in% colnames(DF_tabulate)))) {DF_tabulate$Male <- as.numeric("0")}
  if (isTRUE(!("Female" %in% colnames(DF_tabulate)))) {DF_tabulate$Female <- as.numeric("0")} 
  DF_tabulate$Total <- DF_tabulate$Female + DF_tabulate$Male
  
  DF_tabulate <- DF_tabulate[,c("Age.Cohort","Female","Male","Total")]
  DF_tabulate$Female[which(is.na(DF_tabulate$Female))] <- 0
  DF_tabulate$Male[which(is.na(DF_tabulate$Male))] <- 0
  DF_tabulate$Total[which(is.na(DF_tabulate$Total))] <- 0
  
  cohort.missing <- setdiff(c("Adult (18-64yo)","Child (< 18yo)","Older Adult (>= 65yo)"),
                            unique(DF_tabulate$Age.Cohort))
  if (length(cohort.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Age.Cohort=cohort.missing,
                                    Female=0, Male=0, Total=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate$Age.Cohort <- factor(DF_tabulate$Age.Cohort,
                                   levels = c("Child (< 18yo)","Adult (18-64yo)","Older Adult (>= 65yo)"))
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Age.Cohort),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # Convert to relative frequency
  total.count = sum(DF_tabulate_full$No.Orders)
  DF_tabulate_full$Relative.Frequency <- as.numeric(round((100 * DF_tabulate_full$No.Orders) / total.count,2))
  
  # Y-axis parameters
  ymax <- c()
  for (row_No in 1:length(seq(0,100))) {
    ymax <- append(ymax, 
                   sum(DF_tabulate_full$Relative.Frequency[which(DF_tabulate_full$PatientAge == seq(0,100)[row_No])]))
  }
  ymax <- ceiling(max(ymax)/5)*5
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 10
    }
  } else {y_increment = 250
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$Total) == 1)) {comment = "test order"
  } else {comment = "test orders"
  }
  
  # HISTOGRAM
  #----------------------------------------------
  # Address age >90yo is PHI
  DF_Keep <- DF_tabulate_full[which(DF_tabulate_full$PatientAge < 90),]
  DF_Keep$PatientAge <- formatC(DF_Keep$PatientAge,0,format="f")
  
  DF_Edit <- DF_tabulate_full[which(DF_tabulate_full$PatientAge >= 90),]
  DF_Edit_Summary <- data.frame()
  gender.list <-c("Female","Male")
  for (elem_No in 1:length(gender.list)) {
    gender_id = gender.list[elem_No]
    
    DF_Edit_Summary <- rbind(DF_Edit_Summary,
                             data.frame(PatientAge="90+",
                                        Gender=gender_id,
                                        No.Orders=sum(DF_Edit$No.Orders[which(DF_Edit$Gender==gender_id)]),
                                        Relative.Frequency=sum(DF_Edit$Relative.Frequency[which(DF_Edit$Gender==gender_id)]),
                                        stringsAsFactors = FALSE))
    
    remove(elem_No,gender_id)
  }
  
  DF_tabulate_full_edit <- rbind(DF_Keep,DF_Edit_Summary)
  DF_tabulate_full_edit$PatientAge <- factor(DF_tabulate_full_edit$PatientAge,
                                             levels=c(formatC(seq(0,89,1),0,format="f"),"90+"))
  
  plot <- ggplot(DF_tabulate_full_edit, aes(x=PatientAge, y=Relative.Frequency, fill=Gender)) +
    geom_bar(stat="identity") +
    
    labs(title = "Age and Gender Distribution",
         subtitle = paste("N = ", nrow(DF_subset), " ", comment, sep="")) +
    
    scale_x_discrete(name="Age",breaks=c(as.character(seq(0,89,5)),"90+")) +
    scale_y_continuous(name="Percent of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    scale_fill_discrete(name = "Gender") +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.955, 0.88),
          legend.text=element_text(size=10),
          
          axis.text.x=element_text(size=14,angle = 0),
          axis.text.y=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  # Save to local computer
  #----------------------------------------------
  tiff(filename = paste(outdir, Syapse_Export_timestamp,"_syapse_GenderAgeDistribution_",
                        cohort,".tiff", sep=""),
       width = 15, height = 7.5, units = "in", res = 200)
  grid.arrange(plot, Output.table, 
               heights = c(2, 0.4), ncol = 1, nrow = 2)
  dev.off()
  
  remove(DF_subset,Output.table,plot,row_No,DF_Edit,DF_Edit_Summary,
         DF_Keep,DF_tabulate,DF_tabulate_full,DF_tabulate_full_edit,
         row_append,comment,total.count,cohort.missing,gender.list,
         gender.missing,ymax,y_increment)
}

## STAMP database: pathogenicity status + variant type
#----------------------------------------------
variant_type_distribution_fxn <- function (DF_SNVIndel, cohort, outdir) {
  
  if (nrow(DF_SNVIndel) > 0) {
    DF_SNVIndel_Fxn <- unique(DF_SNVIndel[,c("PatientID","var.type","VariantPathogenicityStatus")])
    
    # Reclassify variant type
    DF_SNVIndel_Fxn$VariantType <- NA
    for (row_No in 1:nrow(DF_SNVIndel_Fxn)) {
      if (isTRUE(DF_SNVIndel_Fxn$var.type[row_No] == "SNV")) {
        DF_SNVIndel_Fxn$VariantType[row_No] <- "SNV"
      } else if (isTRUE(DF_SNVIndel_Fxn$var.type[row_No] %in% 
                        c("Frameshift","Frameshift_Deletion","Frameshift_Delins","Frameshift_Duplication","Frameshift_Insertion"))) {
        DF_SNVIndel_Fxn$VariantType[row_No] <- "Frameshift Indel"
      } else if (isTRUE(DF_SNVIndel_Fxn$var.type[row_No] %in% 
                        c("Deletion","Delins","Duplication","Insertion"))) {
        DF_SNVIndel_Fxn$VariantType[row_No] <- "In-Frame Indel"
      }
    }
    
    # Reclassify pathogenicity status
    DF_SNVIndel_Fxn$PathogenicityStatus <- NA
    for (row_No in 1:nrow(DF_SNVIndel_Fxn)) {
      if (isTRUE(DF_SNVIndel$VariantPathogenicityStatus[row_No] %in% c("Likely Pathogenic","Pathogenic"))) {
        DF_SNVIndel_Fxn$PathogenicityStatus[row_No] <- "Pathogenic"
      } else if (isTRUE(DF_SNVIndel$VariantPathogenicityStatus[row_No] %in% c("Unknown significance","Unknown"))) {
        DF_SNVIndel_Fxn$PathogenicityStatus[row_No] <- "VUS"
      }
    }
    
    DF_SNVIndel_Fxn <- DF_SNVIndel_Fxn[,c("PatientID","VariantType","PathogenicityStatus")]
    
    # TABLE
    #----------------------------------------------
    # Tabulate frequency = age of patient
    DF_tabulate_full <- data.frame(DF_SNVIndel_Fxn %>% group_by(VariantType,PathogenicityStatus) %>% tally())
    colnames(DF_tabulate_full) <- c("VariantType","PathogenicityStatus","No.Occurrences")
    
    vartype.missing <- setdiff(c("SNV_Pathogenic",
                                 "SNV_VUS",
                                 "In-Frame Indel_Pathogenic",
                                 "In-Frame Indel_VUS",
                                 "Frameshift Indel_Pathogenic",
                                 "Frameshift Indel_VUS"),
                               unique(paste(DF_tabulate_full$VariantType,DF_tabulate_full$PathogenicityStatus,sep="_")))
    
    if (length(vartype.missing) > 0) {
      
      for (row_No in 1:length(vartype.missing)) {
        VariantType_add <- gsub("(^[[:alpha:]].*)(_)(.*)","\\1",vartype.missing[row_No])
        PathogenicityStatus_add <- gsub("(^[[:alpha:]].*)(_)([[:alpha:]].*)","\\3",vartype.missing[row_No])
        
        DF_tabulate_full <- rbind(DF_tabulate_full,
                                  data.frame(VariantType=VariantType_add,
                                             PathogenicityStatus=PathogenicityStatus_add,
                                             No.Occurrences=0, 
                                             stringsAsFactors = FALSE))
        
        remove(VariantType_add,PathogenicityStatus_add)
      }
    }
    
    DF_tabulate_full$VariantType <- factor(DF_tabulate_full$VariantType,
                                           levels = c("SNV","In-Frame Indel","Frameshift Indel"))
    DF_tabulate_full$PathogenicityStatus <- factor(DF_tabulate_full$PathogenicityStatus,
                                                   levels = c("Pathogenic","VUS"))
    
    DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$VariantType, decreasing = FALSE),]
    DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$PathogenicityStatus, decreasing = FALSE),]
    
    DF_tabulate_full$Total.No.Occurrences <- NA
    for (row_No in 1:nrow(DF_tabulate_full)) {
      var_id = DF_tabulate_full$VariantType[row_No]
      DF_tabulate_full$Total.No.Occurrences[row_No] <- 
        sum(as.numeric(DF_tabulate_full$No.Occurrences[which(DF_tabulate_full$VariantType == var_id)]))
    }
    
  } else {
    DF_tabulate_full <- data.frame(
      VariantType = c("SNV","In-Frame Indel","Frameshift Indel","SNV","In-Frame Indel","Frameshift Indel"),
      PathogenicityStatus = c("Pathogenic","Pathogenic","Pathogenic","VUS","VUS","VUS"),
      No.Occurrences=0,
      Total.No.Occurrences=0,
      stringsAsFactors = FALSE)
  }
  
  DF_tabulate_full$Variant.Type <- paste(DF_tabulate_full$VariantType,DF_tabulate_full$PathogenicityStatus, sep=" ")
  
  Output.table <- tableGrob(DF_tabulate_full[,1:4], rows = NULL,
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 10))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- c()
  VariantType.list <- as.character(unique(DF_tabulate_full$VariantType))
  for (row_No in 1:length(VariantType.list)) {
    ymax <- append(ymax, 
                   sum(DF_tabulate_full$No.Occurrences[which(DF_tabulate_full$VariantType == VariantType.list[row_No])]))
  }
  remove(VariantType.list)
  
  ymax <- ceiling(max(ymax)/10)*10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 10
    }
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 12
  width = 10
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate_full$No.Occurrences) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  DF_tabulate_full$VariantType <- factor(DF_tabulate_full$VariantType,
                                         levels = c("SNV","In-Frame Indel","Frameshift Indel","CNV","Fusion"))
  
  plot <- ggplot(DF_tabulate_full, aes(x=VariantType, y=No.Occurrences, fill=PathogenicityStatus)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    
    labs(title = "Variant Type Distribution",
         subtitle = paste("N = ", sum(DF_tabulate_full$No.Occurrences), " ", comment, sep="")) +
    
    scale_x_discrete(name="Variant Types") +
    scale_y_continuous(name="Number of Occurrences", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.92, 0.92),
          legend.text=element_text(size=10),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    
    scale_fill_manual(values = custom.hues.4) +
    guides(fill=guide_legend(title="Variant Type"))
  
  # Save to local computer
  file_id = paste("var_type_distribution_", cohort, sep="")
  
  tiff(filename = paste(outdir, file_id,".tiff", sep=""),
       width = width, height = height, units = "in", res = 350)
  grid.arrange(plot, Output.table, heights = c(2, 0.5), ncol = 1, nrow = 2)
  dev.off()
  
  remove(DF_tabulate_full,comment,ymax,y_increment)
}

if (isTRUE(static.plots_FILTER)) {
  
  variant_type_distribution_fxn(DF_SNVIndel = DF,
                                cohort = "all",
                                outdir = outdir)
}

## STAMP database: top primary tumor sites
#----------------------------------------------
top_site_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","PrimaryTumorSite")])
  No.TotalOrders = length(unique(DF_Fxn$PatientID))
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(PrimaryTumorSite) %>% tally())
  colnames(DF_tabulate) <- c("PrimaryTumorSite","No.Orders")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders, decreasing = TRUE),]
  
  DF_tabulate$Percent.Orders <- as.numeric(round(100 * DF_tabulate$No.Orders/No.TotalOrders, 2))
  
  if (isTRUE(nrow(DF_tabulate) > 20)) {
    cutoff = 5*floor(DF_tabulate$No.Orders[[20]]/5) 
    
    if (isTRUE(cutoff < 2)) {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Orders >= 2),]
    } else {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Orders >= cutoff),]
    }
  }
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 18))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$PrimaryTumorSite <- factor(DF_tabulate$PrimaryTumorSite, 
                                         levels = DF_tabulate$PrimaryTumorSite[order(-DF_tabulate$Percent.Orders)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$Percent.Orders)/10) * 10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 10
    }
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 15
  if (nrow(DF_tabulate) <= 2) {width = 12
  } else if (nrow(DF_tabulate) <= 10) {width = 20
  } else {width = 35
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment = "test order"
  } else {comment = "test orders"
  }
  
  plot_jpeg <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=Percent.Orders, fill=PrimaryTumorSite)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Primary Tumor Sites",
         subtitle = paste("N = ", sum(DF_tabulate$No.Orders), " / ", nrow(DF_Fxn), " ", comment, sep="")) +
    
    xlab("Primary Tumor Site") +
    scale_y_continuous(name="Percent of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=30),
          plot.subtitle = element_text(hjust=1, face="bold",size=25),
          
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20,angle = 45, hjust = 1),
          axis.title=element_text(size=20,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.30)
  
  # Save to local computer
  file_id = paste("top_site_count_", cohort, sep="")
  
  tiff(filename = paste(outdir, file_id,".tiff", sep=""),
       width = width, height = height, units = "in", res = 350)
  print(plot_jpeg)
  dev.off()
  
  remove(comment)
}

if (isTRUE(static.plots_FILTER)) {
  
  top_site_count_fxn(DF = DF,
                     cohort = "all",
                     outdir = outdir)
}

## STAMP database: test volume
#----------------------------------------------
# Specify distribution of timeline
#----------------------------------------------
month.number <- c("2014-01","2014-02","2014-03","2014-04","2014-05","2014-06","2014-07","2014-08","2014-09","2014-10","2014-11","2014-12",
                  "2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08","2015-09","2015-10","2015-11","2015-12",
                  "2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08","2016-09","2016-10","2016-11","2016-12",
                  "2017-01","2017-02","2017-03","2017-04","2017-05","2017-06","2017-07","2017-08","2017-09","2017-10","2017-11","2017-12",
                  "2018-01","2018-02","2018-03","2018-04","2018-05","2018-06","2018-07","2018-08","2018-09","2018-10","2018-11","2018-12",
                  "2019-01","2019-02","2019-03","2019-04")

month.alpha <- c("Jan 2014","Feb 2014","Mar 2014","Apr 2014","May 2014","Jun 2014","Jul 2014","Aug 2014","Sep 2014","Oct 2014","Nov 2014","Dec 2014",
                 "Jan 2015","Feb 2015","Mar 2015","Apr 2015","May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
                 "Jan 2016","Feb 2016","Mar 2016","Apr 2016","May 2016","Jun 2016","Jul 2016","Aug 2016","Sep 2016","Oct 2016","Nov 2016","Dec 2016",
                 "Jan 2017","Feb 2017","Mar 2017","Apr 2017","May 2017","Jun 2017","Jul 2017","Aug 2017","Sep 2017","Oct 2017","Nov 2017","Dec 2017",
                 "Jan 2018","Feb 2018","Mar 2018","Apr 2018","May 2018","Jun 2018","Jul 2018","Aug 2018","Sep 2018","Oct 2018","Nov 2018","Dec 2018",
                 "Jan 2019","Feb 2019","Mar 2019","Apr 2019")

test_volume_timeline_fxn <- function (DF, cohort, outdir, width, height, PerSite) {
  DF_Fxn <- unique(DF[,c("PatientID","AssayDateReceived")])
  # Convert to month/year
  DF_Fxn$AssayDateReceived <- gsub("-[[:digit:]]{2}$","", DF_Fxn$AssayDateReceived)
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(AssayDateReceived) %>% tally())
  colnames(DF_tabulate) <- c("AssayDateReceived","No.Orders")
  
  months.missing <- setdiff(month.number, unique(DF_tabulate$AssayDateReceived))
  if (length(months.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(AssayDateReceived=months.missing,No.Orders=0, stringsAsFactors = FALSE))
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$AssayDateReceived, decreasing = FALSE),]
  
  # Start dataframe with first entry
  start_No = min(which(DF_tabulate$No.Orders != "0"))
  DF_tabulate <- DF_tabulate[start_No:nrow(DF_tabulate),]
  
  for (row_No in 1:nrow(DF_tabulate)) {
    if (isTRUE(row_No == 1)) {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Orders[row_No]
    } else {
      DF_tabulate$CumulativeCount[row_No] <- DF_tabulate$No.Orders[row_No] + DF_tabulate$CumulativeCount[row_No -1]
    }
  }
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 12))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$CumulativeCount)/10)*10
  if (isTRUE(ymax < 1000)) {
    if (isTRUE(ymax <= 20)) {y_increment = 1
    } else if (isTRUE(ymax <= 30)) {y_increment = 2
    } else if (isTRUE(ymax <= 100)) {y_increment = 5
    } else if (isTRUE(ymax <= 200)) {y_increment = 10
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 10
    }
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 12
  width = 20
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment = "order"
  } else {comment = "orders"
  }
  
  # Reformat "AssayDateReceived" 
  Month_Key <- data.frame(month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                          digit=seq(1,12,1), 
                          stringsAsFactors = FALSE)
  
  for (row_No in 1:nrow(DF_tabulate)) {
    row_Month <- Month_Key$month[which(Month_Key$digit == as.numeric(gsub("(^[[:digit:]]{4}[-])([[:digit:]]{2})","\\2",DF_tabulate$AssayDateReceived[row_No])))]
    row_Year <- gsub("(^[[:digit:]]{4})(.*)","\\1",DF_tabulate$AssayDateReceived[row_No])
    DF_tabulate$DateReviewed[row_No] <- paste(row_Month,row_Year,sep=" ")
  }
  
  DF_tabulate$DateReviewed <- factor(DF_tabulate$DateReviewed,
                                     levels = month.alpha[start_No:length(month.alpha)])
  
  plot <- ggplot(DF_tabulate[,c(1:3)], aes(x=AssayDateReceived, y=No.Orders, fill=AssayDateReceived)) +
    geom_bar(stat="identity",
             colour = custom.hues.60[1:nrow(DF_tabulate)], 
             fill  = custom.hues.60[1:nrow(DF_tabulate)]) +
    
    geom_point(aes(y=CumulativeCount, color="red")) +
    geom_line(aes(y=CumulativeCount, group=1), color="red") +
    
    labs(title = "Order Volume Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment, sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="none")
  
  # Save to local computer
  file_id = paste("test_volume_", cohort, sep="")
  
  tiff(filename = paste(outdir, file_id,".tiff", sep=""),
       width = width, height = height, units = "in", res = 350)
  grid.arrange(plot, Output.table, widths = c(2, 0.75), ncol = 2, nrow = 1)
  dev.off()
  
  remove(comment)
}

if (isTRUE(static.plots_FILTER)) {
  
  test_volume_timeline_fxn(DF = DF,
                           cohort = "all",
                           outdir = outdir)
}

remove(DF,cohort,variant_type_distribution_fxn,test_volume_timeline_fxn,top_site_count_fxn,
       month.number,month.alpha)
cat(paste("Timestamp of data visualization generation FINISH: ", Sys.time(), sep=""),"\n")
