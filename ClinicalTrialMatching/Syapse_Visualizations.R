# Generate single level pie graph
# Illustrate age and variant type distribution
# DF variable assigned as global variable in pipeline

cat(paste("Timestamp of data visualization generation START: ", Sys.time(), sep=""),"\n")

library("gridExtra")
library("gtable")
library("grid")
library("ggrepel")
library("RColorBrewer")

if (isTRUE(AgePlot.FILTER)) {
  
  #################################
  ## STAMP database: age distribution by gender
  #################################
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
    scale_y_continuous(name="Number of Individuals", breaks = seq(0,150,10)) + 
    scale_fill_discrete(name = "Gender") +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=15),
          
          legend.background = 
            element_rect(color = "black", fill = "white", size = 0.3, linetype = "solid"),
          legend.position = c(0.955, 0.88),
          legend.text=element_text(size=10),
          
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold")) +
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
  colnames(DF_subset_tabulate) <- c("Age.Cohort","Gender","No.Patients")
  
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

if (isTRUE(VariantPlot.FILTER)) {
  
  #################################
  ## STAMP database: pathogenicity status + variant type
  #################################
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
    DF_subset_tabulate$No.Patients[row_No] <- 
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
  
  pie <- ggplot(DF_subset_tabulate) + 
    geom_rect(aes(fill=Variant.Type, ymax=ymax, ymin=ymin, xmax=10, xmin=6)) +
    geom_rect(aes(fill=Mutation.Type, ymax=ymax, ymin=ymin, xmax=6, xmin=3)) +
    geom_rect(aes(fill=Pathogenicity.Status, ymax=ymax, ymin=ymin, xmax=3, xmin=0)) +
    coord_polar(theta="y",start=55, direction = 1) +
    
    labs(title = "STAMP Database: Variant Distribution",
         subtitle = paste("N = ", sum(DF_subset_tabulate$No.Entries), " entries (N = ", 
                          length(unique(DF_subset$PatientID)), " patients; N = ",
                          length(unique(DF_subset$VariantGene)), " genes)", sep="")) +
    
    geom_text_repel(aes(x =9.85, y = pos, label = plotLabel),
                    direction='y', nudge_x = 6.5,
                    segment.size = 0.25, show.legend = TRUE, color = 'black', fontface = "bold") +
    
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
                                                  "No.Entries","No.Genes","No.Patients")], rows = NULL,
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
       width = 14, height = 8, units = "in", res = 150)
  grid.arrange(pie, Output.table, ncol=2, widths=c(1.5, 1))
  dev.off()
  
  remove(DF_subset,DF_subset_tabulate,Output.table,pie,benign,colourCount,
         InFrame,pathogenic,row_No,VUS,getPalette)
}

remove(DF,cohort)
cat(paste("Timestamp of data visualization generation FINISH: ", Sys.time(), sep=""),"\n")
