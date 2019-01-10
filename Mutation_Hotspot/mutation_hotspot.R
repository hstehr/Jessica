setwd("~/Documents/ClinicalDataScience_Fellowship/")

## FUNCTIONS
#----------------------------------------------
# Generate Lollipop Plot
Lollipop_Plot <- function(variant_type, assay,
                          static_plot = TRUE, dynamic_plot = FALSE,
                          loc.freq.max = max(as.numeric(data_var_FULL$var.freq))) {
  # Following variable must exist in global environment: data_var_FULL, DF_gene_INFO
  
  # Specify y-axis parameters = yaxis_max, y_tick, y_legend, y_legend_dynamic
  if (loc.freq.max <= 10) {     # ABCC9
    yaxis_max = 10
    y_tick = 1
    y_legend = -1
    y_legend_dynamic = -0.15
  } else if (loc.freq.max <= 20) {     # BRAF
    yaxis_max = 20
    y_tick = 2
    y_legend = -1.25
    y_legend_dynamic = -0.35
  } else if (loc.freq.max <= 50) {     # HTR1A
    yaxis_max = 50
    y_tick = 5
    y_legend = -5
    y_legend_dynamic = -1
  } else if (loc.freq.max <= 100) {     # EGFR
    yaxis_max = 100
    y_tick = 10
    y_legend = -10
    y_legend_dynamic = -1.75
  } else if (loc.freq.max <= 250) {     # KRAS
    yaxis_max = 250
    y_tick = 25
    y_legend = -25
    y_legend_dynamic = -6
  } else {     # KRAS - STAMP_all
    yaxis_max = loc.freq.max
    y_tick = 50
    y_legend = -50
    y_legend_dynamic =-10
  }
  text_remove <- paste("/>y_legend_dynamic: ", y_legend_dynamic, "<br ", sep="")

  # Specify x-axis parameters = x_tick
  if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 800) {
    x_tick = 25
  } else if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 1400) {
    x_tick = 50
  } else if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 2400) {
    x_tick = 100
  } else if (as.numeric(unique(DF_gene_INFO$gene.end)) > 2400) {
    x_tick = 200
  }
  
  if (isTRUE(static_plot)) {
    
    if (variant_type == "All Variants") {
      point_color = "black"
    } else if (variant_type == "Pathogenic Variant") {
      point_color = "firebrick4"
    } else if (variant_type == "Unknown Significance Variant") {
      point_color = "darkslategrey"
    } else if (variant_type == "Benign Variants") {
      point_color = "darkgreen"
    }
    
    plot_static <- 
      ggplot(data = data_var_FULL[which(data_var_FULL$pathogenicity.status %in% variant_type),]) +
      
      # Add points to indicate freq of mutation at each var.pos
      geom_segment(aes(x = var.pos, xend = var.pos, y = 0, yend = var.freq), 
                   color = "gray88", alpha=0.8) +
      geom_point(aes(x = var.pos, y = var.freq),
                 color = point_color, size=2.5, shape=16, alpha=0.8) +
      
      # Plot length of transcript
      geom_segment(data = DF_gene_INFO, aes(x = 1, xend = unique(DF_gene_INFO$gene.end), 
                                            y = y_legend, yend = y_legend), 
                   color = "gray88", size=6) +
      
      # Add domains of transcript
      geom_segment(data = DF_gene_INFO, aes(x = domain.start, xend = domain.end, 
                                            y = y_legend, yend = y_legend, 
                                            color = domain), size = 6) +
      
      # Labels
      ylab("No. Mutations") +
      labs(title = paste(variant_type, " variants (n=",
                         sum(data_var_FULL$var.freq[data_var_FULL$pathogenicity.status == variant_type]),
                         ")", sep="")) +
      
      # Scaling
      scale_fill_brewer(palette="Set2") +
      scale_x_continuous(breaks=seq(0, unique(DF_gene_INFO$gene.end), x_tick)) +
      scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend, yaxis_max)) +
      
      # Theme
      theme_light() +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.title.x=element_blank(),
            plot.title = element_text(hjust=0, vjust=0, face="bold", size=14)
      )
    
    assign("plot_static", plot_static, envir = .GlobalEnv)
  }
  
  if (isTRUE(dynamic_plot)) {
    
    plot_dynamic <- 
      ggplot(data = data_var_FULL) +
      # Add points to indicate freq of mutation at each var.pos
      geom_segment(aes(x = var.pos, xend = var.pos, y = 0, yend = var.freq, color = pathogenicity.status), 
                   alpha=0.8) +
      geom_point(aes(x = var.pos, y = var.freq,
                     color = pathogenicity.status), size=2.5, shape=16, alpha=0.8) +
      
      # Plot length of transcript
      geom_segment(data = DF_gene_INFO, aes(x = 1, xend = unique(DF_gene_INFO$gene.end), 
                                            y = y_legend_dynamic, yend = y_legend_dynamic), 
                   color = "gray88", size=6) +
      # Add domains of transcript
      geom_segment(data = DF_gene_INFO, aes(x = domain.start, xend = domain.end, 
                                            y = y_legend_dynamic, yend = y_legend_dynamic, 
                                            color = domain), size = 6) +
      
      # Labels
      ylab("No. Mutations") +
      xlab(paste("Variant Position (", unique(DF_gene_INFO$symbol), ")", sep="")) +
      labs(title = paste(unique(DF_gene_INFO$symbol), " variants from ", 
                         assay, " (n=",
                         sum(data_var_FULL$var.freq[data_var_FULL$pathogenicity.status == "All Variants"]),
                         " total)", sep="")) +
      
      # Scaling
      scale_fill_brewer(palette="Set2") +
      scale_x_continuous(breaks=seq(0, unique(DF_gene_INFO$gene.end), x_tick)) +
      scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend_dynamic, yaxis_max)) +
      
      # Theme
      theme_light() +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.title.x=element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0, face="bold", size=14)
      )
    
    
    plot_dynamic_int <- plotly_build(plot_dynamic)
    
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      # Customize name of traces
      if (isTRUE(grepl(",1)", plot_dynamic_int$x$data[[elem_No]]$name))) {
        plot_dynamic_int$x$data[[elem_No]]$name <- 
          gsub("^([(])(.*)", "\\2", plot_dynamic_int$x$data[[elem_No]]$name)
        plot_dynamic_int$x$data[[elem_No]]$name <- 
          gsub("([,]1[)])$", "", plot_dynamic_int$x$data[[elem_No]]$name)
      }
      
      # Customize hover text of domains 
      elem_No_sub <- which(grepl("domain.start", plot_dynamic_int$x$data[[elem_No]]$text))
      if (length(elem_No_sub) > 0) {
          plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
            gsub(text_remove, "", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
      }
      
      # Customize hover text of variants 
      elem_No_sub <- which(grepl("var.pos", plot_dynamic_int$x$data[[elem_No]]$text))
      if (length(elem_No_sub) > 0) {
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub("var.pos", "Amino Acid Position", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub("var.freq", "No. Mutations", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub("pathogenicity.status", "Variant Type", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub(".*y: 0<br />", "", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])  
      }
    }
    
    assign("plot_dynamic_int", plot_dynamic_int, envir = .GlobalEnv)
  }
}

## PIPELINE
#----------------------------------------------
Mutation_Pipeline <- function (gene.list, DF, assay) {
  # Following variables must exist in global environment: Genes, Domains
  
  # Generate info per unique gene
  #----------------------------------------------
  for (gene_num in 1:length(gene.list)) {
    gene_id = gene.list[gene_num]
    
    # Extract gene information
    #----------------------------------------------
    row_id = which(Domains$Name == gene_id)
    DF_gene_INFO = data.frame(
      # Extract start and end AA position of transcript
      symbol=gene_id, gene.start=1, gene.end=as.numeric(Genes[Genes$Name == gene_id, "AA.Length"]),
      # Extract domain information for transcript
      domain=as.factor(Domains[row_id, "Domain"]), domain.start=as.numeric(Domains[row_id, "Domain_Start"]), 
      domain.end=as.numeric(Domains[row_id, "Domain_End"]), 
      stringsAsFactors = FALSE)
    
    # Extract STAMP entries for gene of interest 
    #----------------------------------------------
    DF_patient_INFO <- DF[DF$base.gene == gene_id,]
    
    # Remove STAMP entries with UTR variants i.e. position > AA length of transcript 
    DF_remove <- DF_patient_INFO[which(DF_patient_INFO$var.position > unique(DF_gene_INFO$gene.end)),]
    if (nrow(DF_remove) != 0){
      DF_patient_INFO <- DF_patient_INFO[!DF_patient_INFO$sys.label %in% DF_remove$sys.label,]
    }
    remove(DF_remove,row_id)
    
    # Compile information for plotting STAMP entries 
    #----------------------------------------------
    # Extract variant position from smpl.hgvsProtein
    variant.loc.list <- DF_patient_INFO$var.position
    # Extract pathogenicity status from smpl.pathogenicityStatus
    variant.pathogenicity.list = DF_patient_INFO$smpl.pathogenicityStatus
    
    # Unique pathogenicity statuses to map
    index.sig = unique(variant.pathogenicity.list)
    pList = list()
    
    # Compile full DF per mapped gene 
    allele.count = table(variant.loc.list)
    data_var.all = data.frame(var.pos = as.numeric(names(allele.count)), 
                              var.freq = as.vector(allele.count),
                              pathogenicity.status = "All Variants")
    allele.count = table(variant.loc.list[which(variant.pathogenicity.list %in% pathogenic)])
    if (length(allele.count) > 0) {
      data_var.patho = data.frame(var.pos = as.numeric(names(allele.count)), 
                                  var.freq = as.vector(allele.count),
                                  pathogenicity.status = "Pathogenic Variant")
    }
    allele.count = table(variant.loc.list[which(variant.pathogenicity.list %in% vus)])
    if (length(allele.count) > 0) {
      data_var.vus = data.frame(var.pos = as.numeric(names(allele.count)), 
                                var.freq = as.vector(allele.count),
                                pathogenicity.status = "Unknown Significance Variant")
    }
    allele.count = table(variant.loc.list[which(variant.pathogenicity.list %in% benign)])
    if (length(allele.count) > 0) {
      data_var.benign = data.frame(var.pos = as.numeric(names(allele.count)), 
                                   var.freq = as.vector(allele.count),
                                   pathogenicity.status = "Benign Variants")
    }
    
    data_var_FULL <- data_var.all
    if (exists("data_var.patho")) {
      data_var_FULL <- rbind(data_var_FULL, data_var.patho)
      remove(data_var.patho)
    }
    if (exists("data_var.vus")) {
      data_var_FULL <- rbind(data_var_FULL, data_var.vus)
      remove(data_var.vus)
    }
    if (exists("data_var.benign")) {
      data_var_FULL <- rbind(data_var_FULL, data_var.benign)
      remove(data_var.benign)
    }
    remove(allele.count,data_var.all,DF_patient_INFO)
    data_var_FULL$var.pos <- as.numeric(data_var_FULL$var.pos)
    data_var_FULL$var.freq <- as.numeric(data_var_FULL$var.freq)
    
    assign("data_var_FULL", data_var_FULL, envir = .GlobalEnv)
    assign("DF_gene_INFO", DF_gene_INFO, envir = .GlobalEnv)
    
    # Generate plots for each pathogenicity status
    #----------------------------------------------
    if (as.numeric(table(DF$base.gene == gene_id)["TRUE"]) > 0) {
      Lollipop_Plot(variant_type = "All Variants", assay = assay,
                    dynamic_plot = TRUE)
      pNum = 1
      pList[[pNum]] <- plot_static
      
      if ("Pathogenic" %in% index.sig | "Likely Pathogenic" %in% index.sig) {
        Lollipop_Plot(variant_type = "Pathogenic Variant", assay = assay)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
      if ("Unknown significance" %in% index.sig | "Unknown" %in% index.sig) {
        Lollipop_Plot(variant_type = "Unknown Significance Variant", assay = assay)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
      if ("Likely Benign" %in% index.sig) {
        Lollipop_Plot(variant_type = "Benign Variants", assay = assay)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
      if (isTRUE(assay == "STAMP v1 and v2")) {
        assayName = "STAMP_all"
      } else {
        assayName = gsub("^(STAMP[[:blank:]])(.*)", "STAMP_\\2", assay)
      }
      
      if (isTRUE(saveStaticPlots)) {
        
        # Merge all static plots into single frame
        #----------------------------------------------
        gpanels <- ggarrange(plotlist = pList,
                             nrow = 4, ncol = 1,
                             legend = "bottom", common.legend = TRUE)
        
        # Save plot to local computer
        #----------------------------------------------
        gpanels <- annotate_figure(gpanels, 
                                   top = text_grob(paste(gene_id, " variants from ", assay, 
                                                         " (", unique(DF_gene_INFO$gene.end), " aa)", sep=""), 
                                                   face = "bold", size = 18))
        file_id = paste("Mutation_Hotspot/LollipopPlots/", gene_id, "_", assayName, ".jpg", sep="")
        
        ggexport(gpanels, filename = file_id, width = 4000, height = 4000, res=350)
        remove(gpanels)
      }
      
      # Output dynamic plot to cloud
      #----------------------------------------------
      if (isTRUE(saveDynamicPlots)) {
        p <- ggplotly(plot_dynamic_int)
        file_id_online = paste(gene_id, "_", assayName, sep="")
        
        if (assay == "STAMP v1") {
          filename_Full = paste("STAMP_v1_Lollipop/", file_id_online, sep="")
        } else if (assay == "STAMP v2") {
          filename_Full = paste("STAMP_v2_Lollipop/", file_id_online, sep="")
        } else {
          filename_Full = paste("STAMP_all_Lollipop/", file_id_online, sep="")
        }
        
        api_create(p, filename = filename_Full, 
                   fileopt = "overwrite", sharing = "private")
      }
    }
  }
}

## Read annotation files
#----------------------------------------------
# Import excel file as single list
STAMPv2_Annotation <- 
  import_list(paste("STAMP/", STAMP_annotation_timestamp, "_STAMPv2_Annotation.xlsx", sep=""), setclass = "tbl")

# Generate individual DF per worksheet in excel file
invisible(capture.output(lapply(names(STAMPv2_Annotation), 
                                function(x) assign(x, STAMPv2_Annotation[[x]], envir = .GlobalEnv))))

# Re-format protein domain INFO 
colnames(Domains) <- unlist(Domains[2,])
Domains <- data.frame(Domains[c(3:nrow(Domains)),c(2:ncol(Domains))])

# Re-format protein length INFO 
colnames(Genes) <- unlist(Genes[2,])
Genes <- data.frame(Genes[c(3:nrow(Genes)),c(2:ncol(Genes))])

DF <- read.csv(file = paste("Mutation_Hotspot/", Syapse_Export_timestamp, "_syapse_export_DF_STAMP_4Map.csv", sep=""),
               header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,sep = ",")
col.change <- c("smpl.hasOrderingPhysician", "sys.uniqueId", "smpl.gender", 
                "smpl.assayName", "base.chromosome", "smpl.specimenType", "smpl.specimenSite")
DF[ ,col.change] <- lapply(DF[ ,col.change], factor)

## Check missing gene INFO
#----------------------------------------------
gene.DF = data.frame(gene=sort(unique(Genes$Name)), stringsAsFactors=FALSE)
domain.DF = data.frame(gene=sort(unique(Domains$Name)), stringsAsFactors=FALSE)
gene.list <- data.frame(gene=unique(DF$base.gene), stringsAsFactors=FALSE)

Temp_DF <- anti_join(gene.list, gene.DF, by = "gene")
if (nrow(Temp_DF) > 0) {
  cat(paste("Transcripts with missing amino acid length:", sep=""), "\n")
  print(anti_join(gene.list, gene.DF, by = "gene"))
} else {
  cat(paste("All transcripts have annotated amino acid length.", sep=""), "\n")
}
Temp_DF <- anti_join(gene.list, domain.DF, by = "gene")
if (nrow(Temp_DF) > 0) {
  cat(paste("Transcripts with missing domain information:", sep=""), "\n")
  print(anti_join(gene.list, domain.DF, by = "gene"))
} else {
  cat(paste("All transcripts have annotated domain information.", sep=""), "\n")
}

# Categorize by assay data i.e. STAMP v1 and v2
#----------------------------------------------
DF_STAMP_v1 <- DF[DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" |
                    DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)", ]

DF_STAMP_v2 <- DF[DF$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)", ]

DF_STAMP_all <- DF

remove(Changes,Tiles,STAMPv2_Annotation, col.change, gene.DF, domain.DF,Temp_DF,gene.list,DF)

# Data Review 
#----------------------------------------------
cat(paste("DF_STAMP_v1: ", nrow(DF_STAMP_v1), " entries; ", length(unique(DF_STAMP_v1$base.gene)), 
          " genes; ", length(unique(DF_STAMP_v1$sys.uniqueId)), " patients.", sep=""), "\n")
cat(paste("DF_STAMP_v2: ", nrow(DF_STAMP_v2), " entries; ", length(unique(DF_STAMP_v2$base.gene)), 
          " genes; ", length(unique(DF_STAMP_v2$sys.uniqueId)), " patients.", sep=""), "\n")
cat(paste("DF_STAMP_all: ", nrow(DF_STAMP_all), " entries; ", length(unique(DF_STAMP_all$base.gene)), 
          " genes; ", length(unique(DF_STAMP_all$sys.uniqueId)), " patients.", sep=""), "\n")

# Generate Lollipop Plot
#----------------------------------------------
Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_v1$base.gene)), 
                  DF = DF_STAMP_v1,
                  assay = "STAMP v1")

Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_v2$base.gene)),
                  DF = DF_STAMP_v2,
                  assay = "STAMP v2")

Mutation_Pipeline(gene.list = sort(unique(DF_STAMP_all$base.gene)),
                  DF = DF_STAMP_all,
                  assay = "STAMP v1 and v2")

remove(plot_dynamic_int,plot_static,data_var_FULL,DF_gene_INFO)
remove(DF_STAMP_v1,DF_STAMP_v2,DF_STAMP_all,Lollipop_Plot,Mutation_Pipeline)
