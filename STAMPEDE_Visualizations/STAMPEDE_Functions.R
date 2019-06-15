## FUNCTIONS: incorporates only raw Syapse export for timeline mapping
#----------------------------------------------
# Static plots not generated 
TRF_volume_timeline_fxn <- function (DF, cohort, outdir, width, height) {
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
  
  # Save to cloud: cumulative plot
  if (isTRUE(saveDynamicPlots)) {
    file_id = paste("TRF_volume_cumulative_", cohort, sep="")
    
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$CumulativeCount)/10)*10
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 20)) {y_increment = 2
    } else if (isTRUE(ymax <= 50)) {y_increment = 5
    } else if (isTRUE(ymax <= 100)) {y_increment = 10
    } else if (isTRUE(ymax <= 250)) {y_increment = 25
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 250
    }
    
    plot_cumulative_dynamic <- ggplot(DF_tabulate[,c(4,3)], aes(x=DateReviewed, y=CumulativeCount)) +
      geom_point(aes(y=CumulativeCount, color="red")) +
      geom_line(aes(y=CumulativeCount, group=1), color="red") +
      
      labs(title = "Cumulative Order Volume Distribution",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " orders", sep="")) +
      
      xlab("") +
      scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                         limits=c(0,ymax)) + 
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text.y=element_text(size=14),
            axis.text.x=element_text(size=10,angle = 45, hjust = 1),
            axis.title=element_text(size=14,face="bold"),
            
            legend.position="none")
    
    plot_dynamic_int <- plotly_build(plot_cumulative_dynamic)
    
    # Customize hover text of domains 
    for (i in 1:length(plot_dynamic_int$x$data)) {
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("(^.*<br />colour: red<br />)(.*$)","\\2",plot_dynamic_int$x$data[[i]]$text)
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")  
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
  
  # Save to cloud: monthly plot
  if (isTRUE(saveDynamicPlots)) {
    file_id = paste("test_volume_monthly_", cohort, sep="")
    
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Orders)/10)*10
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 20)) {y_increment = 2
    } else if (isTRUE(ymax <= 50)) {y_increment = 5
    } else if (isTRUE(ymax <= 100)) {y_increment = 10
    } else if (isTRUE(ymax <= 250)) {y_increment = 25
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 250
    }
    
    plot_dynamic_bar <- ggplot(DF_tabulate[,c(4,2)], aes(x=DateReviewed, y=No.Orders)) +
      geom_bar(stat="identity",
               colour = custom.hues.75[1:nrow(DF_tabulate)], 
               fill  = custom.hues.75[1:nrow(DF_tabulate)]) +
      
      labs(title = "Monthly Order Volume Distribution",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " orders", sep="")) +
      
      xlab("") +
      scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                         limits=c(0,ymax)) + 
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text.y=element_text(size=14),
            axis.text.x=element_text(size=10,angle = 45, hjust = 1),
            axis.title=element_text(size=14,face="bold"),
            
            legend.position="none")
    
    plot_dynamic_int <- plotly_build(plot_dynamic_bar)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")  
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## FUNCTIONS: incorporates only Fusion DFs
#----------------------------------------------
fusion_count_fxn <- function (DF_Fusion, cohort, outdir) {
  
  if (nrow(DF_Fusion) > 0) {
    Fusion.1 <- DF_Fusion[,c("PatientID","Gene1")]
    Fusion.2 <- DF_Fusion[,c("PatientID","Gene2")]  
    
    colnames(Fusion.1) <- c("PatientID","VariantGene")
    colnames(Fusion.2) <- c("PatientID","VariantGene")
    
    DF <- rbind(Fusion.1,Fusion.2)
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF %>% group_by(VariantGene) %>% tally())
  colnames(DF_tabulate) <- c("Gene","No.Occurrences")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Occurrences, decreasing = TRUE),]
  
  # Extract listed genes 
  DF_tabulate <- DF_tabulate[which(DF_tabulate$Gene %in% fusion.gene.list.full),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
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
  DF_tabulate$Gene <- factor(DF_tabulate$Gene, 
                             levels = DF_tabulate$Gene[order(-DF_tabulate$No.Occurrences)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Occurrences)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height.table = 30
  width.table = 5
  
  height.plot = 15
  width.plot = 40
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Occurrences) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Gene, y=No.Occurrences, fill=Gene)) +
    geom_bar(stat="identity") +
    
    labs(title = "Quantification of Fusion Genes",
         subtitle = paste("N = ", sum(DF_tabulate$No.Occurrences), " ", comment, sep="")) +
    
    xlab("Gene Name") + 
    scale_y_continuous(name="Number of Occurrences", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 90, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.75)
  
  # Save to local computer
  file_id = paste("fusion_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,"_graph.tiff", sep=""),
         width = width.plot, height = height.plot, units = "in", res = 350)
    grid.arrange(plot)
    dev.off()    
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If all fusions have n=1>> output top 20 variants sorted alphabetically
## If at least 1 variant have (n > 1) 
## If less than 20 unique variants, output all variants
## If more than 20 unique variants
# > apply cutoff = round down value of 20th top variant to nearest 5
top_fusion_count_fxn <- function (DF_Fusion, cohort, outdir) {
  DF_Fusion_Fxn <- 
    DF_Fusion[which(DF_Fusion$Gene1 %in% fusion.gene.list.full | DF_Fusion$Gene2 %in% fusion.gene.list.full),
              c("PatientID","Fusion_Detail")]
  
  if (nrow(DF_Fusion_Fxn) > 0) {
    
    # TABLE
    #----------------------------------------------
    # Tabulate frequency
    DF_tabulate <- data.frame(DF_Fusion_Fxn %>% group_by(Fusion_Detail) %>% tally())
    colnames(DF_tabulate) <- c("Fusion","No.Mutations")
    DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Mutations, decreasing = TRUE),]
    
    continue_checkpoint <- NA
    if (isTRUE(max(DF_tabulate$No.Mutations) == 1)) {
      print(paste(cohort, ": all fusions have a frequency of  n=1", sep=""))
      cat("top_fusion_count plot displays top 20 fusion, as determined alphabetically.","\n")
      continue_checkpoint <- as.logical("FALSE")
      
    } else {
      if (isTRUE(nrow(DF_tabulate) > 20)) {
        cutoff = 5*floor(DF_tabulate$No.Mutations[[20]]/5) 
        
        if (isTRUE(cutoff < 2)) {
          DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= 2),]
        } else {
          DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= cutoff),]
        }
      }
    }
    
    if (isTRUE(continue_checkpoint == "FALSE")) {
      
      if (isTRUE(unique(DF_tabulate$No.Mutations) != "1")) {
        DF_tabulate <- DF_tabulate[order(DF_tabulate$Variant, decreasing = FALSE),]
      }
      
      max_row <- min(20, nrow(DF_tabulate))
      DF_tabulate <- DF_tabulate[1:max_row,]
    }
    
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
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
    DF_tabulate$Fusion <- factor(DF_tabulate$Fusion, 
                                 levels = DF_tabulate$Fusion[order(-DF_tabulate$No.Mutations)])
    
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Mutations)/10) * 10
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 20)) {y_increment = 2
    } else if (isTRUE(ymax <= 50)) {y_increment = 5
    } else if (isTRUE(ymax <= 100)) {y_increment = 10
    } else if (isTRUE(ymax <= 250)) {y_increment = 25
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 250
    }
    
    # Plot parameters
    height = 15
    if (nrow(DF_tabulate) <= 10) {
      if (nrow(DF_tabulate) <= 2) {width = 20
      } else {width = 25
      }
    } else {width = 30
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$No.Mutations) == 1)) {comment = "entry"
    } else {comment = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Fusion, y=No.Mutations, fill=Fusion)) +
      geom_bar(stat="identity") +
      
      labs(title = "Top Fusions",
           subtitle = paste("N = ", sum(DF_tabulate$No.Mutations), " / ", 
                            nrow(DF_Fusion), " total STAMP ", comment, sep="")) +
      
      xlab("Fusion") +
      scale_y_continuous(name="Number of Mutations", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text.y=element_text(size=14),
            axis.text.x=element_text(size=10,angle = 45, hjust = 1),
            axis.title=element_text(size=14,face="bold"),
            
            legend.position = "none") +
      
      scale_fill_manual(values = custom.hues.75)
    
    # Save to local computer
    file_id = paste("top_fusion_count_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      # Structure x-axis
      plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
      plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
      
      # Customize hover text
      for (elem_No in 1:length(plot_dynamic_int$x$data)) {
        plot_dynamic_int$x$data[[elem_No]]$text <- 
          gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
      }
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Top_Fusion_Count/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
}


## FUNCTIONS: incorporates only SNV/Indel DFs
#----------------------------------------------
pt_mutation_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","VariantPathogenicityStatus")])
  
  # TABLE - all variants
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
  DF_tabulate$No.Orders <- as.numeric("0")
  for (row_No in 1:nrow(DF_tabulate)) {
    DF_tabulate$No.Orders[row_No] <- 
      length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
  }
  DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 12))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM - all variants
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Orders)/10)*10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 10
  if (nrow(DF_tabulate) == 1) {width = 5.5
  } else {width = 5 + nrow(DF_tabulate) -1
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment1 = "test order"
  } else {comment1 = "test order"
  }
  
  if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Orders[[1]] == 1)) {comment2 = "entry"
  } else {comment2 = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Orders, fill=Mutation.Count)) +
    geom_bar(stat="identity") +
    
    labs(title = "Mutation Count",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                          nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
    
    xlab("Mutation Count") + 
    scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.5)
  
  # Save to local computer
  file_id = paste("pt_mutation_count_allvariants_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("Mutation.Count","Mutation Count",plot_dynamic_int$x$data[[elem_No]]$text)
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Patient_Mutation_Count/All_Variants/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
  
  # TABLE - all variants stacked pathogenic & VUS
  #----------------------------------------------
  # Reclassify pathogenicity status
  DF_Fxn$PathogenicityStatus <- NA
  for (row_No in 1:nrow(DF_Fxn)) {
    if (isTRUE(DF$VariantPathogenicityStatus[row_No] %in% c("Likely Pathogenic","Pathogenic"))) {
      DF_Fxn$PathogenicityStatus[row_No] <- "Pathogenic"
    } else if (isTRUE(DF$VariantPathogenicityStatus[row_No] %in% c("Unknown significance","Unknown"))) {
      DF_Fxn$PathogenicityStatus[row_No] <- "VUS"
    }
  } 
  
  DF_Fxn <- DF_Fxn[,c("PatientID","PathogenicityStatus")]
  
  DF_tabulate_stacked <- rbind(cbind(DF_tabulate,
                                     data.frame(Pathogenicity="Pathogenic", stringsAsFactors = FALSE)),
                               cbind(DF_tabulate,
                                     data.frame(Pathogenicity="VUS", stringsAsFactors = FALSE)))
  colnames(DF_tabulate_stacked) <- c("Mutation.Count","Total.No.Orders","Pathogenicity.Status")
  DF_tabulate_stacked$No.Orders <- as.numeric("0")
  
  max_No = max(as.numeric(unique(DF_tabulate$Mutation.Count)))
  
  for (mut_No in 1:max_No) {
    mutation_No = seq(1, max_No)[mut_No]
    
    mut.patient.list <- sort(unique(DF_tabulate_pre$PatientID[which(DF_tabulate_pre$n == mutation_No)]))
    
    if (isTRUE(length(mut.patient.list) > 0)) {
      DF_mut.patient <- DF_Fxn[which(DF_Fxn$PatientID %in% mut.patient.list),]
      DF_tabulate_pre_patho <- data.frame(DF_mut.patient %>% group_by(PathogenicityStatus) %>% tally())
      
      if (nrow(DF_tabulate_pre_patho) < 2) {
        DF_tabulate_pre_patho <- rbind(DF_tabulate_pre_patho,
                                       data.frame(PathogenicityStatus=
                                                    setdiff(c("Pathogenic","VUS"),
                                                            unique(DF_tabulate_pre_patho$PathogenicityStatus)),
                                                  n=0, stringsAsFactors = FALSE))  
      }
      
      # Populate DF
      DF_tabulate_stacked$No.Orders[which(DF_tabulate_stacked$Mutation.Count == mutation_No &
                                            DF_tabulate_stacked$Pathogenicity.Status == "Pathogenic")] <-
        DF_tabulate_pre_patho$n[DF_tabulate_pre_patho$PathogenicityStatus == "Pathogenic"]
      
      DF_tabulate_stacked$No.Orders[which(DF_tabulate_stacked$Mutation.Count == mutation_No &
                                            DF_tabulate_stacked$Pathogenicity.Status == "VUS")] <-
        DF_tabulate_pre_patho$n[DF_tabulate_pre_patho$PathogenicityStatus == "VUS"]
    }
  }
  
  # Calculate percentage of occurrence
  for (row_No in 1:nrow(DF_tabulate_stacked)) {
    DF_tabulate_stacked$Percent.Occurrence[row_No] <-
      as.numeric(round(100 * DF_tabulate_stacked$No.Orders[row_No] / 
                         (as.numeric(DF_tabulate_stacked$Mutation.Count[row_No]) * DF_tabulate_stacked$Total.No.Orders[row_No]),2))
  }
  
  # Convert NaN to "0"
  DF_tabulate_stacked$Percent.Occurrence[which(DF_tabulate_stacked$Percent.Occurrence == "NaN")] <- as.numeric(0)
  
  # Reorder columns 
  DF_tabulate_stacked <- DF_tabulate_stacked[,c("Mutation.Count","Pathogenicity.Status","Percent.Occurrence","Total.No.Orders")]
  
  Output.table <- tableGrob(DF_tabulate_stacked, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 12))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM - pathogenic
  #----------------------------------------------
  # Plot parameters
  height = 12
  if (nrow(DF_tabulate) == 1) {width = 5.5
  } else {width = 5 + nrow(DF_tabulate) -1
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment1 = "test order"
  } else {comment1 = "test orders"
  }
  
  if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Orders[[1]] == 1)) {comment2 = "entry"
  } else {comment2 = "entries"
  }
  
  plot <- ggplot(DF_tabulate_stacked, aes(x=Mutation.Count, y=Percent.Occurrence, 
                                          fill=Pathogenicity.Status)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    
    labs(title = "Mutation Count",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                          nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
    
    xlab("Mutation Count") + 
    ylab("Percent of Occurrences") +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    
    scale_fill_manual(values = custom.hues.2) +
    guides(fill=guide_legend(title="Variant Type"))
  
  # Save to local computer
  file_id = paste("pt_mutation_count_stacked_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.75), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      
      # Append total sample size
      elem_No_sub <- which(grepl("Mutation.Count", plot_dynamic_int$x$data[[elem_No]]$text))
      if (length(elem_No_sub) > 0) {
        for (i_sub in 1:length(elem_No_sub)) {
          
          var_id <- as.numeric(gsub("^Mutation.Count: ","",sub("<.*","", plot_dynamic_int$x$data[[elem_No]]$text[i_sub])))
          total.count = as.numeric(unique(DF_tabulate_stacked$Total.No.Orders[which(DF_tabulate_stacked$Mutation.Count == var_id)]))
          
          plot_dynamic_int$x$data[[elem_No]]$text[i_sub] <-
            gsub("$",paste("<br / >Total No.Orders: ",total.count,sep=""),plot_dynamic_int$x$data[[elem_No]]$text[i_sub])
        }
      }
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("Mutation.Count", "Mutation Count", plot_dynamic_int$x$data[[elem_No]]$text)
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("Percent.Occurrence", "Percent of Occurrence", plot_dynamic_int$x$data[[elem_No]]$text)
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<br />Pathogenicity.Status", "%<br />Pathogenicity Status", plot_dynamic_int$x$data[[elem_No]]$text)
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Patient_Mutation_Count/All_Variants_Stacked/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
  
  # TABLE - pathogenic
  #----------------------------------------------
  DF_Fxn <- unique(DF[which(DF$VariantPathogenicityStatus %in% c("Likely Pathogenic","Pathogenic")),
                      c("PatientID","VariantPathogenicityStatus")])
  
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  
  if (nrow(DF_tabulate_pre) > 0) {
    DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
    DF_tabulate$No.Orders <- as.numeric("0")
    for (row_No in 1:nrow(DF_tabulate)) {
      DF_tabulate$No.Orders[row_No] <- 
        length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
    }
    DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
    
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                                     base_size = 12))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 1, l = 1, r = ncol(Output.table))
    
    # HISTOGRAM - pathogenic
    #----------------------------------------------
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Orders)/10)*10
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 20)) {y_increment = 2
    } else if (isTRUE(ymax <= 50)) {y_increment = 5
    } else if (isTRUE(ymax <= 100)) {y_increment = 10
    } else if (isTRUE(ymax <= 250)) {y_increment = 25
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 250
    }
    
    # Plot parameters
    height = 10
    if (nrow(DF_tabulate) == 1) {width = 5.5
    } else {width = 5 + nrow(DF_tabulate) -1
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment1 = "test order"
    } else {comment1 = "test orders"
    }
    
    if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Orders[[1]] == 1)) {comment2 = "entry"
    } else {comment2 = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Orders, fill=Mutation.Count)) +
      geom_bar(stat="identity") +
      
      labs(title = "Mutation Count (Pathogenic)",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                            nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
      
      xlab("Mutation Count") + 
      scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            
            legend.position = "none") +
      
      scale_fill_manual(values = custom.hues.5)
    
    # Save to local computer
    file_id = paste("pt_mutation_count_pathovariants_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      # Structure x-axis
      plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
      plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
      
      # Customize hover text
      for (elem_No in 1:length(plot_dynamic_int$x$data)) {
        
        plot_dynamic_int$x$data[[elem_No]]$text <- 
          gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
        
        plot_dynamic_int$x$data[[elem_No]]$text <- 
          gsub("Mutation.Count","Mutation Count", plot_dynamic_int$x$data[[elem_No]]$text)
      }
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Patient_Mutation_Count/Pathogenic_Variants/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
  
  # TABLE - VUS
  #----------------------------------------------
  DF_Fxn <- unique(DF[which(DF$VariantPathogenicityStatus %in% c("Unknown significance","Unknown")),
                      c("PatientID","VariantPathogenicityStatus")])
  
  # Tabulate frequency
  DF_tabulate_pre <- data.frame(DF_Fxn %>% group_by(PatientID) %>% tally())
  
  if (nrow(DF_tabulate_pre) > 0) {
    DF_tabulate <- data.frame(Mutation.Count = seq(1,max(DF_tabulate_pre$n)))
    DF_tabulate$No.Orders <- as.numeric("0")
    for (row_No in 1:nrow(DF_tabulate)) {
      DF_tabulate$No.Orders[row_No] <- 
        length(DF_tabulate_pre$PatientID[DF_tabulate_pre$n == DF_tabulate$Mutation.Count[row_No]]) 
    }
    DF_tabulate$Mutation.Count <- as.factor(DF_tabulate$Mutation.Count)
    
    Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                              theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
                                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                                     base_size = 12))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
    Output.table <- gtable_add_grob(Output.table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                    t = 1, l = 1, r = ncol(Output.table))
    
    # HISTOGRAM - pathogenic
    #----------------------------------------------
    # Y-axis parameters
    ymax <- ceiling(max(DF_tabulate$No.Orders)/10)*10
    if (isTRUE(ymax <= 10)) {y_increment = 1
    } else if (isTRUE(ymax <= 20)) {y_increment = 2
    } else if (isTRUE(ymax <= 50)) {y_increment = 5
    } else if (isTRUE(ymax <= 100)) {y_increment = 10
    } else if (isTRUE(ymax <= 250)) {y_increment = 25
    } else if (isTRUE(ymax <= 500)) {y_increment = 50
    } else if (isTRUE(ymax <= 1000)) {y_increment = 100
    } else {y_increment = 250
    }
    
    # Plot parameters
    height = 10
    if (nrow(DF_tabulate) == 1) {width = 5.5
    } else {width = 5 + nrow(DF_tabulate) -1
    }
    
    # Subtitle parameters
    if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment1 = "test order"
    } else {comment1 = "test orders"
    }
    
    if (isTRUE(nrow(DF_tabulate) == 1 & DF_tabulate$No.Orders[[1]] == 1)) {comment2 = "entry"
    } else {comment2 = "entries"
    }
    
    plot <- ggplot(DF_tabulate, aes(x=Mutation.Count, y=No.Orders, fill=Mutation.Count)) +
      geom_bar(stat="identity") +
      
      labs(title = "Mutation Count (VUS)",
           subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " ", comment1, "; ",
                            nrow(DF_Fxn), " STAMP ", comment2, sep="")) +
      
      xlab("Mutation Count") + 
      scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
      
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
            plot.subtitle = element_text(hjust=1, face="bold",size=14),
            
            axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            
            legend.position = "none") +
      
      scale_fill_manual(values = custom.hues.5)
    
    # Save to local computer
    file_id = paste("pt_mutation_count_VUS_", cohort, sep="")
    
    if (isTRUE(saveStaticPlots)) {
      tiff(filename = paste(outdir, file_id,".tiff", sep=""),
           width = width, height = height, units = "in", res = 350)
      grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
      dev.off()
    }
    
    # Save to cloud
    if (isTRUE(saveDynamicPlots)) {
      plot_dynamic_int <- ggplotly(plot)
      
      # Autoscale x-axis
      plot_dynamic_int$x$layout$xaxis$autorange = TRUE
      
      # Structure x-axis
      plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
      plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
      
      # Customize hover text
      for (elem_No in 1:length(plot_dynamic_int$x$data)) {
        
        plot_dynamic_int$x$data[[elem_No]]$text <- 
          gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
        
        plot_dynamic_int$x$data[[elem_No]]$text <- 
          gsub("Mutation.Count","Mutation Count", plot_dynamic_int$x$data[[elem_No]]$text)
      }
      
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Patient_Mutation_Count/UnknownSignificance_Variants/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
}


## FUNCTIONS: incorporates SNV/Indels, Fusions, CNVs datasets
#----------------------------------------------
gene_count_fxn <- function (DF_SNVIndel, DF_Fusion, DF_CNV, cohort, outdir) {
  
  # Parse different DFs
  #----------------------------------------------
  colnames(DF_CNV)[which(colnames(DF_CNV) == "CNV_Gene")] <- "VariantGene"
  DF <- rbind(DF_SNVIndel[,c("PatientID","VariantGene")],
              DF_CNV[,c("PatientID","VariantGene")])
  
  # Re-format Fusion DF to have "VariantGene" be the STAMP v2 listed gene
  if (nrow(DF_Fusion) > 0) {
    Fusion.1 <- DF_Fusion[,c("PatientID","Gene1")]
    Fusion.2 <- DF_Fusion[,c("PatientID","Gene2")]  
    
    colnames(Fusion.1) <- c("PatientID","VariantGene")
    colnames(Fusion.2) <- c("PatientID","VariantGene")
    
    DF_Fusion <- rbind(Fusion.1,Fusion.2)
    DF_Fusion <- DF_Fusion[which(DF_Fusion$VariantGene %in% fusion.gene.list.full),]
    
    DF <- rbind(DF, DF_Fusion)
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF %>% group_by(VariantGene) %>% tally())
  colnames(DF_tabulate) <- c("Gene","No.Occurrences")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Occurrences, decreasing = TRUE),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
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
  DF_tabulate$Gene <- factor(DF_tabulate$Gene, 
                             levels = DF_tabulate$Gene[order(-DF_tabulate$No.Occurrences)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Occurrences)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height.table = 30
  width.table = 5
  
  height.plot = 15
  width.plot = 40
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Occurrences) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Gene, y=No.Occurrences, fill=Gene)) +
    geom_bar(stat="identity") +
    
    labs(title = "Quantification of Genes",
         subtitle = paste("N = ", sum(DF_tabulate$No.Occurrences), " ", comment, sep="")) +
    
    xlab("Gene Name") + 
    scale_y_continuous(name="Number of Occurrences", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=10,angle = 90, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") 
  
  # Save to local computer
  file_id = paste("gene_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,"_graph.tiff", sep=""),
         width = width.plot, height = height.plot, units = "in", res = 350)
    grid.arrange(plot) # List of 125 genes does not aesthetically fit into image
    dev.off()    
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If less than 20 unique genes, output all genes
## If more than 20 unique genes 
# > apply cutoff = round down value of 20th top gene to nearest 5
# > if cutoff is less than n=2, remove genes with n=1
top_gene_count_fxn <- function (DF_SNVIndel, DF_Fusion, DF_CNV, cohort, outdir) {
  
  # Parse different DFs
  #----------------------------------------------
  colnames(DF_CNV)[which(colnames(DF_CNV) == "CNV_Gene")] <- "VariantGene"
  DF_Fxn <- rbind(DF_SNVIndel[,c("PatientID","VariantGene")],
                  DF_CNV[,c("PatientID","VariantGene")])
  
  # Re-format Fusion DF to have "VariantGene" be the STAMP v2 listed gene
  if (nrow(DF_Fusion) > 0) {
    Fusion.1 <- DF_Fusion[,c("PatientID","Gene1")]
    Fusion.2 <- DF_Fusion[,c("PatientID","Gene2")]  
    
    colnames(Fusion.1) <- c("PatientID","VariantGene")
    colnames(Fusion.2) <- c("PatientID","VariantGene")
    
    DF_Fusion <- rbind(Fusion.1,Fusion.2)
    DF_Fusion <- DF_Fusion[which(DF_Fusion$VariantGene %in% fusion.gene.list.full),]
    
    DF_Fxn <- rbind(DF_Fxn, DF_Fusion)
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(VariantGene) %>% tally())
  colnames(DF_tabulate) <- c("Gene","No.Mutations")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Mutations, decreasing = TRUE),]
  
  if (isTRUE(nrow(DF_tabulate) > 20)) {
    cutoff = 5*floor(DF_tabulate$No.Mutations[[20]]/5) 
    
    if (isTRUE(cutoff < 2)) {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= 2),]
    } else {
      DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= cutoff),]
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
  DF_tabulate$Gene <- factor(DF_tabulate$Gene, levels = DF_tabulate$Gene[order(-DF_tabulate$No.Mutations)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Mutations)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 15
  if (nrow(DF_tabulate) <= 2) {width = 12
  } else if (nrow(DF_tabulate) <= 10) {width = 15
  } else {width = 20
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Mutations) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Gene, y=No.Mutations, fill=Gene)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Mutated Genes",
         subtitle = paste("N = ", sum(DF_tabulate$No.Mutations), " / ", 
                          (nrow(DF_SNVIndel) + nrow(DF_Fusion)/2 + nrow(DF_CNV)), 
                          " total STAMP ", comment, sep="")) +
    
    xlab("Gene") +
    scale_y_continuous(name="Number of Mutations", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.60)
  
  # Save to local computer
  file_id = paste("top_gene_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Top_Gene_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

variant_type_distribution_fxn <- function (DF_SNVIndel, DF_Fusion, DF_CNV, cohort, outdir) {
  
  if (nrow(DF_SNVIndel) > 0) {
    DF_SNVIndel_Fxn <- unique(DF_SNVIndel[,c("PatientID","var.type","VariantPathogenicityStatus")])
    
    # Reclassify variant type
    DF_SNVIndel_Fxn$VariantType <- NA
    for (row_No in 1:nrow(DF_SNVIndel_Fxn)) {
      if (isTRUE(DF_SNVIndel_Fxn$var.type[row_No] == "SNV")) {
        DF_SNVIndel_Fxn$VariantType[row_No] <- "SNV"
      } else if (isTRUE(DF_SNVIndel_Fxn$var.type[row_No] %in% 
                        c("Frameshift_Deletion","Frameshift_Delins","Frameshift_Duplication","Frameshift_Insertion"))) {
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
  
  # Append CNV data 
  DF_tabulate_full <- rbind(DF_tabulate_full, 
                            data.frame(VariantType="CNV",
                                       PathogenicityStatus="Amplification",
                                       No.Occurrences=0,
                                       Total.No.Occurrences=0,
                                       stringsAsFactors = FALSE))
  
  if (nrow(DF_CNV) > 0) {
    CNV.count = nrow(DF_CNV)
    
    DF_tabulate_full$No.Occurrences[which(DF_tabulate_full$VariantType == "CNV")] <- CNV.count
    DF_tabulate_full$Total.No.Occurrences[which(DF_tabulate_full$VariantType == "CNV")] <- CNV.count
  }
  
  # Append fusion data 
  DF_tabulate_full <- rbind(DF_tabulate_full, 
                            data.frame(VariantType="Fusion",
                                       PathogenicityStatus="Fusion",
                                       No.Occurrences=0,
                                       Total.No.Occurrences=0,
                                       stringsAsFactors = FALSE))
  
  if (nrow(DF_Fusion) > 0) {
    DF_Fusion <- DF_Fusion[which(DF_Fusion$Gene1 %in% fusion.gene.list.full | DF_Fusion$Gene2 %in% fusion.gene.list.full),]
    fusion.count = nrow(DF_Fusion)
    
    DF_tabulate_full$No.Occurrences[which(DF_tabulate_full$VariantType == "Fusion")] <- fusion.count
    DF_tabulate_full$Total.No.Occurrences[which(DF_tabulate_full$VariantType == "Fusion")] <- fusion.count
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
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
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
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    
    scale_fill_manual(values = custom.hues.4) +
    guides(fill=guide_legend(title="Variant Type"))
  
  # Save to local computer
  file_id = paste("var_type_distribution_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.5), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    for (i in 1:length(plot_dynamic_int$x$data)) {
      
      # Append total sample size
      elem_No_sub <- which(grepl("VariantType", plot_dynamic_int$x$data[[i]]$text))
      if (length(elem_No_sub) > 0) {
        for (i_sub in 1:length(elem_No_sub)) {
          
          var_id <- as.character(gsub("^VariantType: ","",sub("<.*","", plot_dynamic_int$x$data[[i]]$text[i_sub])))
          patho_id <- as.character(gsub("(.*>PathogenicityStatus: )(.*$)","\\2", plot_dynamic_int$x$data[[i]]$text[i_sub]))
          occur_id <- gsub("(^VariantType.*)(No.Occurrences.*)(<br />PathogenicityStatus:.*$)","\\2", plot_dynamic_int$x$data[[i]]$text[i_sub])
          
          total.count = as.numeric(unique(DF_tabulate_full$Total.No.Occurrences[which(DF_tabulate_full$VariantType == var_id)]))
          var_type <- paste(patho_id, var_id, sep=" ")
          
          plot_dynamic_int$x$data[[i]]$text[i_sub] <-
            paste("Variant Type: ",var_type,"<br />",occur_id,"<br />Total No.Occurrences: ",total.count,sep="")
          
          plot_dynamic_int$x$data[[i]]$text[i_sub] <- 
            gsub("Fusion Fusion","Fusion", plot_dynamic_int$x$data[[i]]$text[i_sub])
          
          plot_dynamic_int$x$data[[i]]$text[i_sub] <- 
            gsub("Amplification CNV","CNV Amplification", plot_dynamic_int$x$data[[i]]$text[i_sub])
        }
      }
      
      # Customize name of traces
      if (isTRUE(grepl(",1)", plot_dynamic_int$x$data[[i]]$name))) {
        plot_dynamic_int$x$data[[i]]$name <- 
          gsub("^([(])(.*)", "\\2", plot_dynamic_int$x$data[[i]]$name)
        
        plot_dynamic_int$x$data[[i]]$name <- 
          gsub("([,]1[)])$", "", plot_dynamic_int$x$data[[i]]$name)
      }
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Var_Type_Distribution/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If all variants have n=1 >> output top 20 variants sorted alphabetically
## If at least 1 variant have (n > 1) 
## If less than 20 unique variants, output all variants
## If more than 20 unique variants
# > apply cutoff = round down value of 20th top variant to nearest 5
top_variant_count_fxn <- function (DF_SNVIndel, DF_Fusion, DF_CNV, cohort, outdir) {
  
  DF_SNVIndel_Fxn <- DF_SNVIndel[,c("PatientID","VariantGene","VariantHGVSProtein")]
  DF_SNVIndel_Fxn$VariantDetail <- paste(DF_SNVIndel_Fxn$VariantGene, DF_SNVIndel_Fxn$VariantHGVSProtein, sep=" ")
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_SNVIndel_tabulate <- data.frame(DF_SNVIndel_Fxn %>% group_by(VariantDetail) %>% tally())
  colnames(DF_SNVIndel_tabulate) <- c("Variant","No.Mutations")
  DF_SNVIndel_tabulate <- DF_SNVIndel_tabulate[order(DF_SNVIndel_tabulate$No.Mutations, decreasing = TRUE),]
  
  # Append CNVs
  DF_CNV_tabulate <- data.frame(DF_CNV %>% group_by(CNV_Gene) %>% tally())
  colnames(DF_CNV_tabulate) <- c("Variant","No.Mutations")
  if (nrow(DF_CNV_tabulate) > 0) {
    DF_CNV_tabulate$Variant <- paste(DF_CNV_tabulate$Variant, "AMP",sep=" ")
  }
  
  # Append Fusions
  # Re-format Fusion DF to have "VariantGene" be the STAMP v2 listed gene
  if (nrow(DF_Fusion) > 0) {
    Fusion.1 <- DF_Fusion[,c("PatientID","Gene1")]
    Fusion.2 <- DF_Fusion[,c("PatientID","Gene2")]  
    
    colnames(Fusion.1) <- c("PatientID","VariantGene")
    colnames(Fusion.2) <- c("PatientID","VariantGene")
    
    DF_Fusion <- rbind(Fusion.1,Fusion.2)
    DF_Fusion <- DF_Fusion[which(DF_Fusion$VariantGene %in% fusion.gene.list.full),]
    
    DF_Fusion_tabulate <- data.frame(DF_Fusion %>% group_by(VariantGene) %>% tally())
    DF_Fusion_tabulate$VariantGene <- paste(DF_Fusion_tabulate$VariantGene, "Fusion",sep=" ")
    
  } else {
    DF_Fusion_tabulate <- data.frame(DF_Fusion %>% group_by(Gene1) %>% tally())
  }
  colnames(DF_Fusion_tabulate) <- c("Variant","No.Mutations")
  
  # Merge mutations
  DF_tabulate <- rbind(DF_SNVIndel_tabulate,DF_CNV_tabulate,DF_Fusion_tabulate)
  
  continue_checkpoint <- NA
  if (isTRUE(max(DF_tabulate$No.Mutations) == 1)) {
    print(paste(cohort, ": all variants have a frequency of  n=1", sep=""))
    cat("top_variant_count plot displays top 20 variants, as determined alphabetically.","\n")
    continue_checkpoint <- as.logical("FALSE")
    
  } else {
    if (isTRUE(nrow(DF_tabulate) > 20)) {
      cutoff = 5*floor(DF_tabulate$No.Mutations[[20]]/5) 
      
      if (isTRUE(cutoff < 2)) {
        DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= 2),]
      } else {
        DF_tabulate <- DF_tabulate[which(DF_tabulate$No.Mutations >= cutoff),]
      }
    }
  }
  
  if (isTRUE(continue_checkpoint == "FALSE")) {
    DF_tabulate <- DF_tabulate[order(DF_tabulate$Variant, decreasing = FALSE),]
    
    max_row <- min(20, nrow(DF_tabulate))
    DF_tabulate <- DF_tabulate[1:max_row,]
  }
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
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
  DF_tabulate$Variant <- factor(DF_tabulate$Variant, 
                                levels = DF_tabulate$Variant[order(-DF_tabulate$No.Mutations)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Mutations)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 15
  if (nrow(DF_tabulate) <= 10) {
    if (nrow(DF_tabulate) <= 2) {width = 20
    } else {width = 25
    }
  } else {width = 30
  }
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Mutations) == 1)) {comment = "entry"
  } else {comment = "entries"
  }
  
  plot <- ggplot(DF_tabulate, aes(x=Variant, y=No.Mutations, fill=Variant)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Variants",
         subtitle = paste("N = ", sum(DF_tabulate$No.Mutations), " / ", 
                          (nrow(DF_SNVIndel) + nrow(DF_Fusion)/2 + nrow(DF_CNV)), 
                          " total STAMP ", comment, sep="")) +
    
    xlab("Variant") +
    scale_y_continuous(name="Number of Mutations", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=10,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.75)
  
  # Save to local computer
  file_id = paste("top_variant_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Top_Variant_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}


## FUNCTIONS: incorporates SNV/Indels and Fusions datasets into single DF
#----------------------------------------------
# CNVs info is derived from SNV/Indels and Fusions datasets 
# Analysis of test order samples 

site_count_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","PrimaryTumorSite")])
  No.TotalOrders = length(unique(DF_Fxn$PatientID))
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(PrimaryTumorSite) %>% tally())
  colnames(DF_tabulate) <- c("PrimaryTumorSite","No.Orders")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders, decreasing = TRUE),]
  DF_tabulate$Percent.Orders <- as.numeric(round(100 * DF_tabulate$No.Orders/No.TotalOrders, 2))
  
  Output.table <- tableGrob(DF_tabulate[,c("PrimaryTumorSite","Percent.Orders")], rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.05)),
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
  DF_tabulate$PrimaryTumorSite <- factor(DF_tabulate$PrimaryTumorSite, 
                                         levels = DF_tabulate$PrimaryTumorSite[order(-DF_tabulate$Percent.Orders)])
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$Percent.Orders)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height.table = 12
  width.table = 5
  
  height.plot = 15
  width.plot = 25
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate$No.Orders) == 1)) {comment = "test order"
  } else {comment = "test orders"
  }
  
  plot_jpeg <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=Percent.Orders, fill=PrimaryTumorSite)) +
    geom_bar(stat="identity") +
    
    labs(title = "Primary Tumor Sites",
         subtitle = paste("N = ", sum(DF_tabulate$No.Orders), " ", comment, sep="")) +
    
    xlab("Primary Tumor Site") + 
    scale_y_continuous(name="Percent of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=30),
          plot.subtitle = element_text(hjust=1, face="bold",size=25),
          
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20,angle = 45, hjust = 1),
          axis.title=element_text(size=20,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.60)
  
  plot_interactive <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=Percent.Orders, fill=PrimaryTumorSite)) +
    geom_bar(stat="identity") +
    
    labs(title = "Primary Tumor Sites",
         subtitle = paste("N = ", sum(DF_tabulate$No.Orders), " ", comment, sep="")) +
    
    xlab("Primary Tumor Site") + 
    scale_y_continuous(name="Percent of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=10,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.60)
  
  # Save to local computer
  file_id = paste("site_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id, "_table.tiff", sep=""),
         width = width.table, height = height.table, units = "in", res = 350)
    grid.arrange(Output.table)
    dev.off()
    
    tiff(filename = paste(outdir, file_id,"_graph.tiff", sep=""),
         width = width.plot, height = height.plot, units = "in", res = 350)
    grid.arrange(plot_jpeg)
    dev.off()    
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    
    plot_dynamic_int <- ggplotly(plot_interactive)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("PrimaryTumorSite", "Primary Tumor Site", plot_dynamic_int$x$data[[elem_No]]$text)
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("$","%",gsub("Percent.Orders", "Percent of Total Orders", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

## If less than 20 unique sites, output all sites
## If more than 20 unique sites 
# > apply cutoff = round down value of 20th top site to nearest 5
# > if cutoff is less than n=2, remove sites with n=1
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
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
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
  
  plot_interactive <- ggplot(DF_tabulate, aes(x=PrimaryTumorSite, y=Percent.Orders, fill=PrimaryTumorSite)) +
    geom_bar(stat="identity") +
    
    labs(title = "Top Primary Tumor Sites",
         subtitle = paste("N = ", sum(DF_tabulate$No.Orders), " / ", nrow(DF_Fxn), " ", comment, sep="")) +
    
    xlab("Primary Tumor Site") +
    scale_y_continuous(name="Percent of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.30)
  
  # Save to local computer
  file_id = paste("top_site_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot_jpeg, Output.table, widths = c(2, 0.5), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot_interactive)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("PrimaryTumorSite", "Primary Tumor Site", plot_dynamic_int$x$data[[elem_No]]$text)
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("$","%",gsub("Percent.Orders", "Percent of Total Orders", plot_dynamic_int$x$data[[elem_No]]$text))
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Top_Site_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

gender_age_distribution_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","PatientGender","PatientAge")])
  
  # Plot parameters
  height = 7.5
  width = 15
  
  # Subtitle parameters
  if (isTRUE(nrow(DF_Fxn) == 1)) {comment = "test order"
  } else {comment = "test orders"
  }
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency = age of patient
  DF_tabulate_full <- data.frame(DF_Fxn %>% group_by(PatientAge,PatientGender) %>% tally())
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
  DF_Fxn$Age.Cohort <- NA
  for (row_No in 1:nrow(DF_Fxn)) {
    if (isTRUE(DF_Fxn$PatientAge[row_No] < 18)) {DF_Fxn$Age.Cohort[row_No] <- "Child (< 18yo)"
    } else if (isTRUE(DF_Fxn$PatientAge[row_No] < 65)) {DF_Fxn$Age.Cohort[row_No] <- "Adult (18-64yo)"
    } else {DF_Fxn$Age.Cohort[row_No] <- "Older Adult (>= 65yo)"
    }
  }
  
  # Tabulate frequency = age cohorts
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Age.Cohort,PatientGender) %>% tally())
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
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
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
  }
  
  DF_tabulate_full_edit <- rbind(DF_Keep,DF_Edit_Summary)
  DF_tabulate_full_edit$PatientAge <- factor(DF_tabulate_full_edit$PatientAge,
                                             levels=c(formatC(seq(0,89,1),0,format="f"),"90+"))
  
  plot <- ggplot(DF_tabulate_full_edit, aes(x=PatientAge, y=Relative.Frequency, fill=Gender)) +
    geom_bar(stat="identity") +
    
    labs(title = "Age and Gender Distribution",
         subtitle = paste("N = ", nrow(DF_Fxn), " ", comment, sep="")) +
    
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
  file_id = paste("gender_age_distribution_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.4), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    DF_tabulate_full_02 <- data.frame(DF_Fxn %>% group_by(PatientAge) %>% tally())
    
    # Customize hover text
    for (i in 1:length(plot_dynamic_int$x$data)) {
      
      # Append total sample size per age
      # for (age_No in 1:length(plot_dynamic_int$x$data[[i]]$text)) {
      #   age_elem <- gsub("(^PatientAge:[[:blank:]]+)([[:digit:]]+)(.*)", "\\2", plot_dynamic_int$x$data[[i]]$text[[age_No]])
      #   text_add <- paste("<br />Total.No.Orders: ", DF_tabulate_full_02$n[DF_tabulate_full_02$PatientAge == age_elem], sep="")
      #   
      #   plot_dynamic_int$x$data[[i]]$text[[age_No]] <- paste(plot_dynamic_int$x$data[[i]]$text[[age_No]], text_add, sep="")
      # }
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("PatientAge", "Patient Age", plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("Total.No.Orders","Total No.Orders",plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("Relative.Frequency","Percent of Total Orders",plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("<br />Gender","%<br />Gender",plot_dynamic_int$x$data[[i]]$text)
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Gender_Age_Distribution/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

specimen_type_stacked_fxn <- function(DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","smpl.specimenType")])
  colnames(DF_Fxn) <- c("PatientID","Specimen_Type")
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Specimen_Type) %>% tally())
  colnames(DF_tabulate) <- c("Specimen_Type","No.Orders")
  DF_tabulate <- DF_tabulate[order(DF_tabulate$No.Orders, decreasing = TRUE),]
  
  specimen.missing <- setdiff(c("formalin-fixed paraffin embedded tissue (FFPE)",
                                "bone marrow (BM)"),
                              unique(DF_tabulate$Specimen_Type))
  if (length(specimen.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Specimen_Type=specimen.missing,
                                    No.Orders=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate$Specimen_Type <- factor(DF_tabulate$Specimen_Type,
                                      levels = c("formalin-fixed paraffin embedded tissue (FFPE)",
                                                 "bone marrow (BM)"))
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Specimen_Type),]
  
  Output.table <- tableGrob(DF_tabulate, rows = NULL, 
                            theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.025)),
                                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                                   base_size = 11))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 2, b = nrow(Output.table), l = 1, r = ncol(Output.table))
  Output.table <- gtable_add_grob(Output.table,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(Output.table))
  
  # HISTOGRAM
  #----------------------------------------------
  DF_tabulate$Assay <- "STAMP_v2"
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Orders)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 10
  width = 5
  
  plot <- ggplot(DF_tabulate, aes(x=Assay, y=No.Orders, fill=Specimen_Type)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    
    labs(title = "Specimen Type Distribution",
         subtitle = paste("N = ", length(unique(DF$PatientID)), " specimen", sep="")) +
    
    xlab("") +
    scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position="bottom") +
    
    guides(fill=guide_legend(nrow=2,byrow=TRUE,title="Specimen Type")) +
    scale_fill_manual(values = custom.hues.2)
  
  # Save to local computer
  file_id = paste("specimen_type_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,"_stacked.tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, heights = c(2, 0.5), ncol = 1, nrow = 2)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    for (i in 1:length(plot_dynamic_int$x$data)) {
      
      # Customize name of traces
      if (isTRUE(grepl(",1)", plot_dynamic_int$x$data[[i]]$name))) {
        plot_dynamic_int$x$data[[i]]$name <- 
          gsub("^([(])(.*)", "\\2", plot_dynamic_int$x$data[[i]]$name)
        
        plot_dynamic_int$x$data[[i]]$name <- 
          gsub("([,]1[)])$", "", plot_dynamic_int$x$data[[i]]$name)
      }
      
      # Customize hover text of domains 
      for (elem_No in 1:length(plot_dynamic_int$x$data[[i]]$text)) {
        plot_dynamic_int$x$data[[i]]$text[[elem_No]] <- 
          gsub("^Assay: STAMP_v2<br />","", plot_dynamic_int$x$data[[i]]$text[[elem_No]])
        
        text_add <- paste("<br />Total.No.Orders: ", sum(DF_tabulate$No.Orders), sep="")
        
        plot_dynamic_int$x$data[[i]]$text[[elem_No]] <- paste(plot_dynamic_int$x$data[[i]]$text[[elem_No]], text_add, sep="")
        
        plot_dynamic_int$x$data[[i]]$text[[elem_No]] <- 
          gsub("Specimen_Type","Specimen Type", plot_dynamic_int$x$data[[i]]$text[[elem_No]])
      }
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Specimen_Type_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

tumor_purity_count_fxn <- function (DF, cohort, outdir, width, height) {
  DF_Fxn <- unique(DF[,c("PatientID","smpl.percentTumor")])
  
  # Round up tumor purity values 
  DF_Fxn$Tumor.Purity = 5*ceiling(DF_Fxn$smpl.percentTumor/5) 
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency
  DF_tabulate <- data.frame(DF_Fxn %>% group_by(Tumor.Purity) %>% tally())
  colnames(DF_tabulate) <- c("Tumor.Purity","No.Orders")
  
  percent.missing <- setdiff(seq(0,100,5), unique(DF_tabulate$Tumor.Purity))
  
  if (length(percent.missing) > 0) {
    DF_tabulate <- rbind(DF_tabulate,
                         data.frame(Tumor.Purity=percent.missing,
                                    No.Orders=0, stringsAsFactors = FALSE))  
  }
  
  DF_tabulate <- DF_tabulate[order(DF_tabulate$Tumor.Purity, decreasing = FALSE),]
  
  # Convert Tumor.Purity to categorical variable
  DF_tabulate$Tumor.Purity <- factor(DF_tabulate$Tumor.Purity,
                                     levels = seq(0,100,5))
  # HISTOGRAM
  #----------------------------------------------
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate$No.Orders)/10) * 10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 10
  width = 15
  
  plot <- ggplot(DF_tabulate, aes(x=Tumor.Purity, y=No.Orders, fill=Tumor.Purity)) +
    geom_bar(stat="identity") +
    
    labs(title = "Tumor Purity Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " specimen", sep="")) +
    
    xlab("Tumor Purity (rounded up)") +
    scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ymax,y_increment), limits=c(0, ymax)) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none") +
    
    scale_fill_manual(values = custom.hues.30)
  
  # Save to local computer
  file_id = paste("tumor_purity_count_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    print(plot)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- ggplotly(plot)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    # Customize hover text
    for (elem_No in 1:length(plot_dynamic_int$x$data)) {
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<$", "", gsub("(^.*<)(.*)", "\\1", plot_dynamic_int$x$data[[elem_No]]$text))
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("Tumor.Purity", "Tumor Purity", plot_dynamic_int$x$data[[elem_No]]$text)
      
      plot_dynamic_int$x$data[[elem_No]]$text <- 
        gsub("<", "%<", plot_dynamic_int$x$data[[elem_No]]$text)
      
    }
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Tumor_Purity_Count/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

# Missing histological info <NA> marked as "Not Applicable"
histologicaldx_distribution_fxn <- function (DF, cohort, outdir) {
  DF_Fxn <- unique(DF[,c("PatientID","HistologicalDx")])
  
  DF_Fxn$HistologicalDx[is.na(DF_Fxn$HistologicalDx)] <- "Not Applicable"
  
  # TABLE
  #----------------------------------------------
  # Tabulate frequency = age of patient
  DF_tabulate_full <- data.frame(DF_Fxn %>% group_by(HistologicalDx) %>% tally())
  colnames(DF_tabulate_full) <- c("HistologicalDx","No.Orders")
  
  DF_tabulate_full <- DF_tabulate_full[order(DF_tabulate_full$No.Orders, decreasing = TRUE),]
  
  Output.table <- tableGrob(DF_tabulate_full, rows = NULL,
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
  DF_tabulate_full$HistologicalDx <- factor(DF_tabulate_full$HistologicalDx, 
                                            levels = DF_tabulate_full$HistologicalDx[order(-DF_tabulate_full$No.Orders)])
  
  # Y-axis parameters
  ymax <- ceiling(max(DF_tabulate_full$No.Orders)/10)*10
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
  } else {y_increment = 250
  }
  
  # Plot parameters
  height = 7.5
  width = 10
  
  # Subtitle parameters
  if (isTRUE(sum(DF_tabulate_full$No.Orders) == 1)) {comment = "test order"
  } else {comment = "test order"
  }
  
  plot <- ggplot(DF_tabulate_full, aes(x=HistologicalDx, y=No.Orders, fill=HistologicalDx)) +
    geom_bar(stat="identity") +
    
    labs(title = "Histological Diagnosis Distribution",
         subtitle = paste("N = ", sum(DF_tabulate_full$No.Orders), " ", comment, sep="")) +
    
    xlab("Histological Diagnoses") +
    scale_y_continuous(name="Number of Test Orders", breaks = seq(0,ceiling(ymax/y_increment)*y_increment,y_increment), 
                       limits=c(0,ymax)) + 
    
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold",size=20),
          plot.subtitle = element_text(hjust=1, face="bold",size=14),
          
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle = 45, hjust = 1),
          axis.title=element_text(size=14,face="bold"),
          
          legend.position = "none")
  
  # Save to local computer
  file_id = paste("histdx_distribution_", cohort, sep="")
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot)
    
    # Customize hover text of domains 
    for (i in 1:length(plot_dynamic_int$x$data)) {
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("(^HistologicalDx.*)(<br />HistologicalDx:[[:blank:]]+.*$)","\\1", plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("^HistologicalDx","Histological Dx", plot_dynamic_int$x$data[[i]]$text)
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    filename_Full = paste("STAMPEDE/Hist_Dx_distribution/", file_id, sep="")
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}

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
  if (isTRUE(ymax <= 10)) {y_increment = 1
  } else if (isTRUE(ymax <= 20)) {y_increment = 2
  } else if (isTRUE(ymax <= 50)) {y_increment = 5
  } else if (isTRUE(ymax <= 100)) {y_increment = 10
  } else if (isTRUE(ymax <= 250)) {y_increment = 25
  } else if (isTRUE(ymax <= 500)) {y_increment = 50
  } else if (isTRUE(ymax <= 1000)) {y_increment = 100
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
  
  plot_dynamic_bar <- ggplot(DF_tabulate[,c(4,2,3)], aes(x=DateReviewed, y=No.Orders)) +
    geom_bar(stat="identity",
             colour = custom.hues.60[1:nrow(DF_tabulate)], 
             fill  = custom.hues.60[1:nrow(DF_tabulate)]) +
    
    geom_point(aes(y=CumulativeCount, color="red")) +
    geom_line(aes(y=CumulativeCount, group=1), color="red") +
    
    labs(title = "Order Volume Distribution",
         subtitle = paste("N = ", length(unique(DF_Fxn$PatientID)), " orders", sep="")) +
    
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
  
  if (isTRUE(saveStaticPlots)) {
    tiff(filename = paste(outdir, file_id,".tiff", sep=""),
         width = width, height = height, units = "in", res = 350)
    grid.arrange(plot, Output.table, widths = c(2, 0.75), ncol = 2, nrow = 1)
    dev.off()
  }
  
  # Save to cloud
  if (isTRUE(saveDynamicPlots)) {
    plot_dynamic_int <- plotly_build(plot_dynamic_bar)
    
    # Customize hover text of domains 
    for (i in 1:length(plot_dynamic_int$x$data)) {
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("<br />colour: red","",plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("(^CumulativeCount.*)(<br />No.Orders:[[:blank:]]+.*$)","\\1", plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("^CumulativeCount","Cumulative Count", plot_dynamic_int$x$data[[i]]$text)
      
      plot_dynamic_int$x$data[[i]]$text <- 
        gsub("DateReviewed","Date Reviewed", plot_dynamic_int$x$data[[i]]$text)
    }
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
    # Structure x-axis
    plot_dynamic_int$x$layout$xaxis$ticktext <- as.list(plot_dynamic_int$x$layout$xaxis$ticktext)
    plot_dynamic_int$x$layout$xaxis$tickvals <- as.list(plot_dynamic_int$x$layout$xaxis$tickvals)
    
    p <- ggplotly(plot_dynamic_int)
    
    if (isTRUE(PerSite)) {
      filename_Full = paste("STAMPEDE/PerTumor_OrderVolume/", file_id, sep="")
    } else {
      filename_Full = paste("STAMPEDE/", file_id, sep="")  
    }
    
    api_create(p, filename = filename_Full, 
               fileopt = "overwrite", sharing = "public")
  }
}


## FUNCTIONS: incorporates only SNV/Indels datasets to generate lollipop plots
#----------------------------------------------
Lollipop_Plot <- function(variant_type, assay,
                          static_plot = TRUE, dynamic_plot = FALSE,
                          promoter_plot = FALSE) {
  # Following variable must exist in global environment: data_var_FULL, DF_gene_INFO
  
  # Specify y-axis parameters = yaxis_max, y_tick, y_legend, y_legend_dynamic
  # Y-axis parameters
  yaxis_max <- ceiling(max(as.numeric(data_var_FULL$var.freq))/10)*10
  if (isTRUE(yaxis_max <= 10)) {y_tick = 1
  } else if (isTRUE(yaxis_max <= 20)) {y_tick = 2
  } else if (isTRUE(yaxis_max <= 50)) {y_tick = 5
  } else if (isTRUE(yaxis_max <= 100)) {y_tick = 10
  } else if (isTRUE(yaxis_max <= 250)) {y_tick = 25
  } else if (isTRUE(yaxis_max <= 500)) {y_tick = 50
  } else if (isTRUE(yaxis_max <= 1000)) {y_tick = 100
  } else {y_tick = 250
  }
  
  if (yaxis_max <= 10) {     # ABCC9
    y_legend = -1
    y_legend_dynamic = -0.15
  } else if (yaxis_max <= 20) {     # BRAF
    y_legend = -1.25
    y_legend_dynamic = -0.35
  } else if (yaxis_max <= 50) {     # HTR1A
    y_legend = -5
    y_legend_dynamic = -1
  } else if (yaxis_max <= 100) {     # EGFR
    y_legend = -10
    y_legend_dynamic = -1.75
  } else if (yaxis_max <= 250) {     # KRAS
    y_legend = -25
    y_legend_dynamic = -6
  } else {     # KRAS - STAMP_all
    y_legend = -50
    y_legend_dynamic =-10
  }
  
  text_remove <- paste("/>y_legend_dynamic: ", y_legend_dynamic, "<br ", sep="")
  
  # Specify x-axis parameters = x_tick
  if (isTRUE(promoter_plot)) {
    x_tick = 25 
    x_end_promoter = floor(min(data_var_FULL$var.pos)/10)*10
    
  } else {
    if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 800) { x_tick = 25
    } else if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 1400) { x_tick = 50
    } else if (as.numeric(unique(DF_gene_INFO$gene.end)) <= 2400) { x_tick = 100
    } else if (as.numeric(unique(DF_gene_INFO$gene.end)) > 2400) { x_tick = 200
    }
  }
  
  if (isTRUE(static_plot)) {
    
    if (variant_type == "All Variants") { point_color = "black"
    } else if (variant_type == "Pathogenic Variant") { point_color = "firebrick4"
    } else if (variant_type == "Unknown Significance Variant") { point_color = "darkslategrey"
    } else if (variant_type == "Benign Variants") { point_color = "darkgreen"
    }
    
    if (isTRUE(promoter_plot)) {
      plot_static <- 
        ggplot(data = data_var_FULL[which(data_var_FULL$pathogenicity.status %in% variant_type),]) +
        
        # Add points to indicate freq of mutation at each var.pos
        geom_segment(aes(x = var.pos, xend = var.pos, y = 0, yend = var.freq), 
                     color = "gray88", alpha=0.8) +
        geom_point(aes(x = var.pos, y = var.freq),
                   color = point_color, size=2.5, shape=16, alpha=0.8) +
        
        # Plot length of transcript
        geom_segment(data = DF_gene_INFO, aes(x = x_end_promoter, xend = 0, 
                                              y = y_legend, yend = y_legend), 
                     color = "gray88", size=2) +
        
        # Add domains of transcript
        geom_segment(data = DF_gene_INFO, aes(x = x_end_promoter, xend = domain.start, 
                                              y = y_legend, yend = y_legend, 
                                              color = domain), size = 2) +
        
        # Labels
        ylab("No. Mutations") +
        xlab("Amino Acid Position") +
        labs(title = paste(unique(DF_gene_INFO$symbol), " SNV Indels",sep="")) +
        
        # Scaling
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks=seq(x_end_promoter, 0, x_tick)) +
        scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend, yaxis_max)) +
        
        # Theme
        theme_light() +
        theme(legend.position = "bottom",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              
              plot.title = element_text(hjust=0.5, face="bold",size=20),
              
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(angle=0, hjust = 1),
              axis.title=element_text(size=14,face="bold")
        )
      
    } else {
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
                     color = "gray88", size=2) +
        
        # Add domains of transcript
        geom_segment(data = DF_gene_INFO, aes(x = domain.start, xend = domain.end, 
                                              y = y_legend, yend = y_legend, 
                                              color = domain), size = 2) +
        
        # Labels
        ylab("No. Mutations") +
        xlab("Amino Acid Position") +
        labs(title = paste(unique(DF_gene_INFO$symbol), " SNV Indels",sep="")) +
        
        # Scaling
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks=seq(0, unique(DF_gene_INFO$gene.end), x_tick)) +
        scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend, yaxis_max)) +
        
        # Theme
        theme_light() +
        theme(legend.position = "bottom",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              
              plot.title = element_text(hjust=0.5, face="bold",size=20),
              
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(angle=0, hjust = 1),
              axis.title=element_text(size=14,face="bold")
        )
    }
    
    assign("plot_static", plot_static, envir = .GlobalEnv)
  }
  
  if (isTRUE(dynamic_plot)) {
    
    if (isTRUE(promoter_plot)) {
      plot_dynamic <- 
        ggplot(data = data_var_FULL) +
        # Add points to indicate freq of mutation at each var.pos
        geom_segment(aes(x = var.pos, xend = var.pos, y = 0, yend = var.freq, color = pathogenicity.status), 
                     alpha=0.8) +
        geom_point(aes(x = var.pos, y = var.freq,
                       color = pathogenicity.status), size=2.5, shape=16, alpha=0.8) +
        
        # Plot length of transcript
        geom_segment(data = DF_gene_INFO, aes(x = x_end_promoter, xend = 0, 
                                              y = y_legend_dynamic, yend = y_legend_dynamic), 
                     color = "gray88", size=2) +
        # Add domains of transcript
        geom_segment(data = DF_gene_INFO, aes(x = x_end_promoter, xend = domain.start,
                                              y = y_legend_dynamic, yend = y_legend_dynamic, 
                                              color = domain), size = 2) +
        
        # Labels
        ylab("No. Mutations") +
        xlab("Amino Acid Position") +
        labs(title = paste(unique(DF_gene_INFO$symbol), " SNV Indels",sep="")) +
        
        # Scaling
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks=seq(x_end_promoter, 0, x_tick)) +
        scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend_dynamic, yaxis_max)) +
        
        # Theme
        theme_light() +
        theme(legend.position = "bottom",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              
              plot.title = element_text(hjust=0.5, face="bold",size=20),
              
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(angle=0, hjust = 1),
              axis.title=element_text(size=14,face="bold")
        )
      
    } else {
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
                     color = "gray88", size=2) +
        # Add domains of transcript
        geom_segment(data = DF_gene_INFO, aes(x = domain.start, xend = domain.end, 
                                              y = y_legend_dynamic, yend = y_legend_dynamic, 
                                              color = domain), size = 2) +
        
        # Labels
        ylab("No. Mutations") +
        xlab("Amino Acid Position") +
        labs(title = paste(unique(DF_gene_INFO$symbol), " SNV Indels",sep="")) +
        
        # Scaling
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks=seq(0, unique(DF_gene_INFO$gene.end), x_tick)) +
        scale_y_continuous(breaks=seq(0, yaxis_max, y_tick), limits=c(y_legend_dynamic, yaxis_max)) +
        
        # Theme
        theme_light() +
        theme(legend.position = "bottom",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              
              plot.title = element_text(hjust=0.5, face="bold",size=20),
              
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(angle=0, hjust = 1),
              axis.title=element_text(size=14,face="bold")
        )      
    }
    
    plot_dynamic_int <- plotly_build(plot_dynamic)
    
    # Autoscale x-axis
    plot_dynamic_int$x$layout$xaxis$autorange = TRUE
    
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
        
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub("domain.start", "Domain start position", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub(">domain.end", ">Domain end position", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub(">domain", ">Domain name", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
        
        # Remove underscore from domain name abbreviations
        plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub] <- 
          gsub("_", " ", plot_dynamic_int$x$data[[elem_No]]$text[elem_No_sub])
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

# Assume dataset does not contain benign mutations
SNV_lollipop_fxn <- function (DF, gene_id, assay, outdir, width, height,
                              promoter_plot = FALSE) {
  
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
  if (isTRUE(promoter_plot)) {
    # Filter for promoter region domain 
    DF_gene_INFO <- DF_gene_INFO[which(DF_gene_INFO$domain == "Promoter"),]
    
  } else {
    # Remove STAMP entries with UTR variants i.e. position > AA length of transcript 
    DF_remove <- DF[which(as.numeric(DF$var.position) > unique(DF_gene_INFO$gene.end)),]
    if (nrow(DF_remove) != 0){
      DF <- DF[!DF$VariantLabel %in% DF_remove$VariantLabel,]
    }
    remove(DF_remove,row_id)
  }
  
  # Compile information for plotting STAMP entries 
  #----------------------------------------------
  # Extract variant position from smpl.hgvsProtein
  variant.loc.list <- DF$var.position
  # Extract pathogenicity status from smpl.pathogenicityStatus
  variant.pathogenicity.list = DF$VariantPathogenicityStatus
  
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
  
  data_var_FULL <- data_var.all
  
  if (exists("data_var.patho")) {
    data_var_FULL <- rbind(data_var_FULL, data_var.patho)
    remove(data_var.patho)
  }
  
  if (exists("data_var.vus")) {
    data_var_FULL <- rbind(data_var_FULL, data_var.vus)
    remove(data_var.vus)
  }
  
  remove(allele.count,data_var.all)
  
  data_var_FULL$var.pos <- as.numeric(data_var_FULL$var.pos)
  data_var_FULL$var.freq <- as.numeric(data_var_FULL$var.freq)
  
  assign("DF_gene_INFO", DF_gene_INFO, envir = .GlobalEnv)
  assign("data_var_FULL", data_var_FULL, envir = .GlobalEnv)
  
  # Generate plots for each pathogenicity status
  #----------------------------------------------
  if (as.numeric(table(DF$VariantGene == gene_id)["TRUE"]) > 0) {
    
    if (isTRUE(promoter_plot)) {
      Lollipop_Plot(variant_type = "All Variants", 
                    assay = assay,
                    promoter_plot = TRUE,
                    dynamic_plot = TRUE)
      pNum = 1
      pList[[pNum]] <- plot_static
      
      if ("Pathogenic" %in% index.sig | "Likely Pathogenic" %in% index.sig) {
        Lollipop_Plot(variant_type = "Pathogenic Variant", 
                      assay = assay,
                      promoter_plot = TRUE)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
      if ("Unknown significance" %in% index.sig | "Unknown" %in% index.sig) {
        Lollipop_Plot(variant_type = "Unknown Significance Variant", 
                      assay = assay,
                      promoter_plot = TRUE)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
    } else {
      Lollipop_Plot(variant_type = "All Variants",
                    assay = assay,
                    dynamic_plot = TRUE)
      pNum = 1
      pList[[pNum]] <- plot_static
      
      if ("Pathogenic" %in% index.sig | "Likely Pathogenic" %in% index.sig) {
        Lollipop_Plot(variant_type = "Pathogenic Variant", 
                      assay = assay)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
      
      if ("Unknown significance" %in% index.sig | "Unknown" %in% index.sig) {
        Lollipop_Plot(variant_type = "Unknown Significance Variant", 
                      assay = assay)
        pNum = pNum + 1
        pList[[pNum]] <- plot_static
      }
    }
    
    file_id = paste(gene_id, "_", assay, sep="")
    
    # Save to local  >> plot generated in Lollipop_Plot fxn
    if (isTRUE(saveStaticPlots)) {
      
      # Plot parameters
      width.plot = 4000
      width.plot = 4000
      
      # Merge all static plots into single frame
      gpanels <- ggarrange(plotlist = pList, nrow = 3, ncol = 1,legend = "bottom", common.legend = TRUE)
      
      # Save plot to local computer
      gpanels <- annotate_figure(
        gpanels, top = text_grob(paste(gene_id, " SNV/Indel mutations from ", assay, " (", 
                                       unique(DF_gene_INFO$gene.end), " aa)", sep=""), 
                                 face = "bold", size = 18))
      
      ggexport(gpanels, filename =  paste(outdir.lollipop, file_id, ".jpg", sep=""),
               width = width.plot, height = width.plot, res=350)
      remove(gpanels)
    }
    
    # Save to cloud >> plots generated in Lollipop_Plot fxn
    if (isTRUE(saveDynamicPlots)) {
      p <- ggplotly(plot_dynamic_int)
      filename_Full = paste("STAMPEDE/Gene_Lollipop/", file_id, sep="")
      api_create(p, filename = filename_Full, 
                 fileopt = "overwrite", sharing = "public")
    }
  }
}

#################################
## Customized color palattes
#################################
## http://tools.medialab.sciences-po.fr/iwanthue/
## Parameters: n=XX soft (k-means), colorblind-friendly, hard (force vector repulsion algorithm)
## http://tools.medialab.sciences-po.fr/iwanthue/theory.php
custom.hues.75 = c("#c26f9d","#36ec97","#1c1b86","#dcd545","#0037a3","#7fc242","#9039ab","#3cc25d","#a22091","#2ff0ab",
                   "#b60074","#78eb97","#7a56cb","#849e00","#fd90ff","#9fe778","#740063","#01d8a6","#f450a3","#01873c",
                   "#eb58b9","#007738","#eb3f83","#01b088","#b8004f","#00e8fe","#a41f06","#6699ff","#c9a713","#0065c5",
                   "#facf5c","#004193","#bf8300","#0287d7","#ff8b45","#00ac8e","#c4173b","#7ee8bd","#9b0050","#9ada9d",
                   "#890051","#3b9e76","#fe5d5c","#006f45","#ff7bb9","#004d0a","#ffadfd","#414d00","#e4a8ff","#aa7000",
                   "#ffb9f7","#955800","#62255b","#f0d271","#5c0029","#ffb155","#6e001c","#ffb66d","#970033","#855821",
                   "#ff85af","#763700","#ff92a7","#691900","#ff6d8b","#974200","#ff8d8d","#651f11","#ff764c","#7f001a",
                   "#ff9577","#860a00","#ff6875","#d83f3b","#ff7966")

custom.hues.60 = c("#135300","#ff8af4","#51b949","#621488","#d1dd59","#00238a","#b6b520","#3a6ae1","#5d990d","#6f7ffa",
                   "#e0991e","#0060c6","#75ec95","#81007f","#00741f","#c155c4","#008b48","#de51b5","#009f69","#d0237c",
                   "#01e1e8","#c90b55","#01c4bc","#ad0032","#019368","#ff65bc","#094600","#a28eff","#9d8400","#1996ff",
                   "#b86200","#0159b2","#e5d56f","#30165a","#f4d07c","#00347a","#ff7b4c","#0191d2","#863800","#6cb3ff",
                   "#626f00","#be9aff","#8a984f","#930070","#7ea6ff","#8c0039","#0161a5","#ff6f96","#b6abff","#ac005f",
                   "#968cd3","#ff64aa","#654a8a","#ff97b5","#580054","#ffa7f9","#601041","#ff85cd","#78004d","#af535e")

custom.hues.30 = c("#47b041","#4e40b2","#9be876","#830074","#529100","#867cf7","#b7a606","#013c9f","#ded86d","#007ee5",
                   "#ed852d","#006baf","#af3109","#5ceacc","#dd3283","#00b973","#f95fbb","#004d0b","#fe93ff","#738700",
                   "#ce9fff","#685f00","#ff83da","#a96500","#ae75b2","#e9d387","#8a0038","#ff798b","#612500","#c71e3b")

custom.hues.5 = c("#2d1783","#a84c00","#a092ff","#b20049","#00326f")

custom.hues.4 = c("#bb7438","#7f64b9","#72ac5c","#b94b75")

custom.hues.3 = c("#f3a632","#60e9d9","#ba005b")

custom.hues.2 = c("#a2001e","#9acd46")
