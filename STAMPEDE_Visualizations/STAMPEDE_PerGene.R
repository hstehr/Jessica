sink(file = out.output, append = TRUE, split = FALSE)
options(max.print=999999)

# Specify key table
#----------------------------------------------
snv.list <- data.frame(gene=sort(unique(STAMP_DF$VariantGene)), stringsAsFactors = FALSE)

gene.list.total <- append(snv.list$gene,
                          append(STAMP_CNV$CNV_Gene,
                                 append(STAMP_Fusion$Gene1[which(STAMP_Fusion$Gene1 %in% fusion.gene.list.full)], 
                                        STAMP_Fusion$Gene2[which(STAMP_Fusion$Gene2 %in% fusion.gene.list.full)])))

gene.list.total <- sort(unique(gene.list.total))
assign("gene.list.total", gene.list.total, envir = .GlobalEnv)

# Check for missing gene INFO (SNV/Indels)
#----------------------------------------------
gene.DF = data.frame(gene=sort(unique(Genes$Name)), stringsAsFactors=FALSE)
domain.DF = data.frame(gene=sort(unique(Domains$Name)), stringsAsFactors=FALSE)

Temp_DF <- anti_join(snv.list, gene.DF, by = "gene")
if (nrow(Temp_DF) > 0) {
  cat(paste("Transcripts with missing amino acid length:", sep=""), "\n")
  print(anti_join(snv.list, gene.DF, by = "gene"))
} else {cat(paste("All SNV/Indel transcripts have annotated amino acid length.", sep=""), "\n")}

Temp_DF <- anti_join(snv.list, domain.DF, by = "gene")
if (nrow(Temp_DF) > 0) {
  cat(paste("Transcripts with missing domain information:", sep=""), "\n","\n")
  print(anti_join(snv.list, domain.DF, by = "gene"))
} else {cat(paste("All SNV/Indel transcripts have annotated domain information.", sep=""), "\n","\n")}

remove(Temp_DF,gene.DF,domain.DF)

# Iterate through each unique gene
#----------------------------------------------
for (gene_num in 1:length(gene.list.total)) {
  gene_id = gene.list.total[gene_num]
  cat(paste(gene_num,": ", gene_id, sep=""),"\n")
  
  # Subset DFs
  #----------------------------------------------
  gene_DF <- STAMP_DF[which(STAMP_DF$VariantGene == gene_id),]
  gene_Fusion <- STAMP_Fusion[which(grepl(gene_id,STAMP_Fusion$Fusion_Detail) == TRUE),]
  gene_CNV <- STAMP_CNV[which(STAMP_CNV$CNV_Gene == gene_id),]
  
  # Fusion info per site
  #----------------------------------------------
  if (nrow(gene_Fusion) > 0) {
    top_fusion_count_fxn (DF_Fusion = gene_Fusion, cohort = gene_id, outdir = outdir)
  }
  
  # Variant info per gene 
  #----------------------------------------------
  top_variant_count_fxn (DF_SNVIndel = gene_DF,
                         DF_Fusion = gene_Fusion,
                         DF_CNV = gene_CNV,
                         cohort = gene_id, outdir = outdir)
  
  variant_type_distribution_fxn (DF_SNVIndel = gene_DF, 
                                 DF_Fusion = gene_Fusion,
                                 DF_CNV = gene_CNV,
                                 cohort = gene_id, outdir = outdir)
  
  # Sample info per gene
  #----------------------------------------------
  col_extract <- c("PatientID","PrimaryTumorSite")
  top_site_count_fxn (DF = rbind(gene_DF[,col_extract],
                                 gene_Fusion[,col_extract],
                                 gene_CNV[!is.na(gene_CNV$PrimaryTumorSite),col_extract]),
                      cohort = gene_id, outdir = outdir)
  
  # Lollipop info per SNV/Indel gene
  #----------------------------------------------
  if (isTRUE(gene_id %in% unlist(snv.list$gene))) {
    # Map only SNVs to lollipop plots
    gene_DF_SNV <- gene_DF[which(gene_DF$var.type == "SNV"),]
    
    # Remove intronic mutations 
    gene_DF_SNV <- gene_DF_SNV[!(grepl("^c.(-)*[[:digit:]]+[-+]{1}[[:digit:]]+.*", gene_DF_SNV$VariantHGVSCoding)),]
    
    if (nrow(gene_DF_SNV) > 0) {
      
      if (isTRUE(gene_id %in% snv.hotspot.list)) {
        # Extract position info
        gene_DF_SNV$var.position <- gsub("(^c.)([-][[:digit:]]{,3})(.*)","\\2",gene_DF_SNV$VariantHGVSCoding)
        gene_DF_SNV$aa.start <- gsub("(^c.[-][[:digit:]]{,3})([[:alpha:]]{1})(.*)","\\2",gene_DF_SNV$VariantHGVSCoding)
        gene_DF_SNV$aa.end <- gsub("(^c..*[>}])(.*$)","\\2",gene_DF_SNV$VariantHGVSCoding)
        
        SNV_lollipop_fxn (DF = gene_DF_SNV, 
                          gene_id = gene_id,
                          promoter_plot = TRUE,
                          assay = "STAMP_v2",
                          outdir = outdir.lollipop)
        
      } else {
        # Mutation in genes that are !(Genes$Coverage == "promoter") are NOT mapped if outside domain
        SNV_lollipop_fxn (DF = gene_DF_SNV, 
                          gene_id = gene_id,
                          assay = "STAMP_v2",
                          outdir = outdir.lollipop)
      }
    }
    remove(gene_DF_SNV)
  }
  
  remove(gene_id,col_extract,gene_DF,gene_Fusion,gene_CNV)
}
remove(gene_num)
