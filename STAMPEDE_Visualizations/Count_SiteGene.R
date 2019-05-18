# TABLE: PrimaryTumorSite_No.Orders
#----------------------------------------------
Site_List <- data.frame(PrimaryTumorSite = sort(unique(site.list.total$PrimaryTumorSite)),
                        No.Cases = NA, 
                        stringsAsFactors = FALSE)

for (row_No in 1:nrow(Site_List)) {
  site_id = Site_List$PrimaryTumorSite[row_No]
  
  pt.list.SNV = unique(STAMP_DF$PatientID[which(STAMP_DF$PrimaryTumorSite == site_id)])
  pt.list.Fusion = unique(STAMP_Fusion$PatientID[which(STAMP_Fusion$PrimaryTumorSite == site_id)])
  pt.list.CNV = unique(STAMP_CNV$PatientID[which(STAMP_CNV$PrimaryTumorSite == site_id)])
  
  pt.list.total <- append(pt.list.SNV, append(pt.list.Fusion, pt.list.CNV))
  
  Site_List$No.Cases[row_No] = as.numeric(length(unique(pt.list.total)))
  
  remove(site_id,pt.list.SNV,pt.list.Fusion,pt.list.CNV,pt.list.total)
}

Site_List <- Site_List[order(Site_List$PrimaryTumorSite, decreasing = FALSE),]

write.table(Site_List, 
            file="~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/list_cancers.csv",
            append = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

remove(row_No,Site_List)

# TABLE: Gene_No.Occurrences
#----------------------------------------------
Gene_Summary <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/STAMPEDE_GeneName_Description.csv",
                         header = FALSE, stringsAsFactors = FALSE, sep = ",")
colnames(Gene_Summary) <- c("VariantGene","Summary")

GeneName_List <- data.frame(VariantGene = sort(unique(gene.list.total)),
                            No.Occurrences = NA,
                            stringsAsFactors = FALSE)

for (row_No in 1:nrow(GeneName_List)) {
  gene_id = GeneName_List$VariantGene[row_No]
  
  pt.list.SNV = unique(STAMP_DF$PatientID[which(STAMP_DF$VariantGene == gene_id)])
  pt.list.Fusion = unique(STAMP_Fusion$PatientID[grepl(gene_id, STAMP_Fusion$Fusion_Detail)])
  pt.list.CNV = unique(STAMP_CNV$PatientID[which(STAMP_CNV$CNV_Gene == gene_id)])
  
  pt.list.total <- append(pt.list.SNV, append(pt.list.Fusion, pt.list.CNV))
  
  GeneName_List$No.Occurrences[row_No] =  as.numeric(length(unique(pt.list.total)))
  
  remove(gene_id,pt.list.SNV,pt.list.Fusion,pt.list.CNV,pt.list.total)
}

GeneName_List <- left_join(GeneName_List, unique(Gene_Summary[,]), by = "VariantGene")

GeneName_List <- GeneName_List[order(GeneName_List$VariantGene, decreasing = FALSE),]

write.table(GeneName_List, 
            file="~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/list_genes.csv",
            append = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

remove(row_No,GeneName_List)
