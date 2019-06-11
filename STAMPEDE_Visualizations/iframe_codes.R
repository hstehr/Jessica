#################################
## iFRAME CODE TABLE: MANUAL
#################################
# Import most recent mapping from Haik 20190517
#----------------------------------------------
mapping_HK <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/2019-0521_mappings_JC.csv", 
                       header = FALSE, stringsAsFactors = FALSE, sep = ",")
colnames(mapping_HK) <- c("Plot_Name","Folder","Label","iframe")

Folder_root = "STAMPEDE"

# ALL = DF_misc >> Files saved to STAMPEDE directory
#----------------------------------------------
# cohort="all" : site_count_fxn | gene_count_fxn | fusion_count_fxn | test_volume_timeline_fxn | TRF_volume_timeline_fxn
DF_misc <- data.frame(Plot_Name = c("site_count","gene_count","fusion_count",
                                    "test_volume","test_volume_cumulative","test_volume_monthly"),
                      Folder = Folder_root,
                      Label = "all",
                      iframe = c("2549","3328","5088",
                                 "3322","5333","5335"),
                      stringsAsFactors = FALSE)

# Per Gene = DF_PerGene
#----------------------------------------------
# top_site_count_fxn: cohort="all" & gene_DF
DF_site_count <- mapping_HK[which(mapping_HK$Plot_Name == "top_site_count" &
                                    mapping_HK$Folder == paste(Folder_root,"/Top_Site_Count",sep="")),]
DF_site_count <- DF_site_count[which(DF_site_count$Label %in% c("all",gene.list.total)),]
DF_site_count <- DF_site_count[!is.na(DF_site_count$iframe),]

# sort(setdiff(c("all",gene.list.total),DF_site_count$Label))
DF_site_count <- rbind(DF_site_count,
                       data.frame(Plot_Name = "top_site_count",
                                  Folder = paste(Folder_root,"/Top_Site_Count",sep=""),
                                  Label = sort(setdiff(c("all",gene.list.total),DF_site_count$Label)),
                                  iframe = c("5479","5485","5491","5499","5507"),
                                  stringsAsFactors = FALSE))

# variant_type_distribution_fxn: cohort="all" & gene_DF
DF_VariantType_Gene <- mapping_HK[which(mapping_HK$Plot_Name == "var_type_distribution" &
                                          mapping_HK$Folder == paste(Folder_root,"/Var_Type_Distribution",sep="")),]
DF_VariantType_Gene <- DF_VariantType_Gene[which(DF_VariantType_Gene$Label %in% c("all",gene.list.total)),]
DF_VariantType_Gene <- DF_VariantType_Gene[!is.na(DF_VariantType_Gene$iframe),]

# sort(setdiff(c("all",gene.list.total),DF_VariantType_Gene$Label))
DF_VariantType_Gene <- rbind(DF_VariantType_Gene,
                             data.frame(Plot_Name = "var_type_distribution",
                                        Folder = paste(Folder_root,"/Var_Type_Distribution",sep=""),
                                        Label = sort(setdiff(c("all",gene.list.total),DF_VariantType_Gene$Label)),
                                        iframe = c("5474","5483","5489","5497","5505"),
                                        stringsAsFactors = FALSE))

# top_variant_count_fxn: cohort=gene_DF
DF_variant_count <- mapping_HK[which(mapping_HK$Plot_Name == "top_variant_count" &
                                       mapping_HK$Folder == paste(Folder_root,"/Top_Variant_Count",sep="")),]
DF_variant_count <- DF_variant_count[which(DF_variant_count$Label %in% gene.list.total),]
DF_variant_count <- DF_variant_count[!is.na(DF_variant_count$iframe),]

# sort(setdiff(gene.list.total,DF_variant_count$Label))
DF_variant_count <- rbind(DF_variant_count,
                          data.frame(Plot_Name = "top_variant_count",
                                     Folder = paste(Folder_root,"/Top_Variant_Count",sep=""),
                                     Label = sort(setdiff(gene.list.total,DF_variant_count$Label)),
                                     iframe = c("5472","5481","5487","5495","5503"),
                                     stringsAsFactors = FALSE))

# Lollipop_Plot: cohort=gene_DF
snv.list <- sort(unique(STAMP_DF$VariantGene))

DF_lollipop <- mapping_HK[which(mapping_HK$Plot_Name == "gene_lollipop" &
                                  mapping_HK$Folder == paste(Folder_root,"/Gene_Lollipop",sep="")),]
DF_lollipop <- DF_lollipop[which(DF_lollipop$Label %in% snv.list),]
DF_lollipop <- DF_lollipop[!is.na(DF_lollipop$iframe),]

# sort(setdiff(snv.list,DF_lollipop$Label))
DF_lollipop <- rbind(DF_lollipop,
                     data.frame(Plot_Name = "gene_lollipop",
                                Folder = paste(Folder_root,"/Gene_Lollipop",sep=""),
                                Label = sort(setdiff(snv.list,DF_lollipop$Label)),
                                iframe = c("5478",NA,"5493","5501","5509"),
                                stringsAsFactors = FALSE))

# top_fusion_count_fxn: cohort=gene_DF
fusion.list <- sort(unique(append(STAMP_Fusion$Gene1, STAMP_Fusion$Gene2)))

DF_fusion_count <- mapping_HK[which(mapping_HK$Plot_Name == "top_fusion_count" &
                                  mapping_HK$Folder == paste(Folder_root,"/Top_Fusion_Count",sep="")),]
DF_fusion_count <- DF_fusion_count[which(DF_fusion_count$Label %in% fusion.list),]
DF_fusion_count <- DF_fusion_count[!is.na(DF_fusion_count$iframe),]

# sort(setdiff(fusion.list,DF_fusion_count$Label))
# DF_fusion_count <- rbind(DF_fusion_count,
#                          data.frame(Plot_Name = "top_fusion_count", 
#                                     Folder = paste(Folder_root,"/Top_Fusion_Count",sep=""),
#                                     Label = sort(setdiff(fusion.list,DF_fusion_count$Label)),
#                                     iframe = c("5400","5408","5416","5396","5424","5432","5440","5448","5456","5464"),
#                                     stringsAsFactors = FALSE))

DF_PerGene <- rbind(DF_site_count,DF_variant_count,DF_lollipop,DF_VariantType_Gene,DF_fusion_count)

remove(DF_site_count,DF_variant_count,DF_lollipop,DF_VariantType_Gene,DF_fusion_count,
       snv.list,fusion.list)

# Per PrimaryTumorSite
#----------------------------------------------
site.list.total <- site.list.total$PrimaryTumorSite

# top_variant_count_fxn: cohort="all" & site_DF
DF_VariantCount <- mapping_HK[which(mapping_HK$Plot_Name == "top_variant_count" &
                                      mapping_HK$Folder == paste(Folder_root,"/Top_Variant_Count",sep="")),]
DF_VariantCount <- DF_VariantCount[which(DF_VariantCount$Label %in% c("all",site.list.total)),]
DF_VariantCount <- DF_VariantCount[!is.na(DF_VariantCount$iframe),]

# sort(setdiff(c("all",site.list.total),DF_VariantCount$Label))
# DF_VariantCount <- rbind(DF_VariantCount,
#                          data.frame(Plot_Name = "top_variant_count",
#                                     Folder = paste(Folder_root,"/Top_Variant_Count",sep=""),
#                                     Label = sort(setdiff(c("all",site.list.total),DF_VariantCount$Label)),
#                                     iframe = c("5340","5346","5376","5390"),
#                                     stringsAsFactors = FALSE))

# top_gene_count_fxn: cohort="all" & site_DF
DF_GeneCount <- mapping_HK[which(mapping_HK$Plot_Name == "top_gene_count" &
                                   mapping_HK$Folder == paste(Folder_root,"/Top_Gene_Count",sep="")),]
DF_GeneCount <- DF_GeneCount[which(DF_GeneCount$Label %in% c("all",site.list.total)),]
DF_GeneCount <- DF_GeneCount[!is.na(DF_GeneCount$iframe),]

# sort(setdiff(c("all",site.list.total),DF_GeneCount$Label))
# DF_GeneCount <- rbind(DF_GeneCount,
#                       data.frame(Plot_Name = "top_gene_count",
#                                  Folder = paste(Folder_root,"/Top_Gene_Count",sep=""),
#                                  Label = sort(setdiff(c("all",site.list.total),DF_GeneCount$Label)),
#                                  iframe = c("2742","2916","3016","3238"),
#                                  stringsAsFactors = FALSE))

# gender_age_distribution_fxn: cohort="all" & site_DF
DF_GenderAge <- mapping_HK[which(mapping_HK$Plot_Name == "gender_age_distribution" &
                                   mapping_HK$Folder == paste(Folder_root,"/Gender_Age_Distribution",sep="")),]
DF_GenderAge <- DF_GenderAge[which(DF_GenderAge$Label %in% c("all",site.list.total)),]
DF_GenderAge <- DF_GenderAge[!is.na(DF_GenderAge$iframe),]

# sort(setdiff(c("all",site.list.total),DF_GenderAge$Label))
# DF_GenderAge <- rbind(DF_GenderAge,
#                       data.frame(Plot_Name = "gender_age_distribution",
#                                  Folder = paste(Folder_root,"/Gender_Age_Distribution",sep=""),
#                                  Label = sort(setdiff(c("all",site.list.total),DF_GenderAge$Label)),
#                                  iframe = c("2744","2918","3018","3240"),
#                                  stringsAsFactors = FALSE))

# pt_mutation_count_fxn: cohort="all" & site_DF
DF_Mutation.all <- mapping_HK[which(mapping_HK$Plot_Name == "pt_mutation_count_allvariants" &
                                      mapping_HK$Folder == paste(Folder_root,"/Patient_Mutation_Count/All_Variants",sep="")),]
DF_Mutation.all <- DF_Mutation.all[which(DF_Mutation.all$Label %in% c("all",site.list.total)),]
DF_Mutation.all <- DF_Mutation.all[!is.na(DF_Mutation.all$iframe),]

# sort(setdiff(c("all",site.list.total),DF_Mutation.all$Label))
# DF_Mutation.all <- rbind(DF_Mutation.all,
#                          data.frame(Plot_Name = "pt_mutation_count_allvariants",
#                                     Folder = paste(Folder_root,"/Patient_Mutation_Count/All_Variants",sep=""),
#                                     Label = sort(setdiff(c("all",site.list.total),DF_Mutation.all$Label)),
#                                     iframe = c("2746","2920","3020","3242"),
#                                     stringsAsFactors = FALSE))

DF_Mutation.stacked <- mapping_HK[which(mapping_HK$Plot_Name == "pt_mutation_count_allvariantstacked" &
                                      mapping_HK$Folder == paste(Folder_root,"/Patient_Mutation_Count/All_Variants_Stacked",sep="")),]
DF_Mutation.stacked <- DF_Mutation.stacked[which(DF_Mutation.stacked$Label %in% c("all",site.list.total)),]
DF_Mutation.stacked <- DF_Mutation.stacked[!is.na(DF_Mutation.stacked$iframe),]

# sort(setdiff(c("all",site.list.total),DF_Mutation.stacked$Label))
# DF_Mutation.stacked <- rbind(DF_Mutation.stacked,
#                              data.frame(Plot_Name = "pt_mutation_count_allvariantstacked",
#                                         Folder = paste(Folder_root,"/Patient_Mutation_Count/All_Variants_Stacked",sep=""),
#                                         Label = sort(setdiff(c("all",site.list.total),DF_Mutation.stacked$Label)),
#                                         iframe = c("5338","5344","5374","5388"),
#                                         stringsAsFactors = FALSE))

DF_Mutation.patho <- mapping_HK[which(mapping_HK$Plot_Name == "pt_mutation_count_pathovariants" &
                                        mapping_HK$Folder == paste(Folder_root,"/Patient_Mutation_Count/Pathogenic_Variants",sep="")),]
DF_Mutation.patho <- DF_Mutation.patho[which(DF_Mutation.patho$Label %in% c("all",site.list.total)),]
DF_Mutation.patho <- DF_Mutation.patho[!is.na(DF_Mutation.patho$iframe),]

# # NA = back, foot, trachea
# sort(setdiff(c("all",site.list.total),DF_Mutation.patho$Label))
# DF_Mutation.patho <- rbind(DF_Mutation.patho,
#                            data.frame(Plot_Name = "pt_mutation_count_pathovariants",
#                                       Folder = paste(Folder_root,"/Patient_Mutation_Count/Pathogenic_Variants",sep=""),
#                                       Label = sort(setdiff(c("all",site.list.total),DF_Mutation.patho$Label)),
#                                       iframe = c(NA,NA,"2922","3022",NA,"3244"),
#                                       stringsAsFactors = FALSE))

DF_Mutation.VUS <- mapping_HK[which(mapping_HK$Plot_Name == "pt_mutation_count_VUS" &
                                      mapping_HK$Folder == paste(Folder_root,"/Patient_Mutation_Count/UnknownSignificance_Variants",sep="")),]
DF_Mutation.VUS <- DF_Mutation.VUS[which(DF_Mutation.VUS$Label %in% c("all",site.list.total)),]
DF_Mutation.VUS <- DF_Mutation.VUS[!is.na(DF_Mutation.VUS$iframe),]

# # NA = mesenteric mass
# sort(setdiff(c("all",site.list.total),DF_Mutation.VUS$Label))
# DF_Mutation.VUS <- rbind(DF_Mutation.VUS,
#                          data.frame(Plot_Name = "pt_mutation_count_VUS",
#                                     Folder = paste(Folder_root,"/Patient_Mutation_Count/UnknownSignificance_Variants",sep=""),
#                                     Label = sort(setdiff(c("all",site.list.total),DF_Mutation.VUS$Label)),
#                                     iframe = c("2748",NA,"2924","3024","3246"),
#                                     stringsAsFactors = FALSE))

# specimen_type_stacked_fxn: cohort="all" & site_DF
DF_SpecimenCount <- mapping_HK[which(mapping_HK$Plot_Name == "specimen_type_count" &
                                       mapping_HK$Folder == paste(Folder_root,"/Specimen_Type_Count",sep="")),]
DF_SpecimenCount <- DF_SpecimenCount[which(DF_SpecimenCount$Label %in% c("all",site.list.total)),]
DF_SpecimenCount <- DF_SpecimenCount[!is.na(DF_SpecimenCount$iframe),]

# sort(setdiff(c("all",site.list.total),DF_SpecimenCount$Label))
# DF_SpecimenCount <- rbind(DF_SpecimenCount,
#                           data.frame(Plot_Name = "specimen_type_count",
#                                      Folder = paste(Folder_root,"/Specimen_Type_Count",sep=""),
#                                      Label = sort(setdiff(c("all",site.list.total),DF_SpecimenCount$Label)),
#                                      iframe = c("2750","2926","3026","3248"),
#                                      stringsAsFactors = FALSE))

# tumor_purity_count_fxn: cohort="all" & site_DF
DF_TumorPurity <- mapping_HK[which(mapping_HK$Plot_Name == "tumor_purity_count" &
                                     mapping_HK$Folder == paste(Folder_root,"/Tumor_Purity_Count",sep="")),]
DF_TumorPurity <- DF_TumorPurity[which(DF_TumorPurity$Label %in% c("all",site.list.total)),]
DF_TumorPurity <- DF_TumorPurity[!is.na(DF_TumorPurity$iframe),]

# sort(setdiff(c("all",site.list.total),DF_TumorPurity$Label))
# DF_TumorPurity <- rbind(DF_TumorPurity,
#                         data.frame(Plot_Name = "tumor_purity_count",
#                                    Folder = paste(Folder_root,"/Tumor_Purity_Count",sep=""),
#                                    Label = sort(setdiff(c("all",site.list.total),DF_TumorPurity$Label)),
#                                    iframe = c("2752","2928","3028","3250"),
#                                    stringsAsFactors = FALSE))

# histologicaldx_distribution_fxn: cohort="all" & site_DF
DF_HistologicalDx <- mapping_HK[which(mapping_HK$Plot_Name == "hist_dx_distribution" &
                                        mapping_HK$Folder == paste(Folder_root,"/Hist_Dx_distribution",sep="")),]
DF_HistologicalDx <- DF_HistologicalDx[which(DF_HistologicalDx$Label %in% c("all",site.list.total)),]
DF_HistologicalDx <- DF_HistologicalDx[!is.na(DF_HistologicalDx$iframe),]

# sort(setdiff(c("all",site.list.total),DF_HistologicalDx$Label))
# DF_HistologicalDx <- rbind(DF_HistologicalDx,
#                            data.frame(Plot_Name = "hist_dx_distribution",
#                                       Folder = paste(Folder_root,"/Hist_Dx_distribution",sep=""),
#                                       Label = sort(setdiff(c("all",site.list.total),DF_HistologicalDx$Label)),
#                                       iframe = c("5342","5348","5378","5392"),
#                                       stringsAsFactors = FALSE))

# top_fusion_count_fxn: cohort="all" & site_DF
fusion.list <- sort(unique(STAMP_Fusion$PrimaryTumorSite))

DF_fusion_count <- mapping_HK[which(mapping_HK$Plot_Name == "top_fusion_count" &
                                        mapping_HK$Folder == paste(Folder_root,"/Top_Fusion_Count",sep="")),]
DF_fusion_count <- DF_fusion_count[which(DF_fusion_count$Label %in% c("all",fusion.list)),]
DF_fusion_count <- DF_fusion_count[!is.na(DF_fusion_count$iframe),]

# sort(setdiff(c("all",fusion.list),DF_fusion_count$Label))
# DF_fusion_count <- rbind(DF_fusion_count,
#                          data.frame(Plot_Name = "top_fusion_count", 
#                               Folder = paste(Folder_root,"/Top_Fusion_Count",sep=""),
#                               Label = sort(setdiff(c("all",fusion.list),DF_fusion_count$Label)),
#                               iframe = c("5372","5380","5382","5384","5386","5394"),
#                               stringsAsFactors = FALSE))
remove(fusion.list)

# test_volume_timeline_fxn: cohort=site_DF
DF_Tumor_TestVolume <- mapping_HK[which(mapping_HK$Plot_Name == "test_volume" &
                                          mapping_HK$Folder == paste(Folder_root,"/PerTumor_OrderVolume",sep="")),]
DF_Tumor_TestVolume <- DF_Tumor_TestVolume[which(DF_Tumor_TestVolume$Label %in% site.list.total),]
DF_Tumor_TestVolume <- DF_Tumor_TestVolume[!is.na(DF_Tumor_TestVolume$iframe),]

# sort(setdiff(site.list.total,DF_Tumor_TestVolume$Label))
# DF_Tumor_TestVolume <- rbind(DF_Tumor_TestVolume,
#                              data.frame(Plot_Name = "test_volume",
#                                         Folder = paste(Folder_root,"/PerTumor_OrderVolume",sep=""),
#                                         Label = sort(setdiff(site.list.total,DF_Tumor_TestVolume$Label)),
#                                         iframe = c("4166","4190","4204","4234"),
#                                         stringsAsFactors = FALSE))

DF_PerSite <- rbind(DF_VariantCount, DF_GeneCount, DF_GenderAge, DF_HistologicalDx, 
                    DF_Tumor_TestVolume,DF_TumorPurity,DF_SpecimenCount,
                    DF_Mutation.all,DF_Mutation.stacked,DF_Mutation.patho,DF_Mutation.VUS,
                    DF_fusion_count)

DF_iframe_FULL <- rbind(DF_misc,DF_PerGene,DF_PerSite)

# Sort columns 
DF_tumor <- DF_iframe_FULL[which(DF_iframe_FULL$Label %in% site.list.total),]
DF_tumor <- DF_tumor[order(DF_tumor$Label, decreasing = FALSE),]

DF_gene <- DF_iframe_FULL[which(DF_iframe_FULL$Label %in% gene.list.total),]
DF_gene <- DF_gene[order(DF_gene$Label, decreasing = FALSE),]

DF_all <- DF_iframe_FULL[which(DF_iframe_FULL$Label == "all"),]
DF_all <- DF_all[order(DF_all$Folder, decreasing = FALSE),]

isTRUE(nrow(DF_tumor) + nrow(DF_gene) + nrow(DF_all) == nrow(DF_iframe_FULL))

DF_iframe_FULL_ordered <- rbind(DF_all,DF_tumor,DF_gene)
remove(DF_tumor,DF_gene,DF_all,DF_iframe_FULL)

DF_iframe_FULL_ordered <- unique(DF_iframe_FULL_ordered[,])

# Remove rows without iframe code
DF_iframe_FULL_ordered <- DF_iframe_FULL_ordered[!(is.na(DF_iframe_FULL_ordered$iframe)),]

## Write to local computer
#----------------------------------------------
write.table(DF_iframe_FULL_ordered, 
            file="~/Documents/ClinicalDataScience_Fellowship/STAMPEDE_Visualizations/mappings.csv",
            append = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
