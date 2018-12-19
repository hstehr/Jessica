## Clean up 20181018_syapse_export_all_variants_patientNameAndMrnRemoved.csv
## Specific parameters that are quality checked
## smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding
## Output: "20181113_syapse_export_all_variants_QC.csv"

rm(list=ls())
setwd("~/Documents/ClinicalDataScience_Fellowship/STAMP/")

# Load relevant file
#----------------------------------------------
DF_Full <- read.csv(file = "20181018_syapse_export_all_variants_patientNameAndMrnRemoved.csv",
                    header = TRUE,
                    na.strings = c(""," ","NA"),
                    stringsAsFactors = FALSE,
                    sep = ",")
# Shorten column names of DF
colnames(DF_Full) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_Full))

# Correction of assayName and smpl.pipelineVersion
#----------------------------------------------
# Variation in spacing character used
DF_Full$smpl.assayName <- gsub("^Myeloid.+", "Myeloid Panel (54 gene panel)", DF_Full$smpl.assayName)

# Fill in missing assayName >> Helio looked at patient report for assayName
#----------------------------------------------
# which(DF_Full$sys.uniqueId == "TRF-513")
# DF_Full[11024:11025,c(1,39,50,53:61)]
DF_Full$smpl.assayName[c(11024,11025)] <- "Cancer Somatic Mutation Panel (48 gene panel)"

# which(DF_Full$smpl.pipelineVersion == "STAMP  v2.2.1")
DF_Full$smpl.pipelineVersion[c(2063:2070)] <- "STAMP v2.2.1"

# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" &
#         DF_Full$smpl.pipelineVersion == "Myeloid v1.0.2")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" &
#         DF_Full$smpl.pipelineVersion == "New")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" &
#         DF_Full$smpl.pipelineVersion == "STAMP v1.2.4" )
DF_Full$smpl.pipelineVersion[c(1893:1895,9234,926:929,1884:1889,2860:2863,6257:6261,9145,9147,
                               9149,9151,9153,9339:9341,9737:9739,7417:7424)] <- "STAMP v2.0.0"
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)" &
#         DF_Full$smpl.pipelineVersion == "STAMP v.2.2.1")
DF_Full$smpl.pipelineVersion[c(6229:6234)] <- "STAMP v2.2.1"
# which(DF_Full$smpl.pipelineVersion == "v2.2.1" | DF_Full$smpl.pipelineVersion == "V2.2.1" |
#         DF_Full$smpl.pipelineVersion == "STAMPv2.2.1" | DF_Full$smpl.pipelineVersion == "STAMP 2.2.1")
DF_Full$smpl.pipelineVersion[c(645:648,346:360,3721:3727,4756:4758,7809,9211:9215,9668,1252:1256,
                               2221:2229,9525:9529,3593,3594,4841:4845,7696:7702)] <- "STAMP v2.2.1"

# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" &
#         DF_Full$smpl.pipelineVersion == "New")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" &
#         DF_Full$smpl.pipelineVersion == "None")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" &
#         DF_Full$smpl.pipelineVersion == "STAMP v2.0")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" &
#         DF_Full$smpl.pipelineVersion == "STAMP v2.0.0")
DF_Full$smpl.pipelineVersion[c(6460:6463,3971,3839,3841,3843,3845,3847,3849,3851,
                               3853,3855,3857,3859,3861,9502:9504,9931)] <- "STAMP v1.0.0"
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)" &
#         DF_Full$smpl.pipelineVersion == "STAMP v.1.2.4")
# which(DF_Full$smpl.pipelineVersion == "STAMP 1.2.4")
DF_Full$smpl.pipelineVersion[c(6508:6515,1015:1036)] <- "STAMP v1.2.4"

# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)" &
#         DF_Full$smpl.pipelineVersion == "User")
# which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)" &
#         DF_Full$smpl.pipelineVersion == "New")
DF_Full$smpl.pipelineVersion[c(517:531,1470:1475,5594:5596,6166,8373,8374,
                               118:120,11640:11643,11952:11954)] <- "STAMP v1.0.0"
# which(DF_Full$smpl.pipelineVersion == "v1.1.5")
DF_Full$smpl.pipelineVersion[c(1260,1261,11650)] <- "STAMP v1.1.5"
# which(DF_Full$smpl.pipelineVersion == "STAMP v.1.1.7")
DF_Full$smpl.pipelineVersion[c(5221:5225,7328,10089:10092)] <- "STAMP v1.1.7"

assay.noversion <- which(DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)" &
                           DF_Full$smpl.pipelineVersion == "None")
DF_Full[assay.noversion, "smpl.pipelineVersion"]  <- "STAMP v1.0.0"
remove(assay.noversion)

# # Printout of smpl.pipelineVersion 
# #----------------------------------------------
# print("pipelineVersion for STAMP (130 genes) are:")
# print(sort(unique(DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (130 genes)",
#                           "smpl.pipelineVersion"])))
# print("pipelineVersion for STAMP (198 genes) are:")
# print(sort(unique(DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)",
#                           "smpl.pipelineVersion"])))
# print("pipelineVersion for STAMP (200 genes) are:")
# print(sort(unique(DF_Full[DF_Full$smpl.assayName == "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)",
#                           "smpl.pipelineVersion"])))

# Modify gene ID to match DF_domain
#----------------------------------------------
# which(DF_Full$base.gene == "TELOMERASE REVERSE TRANSCRIPTASE; TERT")
DF_Full$base.gene[c(12762,12763,12764)] <- "TERT"
# which(DF_Full$base.gene == "1-Mar")
DF_Full$base.gene[c(251,4305,4306)] <- "MARCH1"

# which(DF_Full$smpl.hasAlteration == "TP53 c.747G>T (p.Arg249Ser)")
# DF_Full[c(4763,4864,4865,10270,11024),c(1,39,50,53:61)]
DF_Full[c(11024,10270),"smpl.hgvsGenomic"] <- "g.7577534C>A"
DF_Full[c(11024,10270),"smpl.genomicDescription"] <- "chr17:g.7577534C>A"

# Convert NA cells = "N/A" and/or smpl.hgvsProtein = c("None,"?","p.?","promoter")
DF_Full[DF_Full == "N/A" | DF_Full == "None"] <- NA
# which(DF_Full$smpl.hgvsProtein == "?")
DF_Full$smpl.hgvsProtein[11724] <- NA
# "p.?" indicates protein has not been analysed, an effect is expected but difficult to predict
phgvs.noeffect <- which(DF_Full$smpl.hgvsProtein == "p.?")
phgvs.sys.label <- unique(DF_Full[phgvs.noeffect, "sys.label"])
for (x in 1:length(phgvs.sys.label)) {
  row.change <- which(DF_Full$sys.label == phgvs.sys.label[x])
  DF_Full[row.change, "smpl.hgvsProtein"]  <- NA
}
remove(phgvs.noeffect, phgvs.sys.label, row.change, x)

# Change to same case for consistency
#----------------------------------------------
# which(DF_Full$smpl.hgvsProtein == "PROMOTER")
DF_Full$smpl.hgvsProtein[c(12743,12762,12763,12764)] <- "Promoter"
phgvs.noeffect <- which(DF_Full$smpl.hgvsProtein == "Promoter")
phgvs.sys.label <- unique(DF_Full[phgvs.noeffect, "sys.label"])
for (x in 1:length(phgvs.sys.label)) {
  row.change <- which(DF_Full$sys.label == phgvs.sys.label[x])
  DF_Full[row.change, "smpl.hgvsProtein"]  <- NA
}
remove(phgvs.noeffect, phgvs.sys.label, row.change, x)

# Correct annotation for Cancer Somatic Mutation Panel
#----------------------------------------------
# which(DF_Full$smpl.hgvsCoding == "c.*27T>C")
DF_Full[c(860,1181,1390,2087,2271,2403,3271,5193,6440,8493,
          9106,9122,9702,9795,11655,12220),"smpl.hgvsCoding"] <- "c.X27T>C"

# which(DF_Full$smpl.hgvsCoding == "c.*36A>G")
DF_Full[c(1188,1397),"smpl.hgvsCoding"] <- "c.X36A>G"

# Correct annotation for Myeloid Panel
#----------------------------------------------
# NM_005157.5(ABL1):c.949T>C (p.Phe317Leu)
# which(DF_Full$smpl.hgvsCoding == "p.Phe317Leu")
# DF_Full[12466,c(1,50,53:61)]
DF_Full[12466,"smpl.hgvsCoding"] <- "c.949T>C"

# COSM211736 (c.4636C>T) (p.Q1546*)
# which(DF_Full$smpl.hgvsProtein == "p.Q1546X")
DF_Full[c(1411,1412),"smpl.hgvsProtein"] <- "p.Gln1546Ter"

# Correct annotation for STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
DF_Full$smpl.hgvsProtein <- gsub("NP.+p.", "p.", DF_Full$smpl.hgvsProtein)
DF_Full$smpl.hgvsCoding <- gsub("NM.+c.", "c.", DF_Full$smpl.hgvsCoding)

# Nucleotides downstream (3’) of the translation termination codon (stop) are marked with an asterisk
# Current nomenclature for termination codon: Ter (3-letter) or X (1-letter)
# which(DF_Full$smpl.hgvsProtein == "p.Arg419Stop")
DF_Full$smpl.hgvsProtein[10124] <- "p.Arg419Ter"
# which(DF_Full$smpl.hgvsProtein == "p.Ser166Stop")
DF_Full$smpl.hgvsProtein[10127] <- "p.Ser166Ter"
# which(DF_Full$smpl.hgvsProtein == "p.Met409Stop")
DF_Full$smpl.hgvsProtein[10128] <- "p.Met409Ter"

# which(DF_Full$smpl.hgvsCoding == "c.586_*3delCGTCAAACGTAAACA")
DF_Full$smpl.hgvsCoding[4607] <- "c.586_X3delCGTCAAACGTAAACA"
# which(DF_Full$smpl.hgvsCoding == "c.5424_*1dup")
DF_Full$smpl.hgvsCoding[7494] <- "c.5424_X1dup"
# which(DF_Full$smpl.hgvsCoding == "c.457-18_*5del")
DF_Full$smpl.hgvsCoding[c(2042,2043)] <- "c.457-18_X5del"
# which(DF_Full$smpl.hgvsCoding == "c.1328_*33del")
DF_Full$smpl.hgvsCoding[5587] <- "c.1328_X33del"
# which(DF_Full$smpl.hgvsCoding == "c.1180_*190del")
DF_Full$smpl.hgvsCoding[10307] <- "c.1180_X190del"
# which(DF_Full$smpl.hgvsCoding == "c.*85_*87delGAA")
DF_Full$smpl.hgvsCoding[c(617,9203)] <- "c.X85_X87delGAA"

# which(DF_Full$smpl.transcript == "NM_000077.4:c.457-18_*5del")
DF_Full[2042,"smpl.transcript"] <- "NM_000077.4"

# chr17:g.7578410T>A = NM_000546.5(TP53):c.520A>G (p.Arg174Gly)
# which(DF_Full$smpl.hgvsProtein == "p.Arg174Trp")
# DF_Full[c(3682,5656),c(1,50,53:61)]
DF_Full$smpl.hgvsCoding[c(3682,5656)] <- "c.520A>G"

# which(DF_Full$sys.label == "CSMD3 p.Trp1403Ter (c.4209G>A)")
# DF_Full[3683,c(1,50,53:61)]
DF_Full$smpl.hgvsCoding[3683] <- "c.4209G>A"
DF_Full$smpl.hgvsProtein[3683] <- "p.Trp1403Ter"

# which(DF_Full$smpl.hgvsProtein == "Glu746_Ala750del")
DF_Full$smpl.hgvsProtein[7212] <- "p.Glu746_Ala750del"
# which(DF_Full$smpl.hgvsProtein == "p. Thr574_Asn587dup")
DF_Full$smpl.hgvsProtein[c(8104,8105)] <- "p.Thr574_Asn587dup"
# which(DF_Full$smpl.hgvsProtein == ":p.Phe221fs")
DF_Full$smpl.hgvsProtein[6200] <- "p.Phe221fs"
# which(DF_Full$smpl.hgvsProtein == "p.A129fs")
DF_Full$smpl.hgvsProtein[4082] <- "p.Ala129fs"
# which(DF_Full$smpl.hgvsProtein == " p.Ile251_Leu252dup")
DF_Full$smpl.hgvsProtein[7080] <- "p.Ile251_Leu252dup"
# which(DF_Full$smpl.hgvsProtein == "p.L747_T751del")
DF_Full$smpl.hgvsProtein[c(10281,10457)] <- "p.Leu747_Thr751del"
# which(DF_Full$smpl.hgvsProtein == "p.Ser3366Asnf")
DF_Full$smpl.hgvsProtein[c(4562,4855,5845)] <- "p.Ser3366Asnfs"
# which(DF_Full$smpl.hgvsCoding == "c.586_X3delCGTCAAACGTAAACA")
DF_Full$smpl.hgvsProtein[4607] <- "p.Arg196_Thr198delinsAlaArgIleLysAsnMetPheProCysLeuSerAspThrSerLeuLeuAspGluAlaArgLysIleTyrMetLysIleLeuLysIleHisIleAlaAspPheMetGluTrpThrSerCysIleSerThrGluLysGlnGlnHisAsnAsnThrLysIleLeuGlyThrLeuLys"

# Mutalyzer NM_005359.5(SMAD4_v001):c.1278_1292delinsGCATATATA = p.Lys428_Asp552delinsIle
# which(DF_Full$smpl.hgvsProtein == "p.K428_I429delinsIX")
DF_Full$smpl.hgvsProtein[9636] <- "p.Lys428_Asp552delinsIle"

# https://www.ncbi.nlm.nih.gov/snp/rs397516975#variant_details
# which(DF_Full$smpl.hgvsProtein == "p.Y772_A775dup")
# DF_Full[c(3678,5023,5475,5476,11511),c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[c(3678,5023,5475,5476,11511)] <- "p.Tyr772_Ala775dup"
DF_Full$smpl.hgvsGenomic[3678] <- "g.37880984_37880995dup"
DF_Full$smpl.genomicDescription[3678] <- "chr17:g.37880984_37880995dup"

# NM_000548.4(TSC2):c.3G>A (p.Met1Ile) -- Met(AUG) >> Ile(AUA, AUU)
# DF_Full[c(1386,4412,4418,1734,10044),c(1,50,53:61)]
# which(DF_Full$sys.label == "TSC2 c.3G>T (p.?)" | DF_Full$sys.label == "MYCN c.3G>T (p.?)" |
#         DF_Full$sys.label == "CDKN1B c.3G>A (p.?)" | DF_Full$sys.label == "BAP1 c.3G>A (p.?)")
DF_Full$smpl.hgvsProtein[c(1636,2063,5334,5340,12372)] <- "p.Met1Ile"

# c.2T>C -- Met(AUG) >> Thr(ACG)
# which(DF_Full$smpl.hgvsProtein == "p.M1T")
DF_Full$smpl.hgvsProtein[c(1076,12937)] <- "p.Met1Thr"

# Mutalyzer: NM_001042705.2(IQCJ_i001):p.(Arg152Gly)
# which(DF_Full$smpl.hgvsCoding == "c.454A>G")
# DF_Full[2770,c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[2770] <- "p.Arg152Gly"

# Translation termination codon (stop codon) is indicated by "Ter" or "*"
# which(DF_Full$smpl.hgvsProtein == "p.*394Gluext*4")
# DF_Full[10307,c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[10307] <- "p.Ter394GluextTer4"
DF_Full$smpl.genomicDescription[10307] <- "chr17:g.7572926_7572929delGTCA"

# which(DF_Full$smpl.hgvsProtein == "p.Pro408Alafs*99")
DF_Full$smpl.hgvsProtein[12365] <- "p.Pro408AlafsTer99"
# which(DF_Full$smpl.hgvsProtein == "p.Ser437Lysfs*70")
DF_Full$smpl.hgvsProtein[4500] <- "p.Ser437LysfsTer70"
# which(DF_Full$smpl.hgvsProtein == "p.Ter444Valext*49")
DF_Full$smpl.hgvsProtein[5587] <- "p.Ter444ValextTer49"

# gnomAD: 12:25362768 CTTT / C == p.Lys180del
# Actual correct annotation [Mutalyzer]: NM_033360.2:c.*80_*82delAAG
# which(DF_Full$sys.label == "KRAS c.*85_*87delGAA (N/A)")
# DF_Full[c(617,9203),c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[c(617,9203)] <- "p.Lys180del"

# NM_058195.3(CDKN2A):c.92C>T (p.Thr31Met) -- Thr(ACG) >> Met(AUG) >> position 92 is aa 2
# c.92C>G == Thr(ACG) >> Arg(AGG)
# which(DF_Full$sys.label == "CDKN2A c.92C>G (N/A)")
# DF_Full[c(3186,5775),c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[c(3186,5775)] <- "p.Thr31Arg"

# NM_001195132.1(CDKN2A):c.496C>T (p.His166Tyr) AND not specified
# which(DF_Full$sys.label == "CDKN2A c.496C>T (N/A)")
# DF_Full[c(5489,6166),c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[c(5489,6166)] <- "p.His166Tyr"

# gnomAD: 2:80529645 C / T (rs756853344) == p.Ala434Thr
# which(DF_Full$sys.label == "CTNNA2 c.1300G>A (N/A)")
# DF_Full[c(701,702,6290),c(1,50,53:61)]
DF_Full$smpl.hgvsProtein[c(701,702,6290)] <- "p.Ala434Thr"

# Backcalculation and similarity approximation
#----------------------------------------------
# Upstream "smpl.hgvsCoding" coordinates assumes no breaks in non-coding region

# TERT (chr5:g.1295105C>T c.-1G>A) & (chr5:g.1295322G>C c.-218C>G)
# which(DF_Full$smpl.genomicDescription == "chr5:g.1295276C>T")
# DF_Full[3635,c(1,50,53:61)]
DF_Full[3635,"smpl.hgvsCoding"] <- "c.-172G>A"
# which(DF_Full$smpl.genomicDescription == "chr5:g.1295315G>T")
# DF_Full[1637,c(1,50,53:61)]
DF_Full[1637,"smpl.hgvsCoding"] <- "c.-211C>A"
# which(DF_Full$smpl.genomicDescription == "chr5:g.1295295G>A")
# DF_Full[6469,c(1,50,53:61)]
DF_Full[6469,"smpl.hgvsCoding"] <- "c.-191C>T"
# which(DF_Full$smpl.genomicDescription == "chr5:g.1295320C>T")
# DF_Full[4126,c(1,50,53:61)]
DF_Full[4126,"base.gene"] <- "TERT"
DF_Full[4126,"smpl.hgvsCoding"] <- "c.-216G>A"
DF_Full[4126,"smpl.transcript"] <- "NM_198253.2"
# which(DF_Full$smpl.genomicDescription == "chr5:g.1295258G>A")
# DF_Full[8129,c(1,50,53:61)]
DF_Full[8129,"base.gene"] <- "TERT"
DF_Full[8129,"smpl.hgvsCoding"] <- "c.-154C>T"
DF_Full[8129,"smpl.transcript"] <- "NM_198253.2"

## Remove duplicates of entries based on multiple unique sys.date_created.1
## Include cases of smpl.amendedString == "AMENDED"
#----------------------------------------------
DF_Full_AMENDED <- data.frame()
patient_list <- unique(DF_Full$sys.uniqueId)

for (id_num in 1:length(patient_list)) {
  patient_id <- patient_list[id_num]
  
  assay_list <- unique(DF_Full$smpl.assayName[which(DF_Full$sys.uniqueId == patient_id)])
  
  for (assay_num in 1:length(assay_list)) {
    assay_id <- assay_list[assay_num]
    
    if (length(unique(DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id), 11])) > 1) {
      created_list <- unique(DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id), 11])
      
      ## Retain most recent version if multiple sys.date_created.1 exist in report 
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id),]
      DF_Full_AMENDED_pre <- DF_Full_AMENDED_pre[DF_Full_AMENDED_pre[[11]] == max(created_list),]
      
    } else {
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id &
                                             DF_Full$smpl.assayName == assay_id),]
    }
    
    DF_Full_AMENDED <- rbind(DF_Full_AMENDED, DF_Full_AMENDED_pre)
  }
}

DF_Full <- DF_Full_AMENDED

remove(DF_Full_AMENDED, DF_Full_AMENDED_pre,assay_list,assay_id,assay_num,
       patient_list,patient_id,id_num,created_list)

# ## Simple review of data
# #----------------------------------------------
# print(paste("Total number of unique patient_id: ", length(unique(DF_Full$sys.uniqueId)), sep=""))
# print(paste("Total number of unique specimen sites: ", length(sort(unique(DF_Full$smpl.specimenSite))), sep=""))
# sort(unique(DF_Full$smpl.specimenSite))
# print(paste("Total number of unique Dx: ", length(sort(unique(DF_Full$smpl.ppDiagnosticSummary))), sep=""))
# sort(unique(DF_Full$smpl.ppDiagnosticSummary))
# print(paste("Total number of unique genes: ", length(sort(unique(DF_Full$base.gene))), sep=""))
# sort(unique(DF_Full$base.gene))
# print(paste("Total number of unique pathogenicity statuses: ", length(sort(unique(DF_Full$smpl.pathogenicityStatus))), sep=""))
# table(sort(DF_Full$smpl.pathogenicityStatus))

## Write to local computer
#----------------------------------------------
write.csv(DF_Full, file = "20181113_syapse_export_all_variants_QC.csv",
          na = "NA", row.names = FALSE)
