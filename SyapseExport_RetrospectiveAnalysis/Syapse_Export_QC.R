## QC-parameters: smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding, smpl.hgvsGenomic
## Calculate patient age > merge with primary tumor site data > filter for STAMP entries 
## Output: "syapse_export_all_variants_QC.csv"

cat(paste("Timestamp of Syapse_Report: ", Syapse_Export_timestamp, sep=""),"\n","\n")

# Load relevant file from global environment 
#----------------------------------------------
DF_Full <- STAMP_DF

# Shorten column names of DF
colnames(DF_Full) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_Full))

# Subset columns of interest
colnames_extract <- c("sys.uniqueId","base.dob","smpl.gender",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","smpl.hasOrderingPhysician",
                      "sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","base.gene","smpl.pathogenicityStatus",
                      "smpl.hgvsCoding","smpl.hgvsProtein","smpl.transcript",
                      "smpl.genomicDescription","base.chromosome","smpl.hgvsGenomic","smpl.chromosomePositionStart")

DF_Full <- DF_Full[,colnames_extract]

# General structurization 
#----------------------------------------------
# Convert NA cells = "N/A" and/or smpl.hgvsProtein = c("None,"?","p.?","promoter")
DF_Full[DF_Full == "N/A" | DF_Full == "None" | DF_Full == "?"] <- NA

# "p.?" indicates protein has not been analysed, an effect is expected but difficult to predict
phgvs.noeffect <- which(DF_Full$smpl.hgvsProtein == "p.?")
phgvs.sys.label <- unique(DF_Full[phgvs.noeffect, "sys.label"])
for (x in 1:length(phgvs.sys.label)) {
  row.change <- which(DF_Full$sys.label == phgvs.sys.label[x])
  DF_Full[row.change, "smpl.hgvsProtein"]  <- NA
}

# HGVS nomenclature
DF_Full$smpl.hgvsProtein <- gsub("NP.+p.", "p.", DF_Full$smpl.hgvsProtein)
DF_Full$smpl.hgvsCoding <- gsub("NM.+c.", "c.", DF_Full$smpl.hgvsCoding)

# Nucleotides downstream (3’) of the translation termination codon (stop) are marked with an asterisk
# Use updated nomenclature for termination codon: Ter (3-letter) or X (1-letter)
DF_Full$smpl.hgvsProtein <- gsub("Stop$", "Ter", DF_Full$smpl.hgvsProtein)
DF_Full$smpl.hgvsProtein <- gsub("\\*", "Ter", DF_Full$smpl.hgvsProtein)

DF_Full$smpl.hgvsCoding <- gsub("\\*", "X", DF_Full$smpl.hgvsCoding)

# Remove trailing whitespace 
DF_Full$smpl.hgvsProtein <- gsub("^[[:space:]]*p.[[:space:]]*", "p.", DF_Full$smpl.hgvsProtein)
DF_Full$smpl.hgvsProtein <- gsub("[[:space:]]+$", "", DF_Full$smpl.hgvsProtein)

DF_Full$smpl.hgvsCoding <- gsub("^[[:space:]]*c.[[:space:]]*", "c.", DF_Full$smpl.hgvsCoding)
DF_Full$smpl.hgvsCoding <- gsub("[[:space:]]+$", "", DF_Full$smpl.hgvsCoding)

###################################
## Manual edit BEGIN
###################################
if (isTRUE(Syapse_Export_timestamp  == "2019-04-30")) {
  
  # Modify gene ID to match DF_domain
  #----------------------------------------------
  row.change = which(DF_Full$base.gene == "TELOMERASE REVERSE TRANSCRIPTASE; TERT")
  DF_Full$base.gene[row.change] <- "TERT"
  
  row.change = which(DF_Full$base.gene == "1-Mar")
  DF_Full$base.gene[row.change] <- "MARCH1"
  
  row.change = which(DF_Full$sys.label == "TP53 c.747G>T (p.Arg249Ser)")
  DF_Full$smpl.hgvsGenomic[row.change] <- "g.7577534C>A"
  DF_Full$smpl.genomicDescription[row.change] <- "chr17:g.7577534C>A"
  
  # Change to same capitalization scheme for consistency
  #----------------------------------------------
  row.change = which(tolower(DF_Full$smpl.hgvsProtein) == "promoter")
  DF_Full$smpl.hgvsProtein[row.change] <- "Promoter"
  phgvs.noeffect <- which(DF_Full$smpl.hgvsProtein == "Promoter")
  phgvs.sys.label <- unique(DF_Full[phgvs.noeffect, "sys.label"])
  for (x in 1:length(phgvs.sys.label)) {
    row.change <- which(DF_Full$sys.label == phgvs.sys.label[x])
    DF_Full$smpl.hgvsProtein[row.change]  <- NA
  }
  
  # Correct annotation 
  #----------------------------------------------
  row.change = which(DF_Full$smpl.hgvsCoding == ":c.1786_1806dup")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.1786_1806dup"
  
  row.change = which(DF_Full$smpl.hgvsProtein == "Met9fs")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Met9fs"
  
  # NM_005157.5(ABL1):c.949T>C (p.Phe317Leu)
  row.change = which(DF_Full$smpl.hgvsCoding == "p.Phe317Leu")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.949T>C"
  
  # COSM211736 (c.4636C>T) (p.Q1546*)
  row.change = which(DF_Full$smpl.hgvsProtein == "p.Q1546X")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Gln1546Ter"
  
  row.change = which(DF_Full$smpl.hgvsCoding == "c.457-18_X5del")
  DF_Full$smpl.transcript[row.change] <- "NM_000077.4"
  
  # chr17:g.7578410T>A = NM_000546.5(TP53):c.520A>G (p.Arg174Gly)
  row.change = which(DF_Full$smpl.hgvsProtein == "p.Arg174Trp")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.520A>G"
  
  row.change = which(DF_Full$sys.label == "CSMD3 p.Trp1403Ter (c.4209G>A)")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.4209G>A"
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Trp1403Ter"
  
  row.change = which(DF_Full$smpl.hgvsProtein == "Glu746_Ala750del")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Glu746_Ala750del"
  
  row.change = which(DF_Full$smpl.hgvsProtein == ":p.Phe221fs")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Phe221fs"
  
  row.change = which(DF_Full$smpl.hgvsProtein == "p.A129fs")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Ala129fs"
  
  row.change = which(DF_Full$smpl.hgvsProtein == "p.L747_T751del")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Leu747_Thr751del"
  
  row.change = which(DF_Full$smpl.hgvsProtein == "p.Ser3366Asnf")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Ser3366Asnfs"
  
  row.change = which(DF_Full$smpl.hgvsCoding == "c.586_X3delCGTCAAACGTAAACA")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Arg196_Thr198delinsAlaArgIleLysAsnMetPheProCysLeuSerAspThrSerLeuLeuAspGluAlaArgLysIleTyrMetLysIleLeuLysIleHisIleAlaAspPheMetGluTrpThrSerCysIleSerThrGluLysGlnGlnHisAsnAsnThrLysIleLeuGlyThrLeuLys"
  
  # Mutalyzer NM_005359.5(SMAD4_v001):c.1278_1292delinsGCATATATA = p.Lys428_Asp552delinsIle
  row.change = which(DF_Full$smpl.hgvsProtein == "p.K428_I429delinsIX")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Lys428_Asp552delinsIle"
  
  # https://www.ncbi.nlm.nih.gov/snp/rs397516975#variant_details
  row.change = which(DF_Full$smpl.hgvsProtein == "p.Y772_A775dup")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Tyr772_Ala775dup"
  DF_Full$smpl.hgvsGenomic[row.change] <- "g.37880984_37880995dup"
  DF_Full$smpl.genomicDescription[row.change] <- "chr17:g.37880984_37880995dup"
  
  # NM_000548.4(TSC2):c.3G>A (p.Met1Ile) -- Met(AUG) >> Ile(AUA, AUU)
  # DF_Full[c(1386,4412,4418,1734,10044),c(1,50,53:61)]
  row.change = which(DF_Full$sys.label == "TSC2 c.3G>T (p.?)" | DF_Full$sys.label == "MYCN c.3G>T (p.?)" |
                       DF_Full$sys.label == "CDKN1B c.3G>A (p.?)" | DF_Full$sys.label == "BAP1 c.3G>A (p.?)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Met1Ile"
  
  # c.2T>C -- Met(AUG) >> Thr(ACG)
  row.change = which(DF_Full$smpl.hgvsProtein == "p.M1T")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Met1Thr"
  
  # Mutalyzer: NM_001042705.2(IQCJ_i001):p.(Arg152Gly)
  row.change = which(DF_Full$smpl.hgvsCoding == "c.454A>G")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Arg152Gly"
  
  row.change = which(DF_Full$sys.label == "TP53 c.1180_*190del (p.*394Gluext*4)")
  DF_Full$smpl.genomicDescription[row.change] <- "chr17:g.7572926_7572929delGTCA"
  
  # gnomAD: 12:25362768 CTTT / C == p.Lys180del
  # Actual correct annotation [Mutalyzer]: NM_033360.2:c.*80_*82delAAG
  row.change = which(DF_Full$sys.label == "KRAS c.*85_*87delGAA (N/A)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Lys180del"
  
  # NM_058195.3(CDKN2A):c.92C>T (p.Thr31Met) -- Thr(ACG) >> Met(AUG) >> position 92 is aa 2
  # c.92C>G == Thr(ACG) >> Arg(AGG)
  row.change = which(DF_Full$sys.label == "CDKN2A c.92C>G (N/A)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Thr31Arg"
  
  # NM_001195132.1(CDKN2A):c.496C>T (p.His166Tyr) AND not specified
  row.change = which(DF_Full$sys.label == "CDKN2A c.496C>T (N/A)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.His166Tyr"
  
  # gnomAD: 2:80529645 C / T (rs756853344) == p.Ala434Thr
  row.change = which(DF_Full$sys.label == "CTNNA2 c.1300G>A (N/A)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Ala434Thr"
  
  # https://www.ncbi.nlm.nih.gov/clinvar/variation/133881/
  row.change = which(DF_Full$sys.label == "CDKN2A c.458-492G>C (N/A)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Arg165Ser"
  
  # https://www.ncbi.nlm.nih.gov/clinvar/variation/90833/ 
  row.change = which(DF_Full$sys.label == "MSH2 c.1A>G (p.?)")
  DF_Full$smpl.hgvsProtein[row.change] <- "p.Met1Val"
  
  # Backcalculation and similarity approximation
  #----------------------------------------------
  # Upstream "smpl.hgvsCoding" coordinates assumes no breaks in non-coding region
  
  # TERT (chr5:g.1295105C>T c.-1G>A) & (chr5:g.1295322G>C c.-218C>G)
  row.change = which(DF_Full$smpl.genomicDescription == "chr5:g.1295276C>T")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.-172G>A"
  
  row.change = which(DF_Full$smpl.genomicDescription == "chr5:g.1295315G>T")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.-211C>A"
  
  row.change = which(DF_Full$smpl.genomicDescription == "chr5:g.1295295G>A")
  DF_Full$smpl.hgvsCoding[row.change] <- "c.-191C>T"
  
  row.change = which(DF_Full$smpl.genomicDescription == "chr5:g.1295320C>T")
  DF_Full$base.gene[row.change] <- "TERT"
  DF_Full$smpl.hgvsCoding[row.change] <- "c.-216G>A"
  DF_Full$smpl.transcript[row.change] <- "NM_198253.2"
  
  row.change = which(DF_Full$smpl.genomicDescription == "chr5:g.1295258G>A")
  DF_Full$base.gene[row.change] <- "TERT"
  DF_Full$smpl.hgvsCoding[row.change] <- "c.-154C>T"
  DF_Full$smpl.transcript[row.change] <- "NM_198253.2"
  
  row.change = which(DF_Full$sys.label == "BRAF c.1798_1799delinsAA (p.Val600Lys)")
  DF_Full$smpl.genomicDescription[row.change] <- "chr7:g.140453136_140453137delACinsTT"
  DF_Full$base.chromosome[row.change]  <- "chr7"
  
  row.change = which(DF_Full$sys.label == "DNMT3A c.2645G>C (p.Arg882Pro)")
  DF_Full$base.chromosome[row.change]  <- "chr2"
  DF_Full$smpl.genomicDescription[row.change] <- "chr2:g.25457242C>G"
  
  row.change = which(DF_Full$sys.label == "IDH2 c.419G>T (p.Arg140Leu)")
  DF_Full$base.chromosome[row.change]  <- "chr15"
  DF_Full$smpl.genomicDescription[row.change] <- "chr15:g.90631934C>A"
  
  row.change = which(DF_Full$sys.label == "TP53 c.265C>G (p.Pro89Ala)")
  DF_Full$base.chromosome[row.change]  <- "chr17"
  DF_Full$smpl.genomicDescription[row.change] <- "chr17:g.7579422G>C"
  
  row.change = which(DF_Full$sys.label == "VHL c.472C>G (p.Leu158Val)")
  DF_Full$smpl.genomicDescription[row.change] <- "chr3:g.10191479C>G"
  
  row.change = which(DF_Full$sys.label == "TET2 c.1648C>T (p.Arg550Ter)")
  DF_Full$smpl.transcript[row.change]  <- "NM_001127208.2"
  DF_Full$base.chromosome[row.change]  <- "chr4"
  DF_Full$smpl.genomicDescription[row.change] <- "chr4:g.106156747C>T"
  
  row.change = which(DF_Full$sys.label == "TET2 c.2645_2646delinsAA (p.Cys882Ter)")
  DF_Full$base.chromosome[row.change]  <- "chr4"
  DF_Full$smpl.genomicDescription[row.change] <- "chr4:g.106157744_106157745delinsAA"
  
  row.change = which(DF_Full$sys.label == "EGFR c.89-7736_889+41del (p.Leu30_Arg297delinsGly)")
  DF_Full$base.chromosome[row.change]  <- "chr7"
  DF_Full$smpl.genomicDescription[row.change] <- "chr7:g.55202243_55221886del"
  
  row.change = which(DF_Full$sys.label == "IDH2 c.419G>A (p.Arg140Gln)")
  DF_Full$base.chromosome[row.change]  <- "chr15"
  DF_Full$smpl.genomicDescription[row.change] <- "chr15:g.90631934C>T"
  
  row.change = which(DF_Full$sys.label == "PIK3CA c.1348_1374del (p.His450_Pro458del)")
  DF_Full$base.chromosome[row.change]  <- "chr3"
  DF_Full$smpl.genomicDescription[row.change] <- "chr3:g.1505_1531del"
  
  row.change = which(DF_Full$sys.label == "PTCH1 c.3947A>G (p.Tyr1316Cys)")
  DF_Full$base.chromosome[row.change]  <- "chr9"
  DF_Full$smpl.genomicDescription[row.change] <- "chr9:g.98209591T>C"
  
  row.change = which(DF_Full$sys.label == "TERT c.-138_139delinsTT ()")
  DF_Full$base.chromosome[row.change]  <- "chr5"
  DF_Full$smpl.genomicDescription[row.change] <- "chr5:g.1295242_1295243delinsAA"
  
  row.change = which(DF_Full$sys.label == "TET2 c.1072A>G (p.Ser358Gly)")
  DF_Full$base.chromosome[row.change]  <- "chr4"
  DF_Full$smpl.genomicDescription[row.change] <- "chr4:g.106156171A>G"

  row.change = which(DF_Full$sys.label == "TP53 NM_000546.5:c.856G>A (NP_000537.3:p.Glu286Lys)")
  DF_Full$smpl.genomicDescription[row.change] <- "chr17:g.7577082C>T"
  DF_Full$smpl.hgvsGenomic[row.change] <- "g.7577082C>T"
  
  # https://mutalyzer.nl/name-checker?description=NM_000077.4%3Ac.26delinsTT+
  row.change = which(DF_Full$sys.label == "CDKN2A NM_000077.4:c.26delinsTT (Met9fs)")
  DF_Full$smpl.pathogenicityStatus[row.change] <- "Pathogenic"
  DF_Full$smpl.hgvsCoding[row.change] <- "c.26dupT"
  DF_Full$smpl.genomicDescription[row.change] <- "chr9:g.21974801dupA"
  DF_Full$smpl.hgvsGenomic[row.change] <- "g.21974801dupA"
  
  row.change = which(DF_Full$sys.label == "KIT c.1612A>G (p.Ile538Val)")
  DF_Full$base.chromosome[row.change] <- "chr4"
  DF_Full$smpl.genomicDescription[row.change] <- "chr4:g.55593455A>G"
  
  row.change = which(DF_Full$sys.label == "KIT c.1594G>A (p.Val532Ile)")
  DF_Full$smpl.genomicDescription[row.change] <- "chr4:g.55593437G>A"
  
  row.change = which(DF_Full$sys.label == "TERT c.-146C>T (N/A)")
  DF_Full$smpl.pathogenicityStatus[row.change] <- "Pathogenic"
  
  row.change = which(DF_Full$sys.label == "PIK3CA c.1348_1374del (p.His450_Pro458del)")
  DF_Full$smpl.genomicDescription[row.change] <- NA
  DF_Full$smpl.hgvsGenomic[row.change] <- NA
  
  row.change = which(DF_Full$sys.label == "RB1 c.1592_1634delinsGCTTTT (p.Val531fs)")
  DF_Full$smpl.genomicDescription[row.change] <- NA
  DF_Full$smpl.hgvsGenomic[row.change] <- NA
  
  row.change = which(grepl("^chr[[:digit:]]+:",DF_Full$smpl.hgvsGenomic) == TRUE) 
  DF_Full$smpl.hgvsGenomic[row.change] <- gsub("^chr[[:digit:]]+:","",DF_Full$smpl.hgvsGenomic[row.change])

  # VariantHGVSGenomic correction
  #----------------------------------------------
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "PTCH1 c.3947A>G (p.Tyr1316Cys)")] <-"g.98209591T>C"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "KIT c.1721_1762dup (p. Thr574_Asn587dup)")] <- "g.55593655_55593696dup"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "TERT c.-138_139delinsTT ()")] <- "g.1295242_1295243delinsAA"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "EGFR c.89-7736_889+41del (p.Leu30_Arg297delinsGly)")] <- "g.55202243_55221886del"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "MET c.2942-18_2953del (p.?)")] <- "g.116411885_116411914del"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "DNMT3A c.2645G>C (p.Arg882Pro)")] <- "g.25457242C>G"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "TET2 c.1072A>G (p.Ser358Gly)")] <- "g.106156171A>G"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "IDH2 c.419G>T (p.Arg140Leu)")] <- "g.90631934C>A"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$sys.label == "TET2 c.2645_2646delinsAA (p.Cys882Ter)")] <- "g.106157744_106157745delinsAA"
  
  DF_Full$smpl.hgvsGenomic[which(DF_Full$smpl.hgvsGenomic == ":g.10602467 _10602472dup")] <- "g.10602467_10602472dup"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$smpl.hgvsGenomic == "g.10191479")] <- "g.10191479C>G"
  DF_Full$smpl.hgvsGenomic[which(DF_Full$smpl.hgvsGenomic == "NC_000016.9:g.68842671_68842677dup")] <- "g.68842671_68842677dup"
  
  row.change = which(is.na(DF_Full$smpl.genomicDescription))
  for (row_No in 1:length(row.change)) {
    row.change_sub <- row.change[row_No]
    
    chr_No = DF_Full$base.chromosome[row.change_sub]
    gen_No = DF_Full$smpl.hgvsGenomic[row.change_sub]
    
    if (isTRUE(!is.na(chr_No) & !is.na(gen_No))) {
      DF_Full$smpl.genomicDescription[row.change_sub] <- paste(chr_No,":",gen_No,sep="")
    }
  }
  
  sort(unique(DF_Full$smpl.genomicDescription[row.change]))

  DF_Full$smpl.histologicalDiagnosis[which(DF_Full$smpl.histologicalDiagnosis == "Sqaumous cell carcinoma")] <- "Squamous cell carcinoma"
  DF_Full$smpl.histologicalDiagnosis[which(DF_Full$smpl.histologicalDiagnosis == "Carcinoma:Other (specify)")] <- "Carcinoma:Other carcinoma (specify)"
  
  DF_Full$smpl.primaryTumorSite[which(DF_Full$smpl.primaryTumorSite == "Other Primary Site (Including Pediatric)")] <- "Other Primary Site"
  
  DF_Full$smpl.pathogenicityStatus[which(DF_Full$sys.label == "TERT c.-146C>T (N/A)")] <- "Pathogenic"
  DF_Full$smpl.pathogenicityStatus[which(DF_Full$sys.label == "CDKN2A NM_000077.4:c.26delinsTT (Met9fs)")] <- "Pathogenic"
}

# ## Compare to new terms added since mapping
# dx_new <- data.frame(HistologicalDx=sort(unique(tolower(DF_Full$smpl.histologicalDiagnosis))), stringsAsFactors = FALSE)
# dx_mapped <- read.csv(file = "~/Desktop/HistologicalDx.tsv", header = TRUE, stringsAsFactors = FALSE,sep = ",")
# dx_mapped$HistologicalDx <- tolower(dx_mapped$HistologicalDx)
# dx_missing <- rbind(anti_join(dx_mapped,dx_new),anti_join(dx_new,dx_mapped))
# 
# dx_new <- data.frame(PrimaryTumorSite=sort(unique(tolower(DF_Full$smpl.primaryTumorSite))), stringsAsFactors = FALSE)
# dx_mapped <- read.csv(file = "~/Desktop/PrimaryTumorSite.tsv", header = TRUE, stringsAsFactors = FALSE,sep = ",")
# dx_mapped$PrimaryTumorSite <- tolower(dx_mapped$PrimaryTumorSite)
# dx_missing <- rbind(anti_join(dx_mapped,dx_new),anti_join(dx_new,dx_mapped))
# remove(dx_new,dx_mapped,dx_missing)
# 
# ## Check columns for obvious errors 
# sort(unique(DF_Full$smpl.gender))
# sort(unique(DF_Full$smpl.specimenType))
# sort(unique(DF_Full$smpl.percentTumor))
# sort(unique(DF_Full$smpl.assayName))
# sort(unique(DF_Full$base.gene))
# sort(unique(DF_Full$smpl.pathogenicityStatus))
# data.frame(coding=sort(unique(DF_Full$smpl.hgvsCoding)),stringsAsFactors = FALSE) -> A
# data.frame(protein=sort(unique(DF_Full$smpl.hgvsProtein)),stringsAsFactors = FALSE) -> A
# for (row_No in 1:nrow(A)) {
#   if (isTRUE(grepl("^p.[[:alpha:]]{3}",A$protein[row_No]) == "FALSE")) {print(A$protein[row_No])
#   }
# }
# sort(unique(DF_Full$base.chromosome))
# remove(A,row_No,row.change)

###################################
## Manual edit END
###################################

## Remove duplicates of entries based on multiple unique sys.date_created
## Include cases of smpl.amendedString == "AMENDED"
#----------------------------------------------
DF_Full_AMENDED <- data.frame()
patient_list <- unique(DF_Full$sys.uniqueId)

for (id_num in 1:length(patient_list)) {
  patient_id <- patient_list[id_num]
  
  assay_list <- unique(DF_Full$smpl.assayName[which(DF_Full$sys.uniqueId == patient_id)])
  
  for (assay_num in 1:length(assay_list)) {
    assay_id <- assay_list[assay_num]
    
    if (length(unique(DF_Full$sys.date_created[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id)])) > 1) {
      created_list <- unique(DF_Full$sys.date_created[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id)])
      
      ## Retain most recent version if multiple sys.date_created.1 exist in report 
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id),]
      DF_Full_AMENDED_pre <- DF_Full_AMENDED_pre[DF_Full_AMENDED_pre$sys.date_created == max(created_list),]
      
      remove(created_list)
      
    } else {
      DF_Full_AMENDED_pre <- DF_Full[which(DF_Full$sys.uniqueId == patient_id & DF_Full$smpl.assayName == assay_id),]
    }
    
    DF_Full_AMENDED <- rbind(DF_Full_AMENDED, DF_Full_AMENDED_pre)
    remove(assay_id,assay_num)
  }
  remove(patient_id,assay_list,id_num)
}

# Overwrite existing dataframe
DF_Full <- DF_Full_AMENDED
remove(DF_Full_AMENDED,DF_Full_AMENDED_pre)

# Use sys.date_created as alternative for missing smpl.dateReceived
date_NA <- which(is.na(DF_Full$smpl.dateReceived))
for (i in 1:length(date_NA)) {
  row_id = date_NA[i]
  patient_id = DF_Full$sys.uniqueId[row_id]
  DF_Full$smpl.dateReceived[row_id] <- unique(DF_Full$sys.date_created[DF_Full$sys.uniqueId == patient_id])
  
  remove(i,row_id,patient_id)
}
remove(date_NA)

## Remove entries with missing information
#----------------------------------------------
## Remove entries without gender
DF_Full <- DF_Full[which(!is.na(DF_Full$smpl.gender) & DF_Full$smpl.gender != "Unknown"),]

# Remove patients without DOB
DF_Full <- DF_Full[!is.na(DF_Full$base.dob),]

# Remove patients without pathogenicity
DF_Full <- DF_Full[!is.na(DF_Full$smpl.pathogenicityStatus),]

## Structure patient DOB and input current age
#----------------------------------------------
curr_year <- as.numeric(gsub("^([[:digit:]]{2})([[:digit:]]{2})", "\\2", format(as.Date(Sys.Date(), format="%Y-%m-%d"),"%Y")))

# Extract components of DOB
DF_Full$month = gsub("^0","",gsub("(^[[:digit:]]{,2})(.*)", "\\1", DF_Full$base.dob))
DF_Full$day = gsub("^0","",gsub("(^[[:digit:]]{,2}[/])([[:digit:]]{,2})(.*)", "\\2", DF_Full$base.dob))
DF_Full$year = sprintf("%02d", as.numeric(gsub("(^[[:digit:]]{,2}[/][[:digit:]]{,2}[/])([[:digit:]]{2}$)", "\\2", DF_Full$base.dob)))

# Assume no individual is >= 100yo
for (row_No in 1:nrow(DF_Full)) {
  if (isTRUE(DF_Full$year[row_No] > curr_year)) {
    DF_Full$year[row_No] <- paste("19", DF_Full$year[row_No], sep="")
  } else if (isTRUE(DF_Full$year[row_No] <= curr_year)) {
    DF_Full$year[row_No] <- paste("20", DF_Full$year[row_No], sep="")
  }
}

for (row_No in 1:nrow(DF_Full)) {
  DF_Full$patient.dob[row_No] <- paste(DF_Full$month[row_No], DF_Full$day[row_No], DF_Full$year[row_No], sep="/")
}
DF_Full$patient.dob <-  as.Date(DF_Full$patient.dob, "%m/%d/%Y")

# Age rounded down to nearest integer -- relative to smpl.dateReceived
DF_Full$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                          DF_Full$smpl.dateReceived), "%Y-%m-%d")

DF_Full$Age <- as.numeric(floor(age_calc(dob = DF_Full$patient.dob, 
                                         enddate = DF_Full$smpl.dateReceived, units = "years")))

# Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
#----------------------------------------------
DF_Full$smpl.assayName[which(DF_Full$smpl.assayName == "Myeloid Panel  (54 gene panel)")] <- "Myeloid Panel (54 gene panel)"
DF_Full$smpl.assayName[which(DF_Full$smpl.assayName == "Myeloid panel (54 gene panel)")] <- "Myeloid Panel (54 gene panel)"

A <- DF_Full[which(DF_Full$smpl.assayName == "Cancer Somatic Mutation Panel (48 gene panel)"), ]
cat("Cancer Somatic Mutation Panel (48 gene panel)","\n")
cat(paste(nrow(A), " total entries and ", length(unique(A[[1]])), " total patients", sep=""),"\n","\n")

A <- DF_Full[which(DF_Full$smpl.assayName == "Myeloid Panel (54 gene panel)"), ]
cat("Myeloid Panel (54 gene panel)","\n")
cat(paste(nrow(A), " total entries and ", length(unique(A[[1]])), " total patients", sep=""),"\n","\n")
remove(A)

DF_Full <- DF_Full[DF_Full$smpl.assayName %in% c("STAMP - Solid Tumor Actionable Mutation Panel (130 genes)",
                                                 "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)",
                                                 "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)"), ]
## Format dataframe 
#----------------------------------------------
# Rename columns
colnames_extract <- c("sys.uniqueId","smpl.gender","patient.dob","Age",
                      "smpl.histologicalDiagnosis","smpl.primaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "smpl.assayName","smpl.hasOrderingPhysician",
                      "sys.date_changed","sys.date_created","smpl.reportDateReviewed","smpl.dateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "sys.label","base.gene","smpl.pathogenicityStatus",
                      "smpl.hgvsCoding","smpl.hgvsProtein","smpl.transcript",
                      "smpl.genomicDescription","base.chromosome","smpl.hgvsGenomic","smpl.chromosomePositionStart")
DF_Full <- DF_Full[,colnames_extract]

colnames_generic <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                      "HistologicalDx","PrimaryTumorSite","smpl.specimenType","smpl.percentTumor",
                      "AssayName","AssayOrderingPhysician",
                      "sys.date_changed","sys.date_created","AssayReportDateReviewed","AssayDateReceived",
                      "smpl.amendedString","smpl.amendmentReason",
                      "VariantLabel","VariantGene","VariantPathogenicityStatus",
                      "VariantHGVSCoding","VariantHGVSProtein","VariantNMAccession",
                      "VariantHGVSGenomic","VariantCHR","smpl.hgvsGenomic","smpl.chromosomePositionStart")
colnames(DF_Full) <- colnames_generic

# General structurization
#----------------------------------------------
# Convert to lowercase 
DF_Full$HistologicalDx <- tolower(DF_Full$HistologicalDx)
DF_Full$PrimaryTumorSite <- tolower(DF_Full$PrimaryTumorSite)

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_DF", DF_Full, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(DF_Full, file = paste(tempdir, Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(phgvs.noeffect, phgvs.sys.label, row.change, x,patient_list,colnames_generic,colnames_extract,DF_Full,curr_year,row_No)
