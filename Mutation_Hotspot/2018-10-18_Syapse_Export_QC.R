## Confirmed for timestamp = c("2018-10-18")

## QC-parameters: smpl.assayName, smpl.pipelineVersion, base.gene, smpl.hgvsProtein, smpl.hgvsCoding
## Calculate patient age > merge with primary tumor site data > filter for STAMP entries 
## Output: "Syapse_Export_QC.tsv"

cat(paste("Timestamp of Syapse_Report: ", Syapse_Export_timestamp, sep=""),"\n")

# Load relevant file
#----------------------------------------------
DF_Full <- 
  read.csv(file = paste("~/Documents/ClinicalDataScience_Fellowship/STAMP/", 
                        Syapse_Export_timestamp, "_syapse_export_all_variants_patientNameAndMrnRemoved.csv", sep=""),
                    header = TRUE, na.strings = c(""," ","NA"), stringsAsFactors = FALSE, sep = ",")

# Shorten column names of DF
colnames(DF_Full) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_Full))

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
if (Syapse_Export_timestamp  == "2018-10-18") {
  
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
  
  # Change to same capitalization scheme for consistency
  #----------------------------------------------
  # which(DF_Full$smpl.hgvsProtein == "PROMOTER")
  DF_Full$smpl.hgvsProtein[c(12743,12762,12763,12764)] <- "Promoter"
  phgvs.noeffect <- which(DF_Full$smpl.hgvsProtein == "Promoter")
  phgvs.sys.label <- unique(DF_Full[phgvs.noeffect, "sys.label"])
  for (x in 1:length(phgvs.sys.label)) {
    row.change <- which(DF_Full$sys.label == phgvs.sys.label[x])
    DF_Full[row.change, "smpl.hgvsProtein"]  <- NA
  }
  
  # Correct annotation 
  #----------------------------------------------
  # NM_005157.5(ABL1):c.949T>C (p.Phe317Leu)
  # which(DF_Full$smpl.hgvsCoding == "p.Phe317Leu")
  # DF_Full[12466,c(1,50,53:61)]
  DF_Full[12466,"smpl.hgvsCoding"] <- "c.949T>C"
  
  # COSM211736 (c.4636C>T) (p.Q1546*)
  # which(DF_Full$smpl.hgvsProtein == "p.Q1546X")
  DF_Full[c(1411,1412),"smpl.hgvsProtein"] <- "p.Gln1546Ter"
  
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
  # which(DF_Full$smpl.hgvsProtein == ":p.Phe221fs")
  DF_Full$smpl.hgvsProtein[6200] <- "p.Phe221fs"
  # which(DF_Full$smpl.hgvsProtein == "p.A129fs")
  DF_Full$smpl.hgvsProtein[4082] <- "p.Ala129fs"
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
  
  DF_Full$smpl.genomicDescription[10307] <- "chr17:g.7572926_7572929delGTCA"
  
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
}
###################################
## Manual edit END
###################################

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

## Subset columns of interest
#----------------------------------------------
colnames_keep <- c("sys.uniqueId","smpl.gender","base.dob","smpl.dateReceived",
                   "smpl.hasOrderingPhysician","smpl.csmAssay","smpl.dateCollected","smpl.reportDateReviewed",
                   "smpl.specimenSite","smpl.ppDiagnosticSummary",
                   "smpl.assayName","smpl.transcript","base.chromosome","smpl.hgvsGenomic",
                   "sys.label","base.gene","smpl.hgvsCoding","smpl.hgvsProtein","smpl.pathogenicityStatus")
DF_Full_Filter <- DF_Full[,colnames_keep]

# Use sys.date_created as alternative for missing smpl.dateReceived
date_NA <- which(is.na(DF_Full_Filter$smpl.dateReceived))
for (i in 1:length(date_NA)) {
  row_id = date_NA[i]
  DF_Full_Filter$smpl.dateReceived[row_id] <- 
    unique(DF_Full$sys.date_created[DF_Full$sys.uniqueId == DF_Full_Filter$sys.uniqueId[row_id]])
}
DF_Full <- DF_Full_Filter

## Remove entries with missing information
#----------------------------------------------
## Remove entries without gender = 23 entries
DF_Full <- DF_Full[which(!is.na(DF_Full$smpl.gender) & 
                           DF_Full$smpl.gender != "Unknown"),]

# Remove patients without DOB = 4 entries
DF_Full <- DF_Full[!is.na(DF_Full$base.dob),]

## Structure patient DOB and input current age
#----------------------------------------------
curr_year <- as.numeric(gsub("^([[:digit:]]{2})([[:digit:]]{2})", "\\2", 
                             format(as.Date(Sys.Date(), format="%Y-%m-%d"),"%Y")))

# Extract components of DOB
DF_Full$month = gsub("(^[[:digit:]]{,2})(.*)", "\\1", DF_Full$base.dob)
DF_Full$day=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})(.*)", "\\3", DF_Full$base.dob)
DF_Full$year=gsub("(^[[:digit:]]{,2})([/])([[:digit:]]{,2})([/])([[:digit:]]{2})[[:blank:]]*$", "\\5", 
                  DF_Full$base.dob)

# Assume no individual is >= 100yo
for (row_No in 1:nrow(DF_Full)) {
  if (DF_Full$year[row_No] > curr_year) {
    DF_Full$year[row_No] <- paste("19", DF_Full$year[row_No], sep="")
  } else if (DF_Full$year[row_No] <= curr_year) {
    DF_Full$year[row_No] <- paste("20", DF_Full$year[row_No], sep="")
  }}

for (row_No in 1:nrow(DF_Full)) {
  DF_Full$patient.dob <- as.Date(paste(DF_Full$month, DF_Full$day, DF_Full$year, sep="/"), "%m/%d/%Y")
}

# Age rounded down to nearest integer -- relative to smpl.dateReceived
DF_Full$smpl.dateReceived <- as.Date(gsub("(^[[:digit:]]{4}[-][[:digit:]]{2}[-][[:digit:]]{2})(.*)", "\\1", 
                                          DF_Full$smpl.dateReceived), "%Y-%m-%d")

DF_Full$Age <- as.numeric(floor(age_calc(dob = DF_Full$patient.dob, 
                                         enddate = DF_Full$smpl.dateReceived, units = "years")))

## Merge with primary tumor site data
#----------------------------------------------
if (Syapse_Export_timestamp  == "2018-10-18") {
  
  DF_TumorSite <- read.csv(file = "~/Documents/ClinicalDataScience_Fellowship/STAMP/2018-12-17_syapse_primary_tumor_site.csv",
                           header = TRUE, na.strings = c("NA","None"), stringsAsFactors = FALSE,sep = ",")
  # Remove extraneous columns
  DF_TumorSite <- DF_TumorSite[,c(1,2,7)]
  # Shorten column names of DF
  colnames(DF_TumorSite) <- gsub("smpl.[a-zA-Z]+[.]{3}", "", colnames(DF_TumorSite))
  
  # Remove duplicate rows   
  DF_TumorSite <- DF_TumorSite %>% dplyr::distinct(sys.uniqueId, .keep_all = TRUE)
  
  # Merge with STAMP entries 
  DF_Full <- left_join(DF_Full, DF_TumorSite, by = c("sys.uniqueId"))
  
  # Filter entries from STAMP - Solid Tumor Actionable Mutation Panel
  #----------------------------------------------
  DF_Full <- DF_Full[DF_Full$smpl.assayName %in% c("STAMP - Solid Tumor Actionable Mutation Panel (130 genes)",
                                                   "STAMP - Solid Tumor Actionable Mutation Panel (198 genes)",
                                                   "STAMP - Solid Tumor Actionable Mutation Panel (200 genes)"), ]
  cat(paste("Total Count: n=",(nrow(DF_Full)), 
            " STAMP entries (n=", length(unique(DF_Full$sys.uniqueId)), " patients)",sep=""),"\n")
}

## Format dataframe 
#----------------------------------------------
# Select columns of interest
colnames_keep <- c("sys.uniqueId","smpl.gender","patient.dob","Age","smpl.dateReceived",
                   "smpl.hasOrderingPhysician","smpl.csmAssay","smpl.dateCollected","smpl.reportDateReviewed",
                   "smpl.specimenSite","smpl.ppDiagnosticSummary",
                   "smpl.assayName","smpl.transcript","base.chromosome","smpl.hgvsGenomic",
                   "sys.label","base.gene","smpl.hgvsCoding","smpl.hgvsProtein","smpl.pathogenicityStatus",
                   "smpl.primaryTumorSite","smpl.histologicalDiagnosis")
DF_Full <- DF_Full[,colnames_keep]

# Rename columns
colnames_generic <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                      "AssayDateReceived","AssayOrderingPhysician","AssayType",
                      "SpecimenDateCollected","AssayReportDateReviewed","SpecimenSite","DxSummary","AssayName",
                      "VariantNMAccession","VariantCHR","VariantHGVSGenomic","VariantLabel","VariantGene",
                      "VariantHGVSCoding","VariantHGVSProtein","VariantPathogenicityStatus",
                      "PrimaryTumorSite","HistologicalDx")
colnames(DF_Full) <- colnames_generic

colnames_order <- c("PatientID","PatientGender","PatientDOB","PatientAge",
                    "AssayOrderingPhysician","AssayDateReceived","AssayType","AssayName","AssayReportDateReviewed",
                    "SpecimenDateCollected","SpecimenSite","DxSummary","HistologicalDx","PrimaryTumorSite",
                    "VariantNMAccession","VariantCHR","VariantHGVSGenomic","VariantLabel","VariantGene",
                    "VariantHGVSCoding","VariantHGVSProtein","VariantPathogenicityStatus")
DF_Full <- DF_Full[,colnames_order]

# General structurization
#----------------------------------------------
# Convert to lowercase 
DF_Full$PrimaryTumorSite <- tolower(DF_Full$PrimaryTumorSite)

## Overwrite variable in global environment
#----------------------------------------------
assign("STAMP_DF", DF_Full, envir = .GlobalEnv)

## Write to local computer
#----------------------------------------------
write.table(DF_Full, file = paste(getwd(), "/", Syapse_Export_timestamp, "_Syapse_Export_QC.tsv", sep=""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

remove(phgvs.noeffect, phgvs.sys.label, row.change, x,DF_Full_AMENDED, 
       DF_Full_AMENDED_pre,assay_list,assay_id,assay_num,patient_list,patient_id,id_num,
       created_list,date_NA,row_id,i,DF_Full_Filter,DF_TumorSite,colnames_generic,
       colnames_keep,colnames_order,DF_Full,curr_year,row_No)
