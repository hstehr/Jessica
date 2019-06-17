HistologicalDxCategory <-
  read.csv(file = histoDx.key, header = TRUE, na.strings = NA, stringsAsFactors = FALSE, sep = "\t")

# Remove NA histological dx columns
HistologicalDxCategory <- HistologicalDxCategory[which(!is.na(HistologicalDxCategory$histologicalDiagnosis)),]

# Convert all text to lowercase
HistologicalDxCategory$histologicalDiagnosis <- tolower(HistologicalDxCategory$histologicalDiagnosis)
HistologicalDxCategory$primaryTumorSite <- tolower(HistologicalDxCategory$primaryTumorSite)
HistologicalDxCategory$CTEP.CATEGORY <- tolower(HistologicalDxCategory$CTEP.CATEGORY)
HistologicalDxCategory$CTEP.SUBCATEGORY <- tolower(HistologicalDxCategory$CTEP.SUBCATEGORY)
HistologicalDxCategory$CTEP.TERM <- tolower(HistologicalDxCategory$CTEP.TERM)
HistologicalDxCategory$SHORT.NAME <- tolower(HistologicalDxCategory$SHORT.NAME)

# Assign DF to global environment
assign("HistologicalDxCategory", HistologicalDxCategory, envir = .GlobalEnv)
