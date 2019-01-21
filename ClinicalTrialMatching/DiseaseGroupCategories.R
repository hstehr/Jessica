print("Melanoma and Blook have NOT been classified into a primaryTumorSite")

#################################
## Customize disease group categories 
#################################
# primaryTumorSite and Disease.Group classified according to medical specializations 
# all unique Disease.Site entries should be included in primaryTumorSite column 
Breast <- data.frame(Disease.Group.category = "Breast",
                     primaryTumorSite =  c("Breast"),
                     Disease.Group = NA,
                     stringsAsFactors = FALSE)

Dermatology <- data.frame(Disease.Group.category = "Dermatology",
                          primaryTumorSite =  c("Back","Foot","Scalp","Skin"),
                          Disease.Group = c("Cutaneous Oncology"),
                          stringsAsFactors = FALSE)

Endocrinology <- data.frame(Disease.Group.category = "Endocrinology",
                 primaryTumorSite =  c("Adrenal"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

Gastroenterology <- data.frame(Disease.Group.category = "Gastroenterology",
                 primaryTumorSite =  c("Ampulla of Vater","Anus","Appendix","Colon","Colon and Rectum",
                                       "Gallbladder","Hepatocellular (Liver)","Intrahepatic Bile Ducts",
                                       "Liver","Mesenteric Mass","Pancreas","Perihilar Bile Ducts",
                                       "Peritoneum","Small Intestine","Stomach"),
                 Disease.Group = c("Gastrointestinal Oncology"),
                 stringsAsFactors = FALSE)

Genitourinary <- data.frame(Disease.Group.category = "Genitourinary",
                 primaryTumorSite =  c("Foreskin","Kidney","Penis","Prostate Gland","Testes","Testis",
                                       "Urethra","Urinary Bladder"),
                 Disease.Group = c("Genitourinary Oncology"),
                 stringsAsFactors = FALSE)

Gynecology <- data.frame(Disease.Group.category = "Gynecology",
                 primaryTumorSite =  c("Cervix","Ovary","Uterus","Vulva"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

HematoLymphoid <- data.frame(Disease.Group.category = "HematoLymphoid",
                 primaryTumorSite =  c("Bone Marrow","Hematologic and Lymphatic Neoplasm","Lymph node",
                                       "Peripheral blood"),
                 Disease.Group = c("Hematology"),
                 stringsAsFactors = FALSE)

Melanoma <- data.frame(Disease.Group.category = "Melanoma",
                 primaryTumorSite =  NA,
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

Neurology <- data.frame(Disease.Group.category = "Neurology",
                 primaryTumorSite =  c("Central Nervous System (Brain/Spinal Cord)"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

Ophthalmology <- data.frame(Disease.Group.category = "Ophthalmology",
                 primaryTumorSite =  c("Ocular Tumors","Uvea"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

Otolaryngology <- data.frame(Disease.Group.category = "Otolaryngology",
                 primaryTumorSite =  c("Esophagus","Larynx","Lip and Oral Cavity",
                                       "Nasal Cavity and Paranasal Sinuses","Nasopharynx","Neck",
                                       "Pharynx","Salivary Glands","Thymus","Thyroid"),
                 Disease.Group = c("Head & Neck Oncology"),
                 stringsAsFactors = FALSE)

Pulmonology <- data.frame(Disease.Group.category = "Pulmonology",
                 primaryTumorSite =  c("Lung","Pleura","Trachea"),
                 Disease.Group = c("Thoracic Oncology"),
                 stringsAsFactors = FALSE)

Sarcoma <- data.frame(Disease.Group.category = "Sarcoma",
                 primaryTumorSite =  c("Foot","Back", "Bone","Soft Tissue"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

Unknown <- data.frame(Disease.Group.category = "Unknown",
                 primaryTumorSite =  c("Other Primary Site","Unknown"),
                 Disease.Group = NA,
                 stringsAsFactors = FALSE)

AnySite <- data.frame(Disease.Group.category = "Any Site",
                      primaryTumorSite =  NA,
                      Disease.Group = c("Non-CRG Specific", "Developmental Therapeutics"),
                      stringsAsFactors = FALSE)

DF_DiseaseGroupCategory <- rbind(Breast,Dermatology,Endocrinology,Gastroenterology,
                                 Genitourinary,Gynecology,HematoLymphoid,Melanoma,Neurology,
                                 Ophthalmology,Otolaryngology,Pulmonology,Sarcoma,Unknown,AnySite)
remove(Breast,Dermatology,Endocrinology,Gastroenterology,Genitourinary,Gynecology,HematoLymphoid,
       Melanoma,Neurology,Ophthalmology,Otolaryngology,Pulmonology,Sarcoma,Unknown,AnySite)

# Remove categories with NA columns
DF_DiseaseGroupCategory <- DF_DiseaseGroupCategory[!(is.na(DF_DiseaseGroupCategory$primaryTumorSite) & 
                                                       is.na(DF_DiseaseGroupCategory$Disease.Group)),]

# Convert all text to lowercase
DF_DiseaseGroupCategory$Disease.Group.category <- tolower(DF_DiseaseGroupCategory$Disease.Group.category)
DF_DiseaseGroupCategory$primaryTumorSite <- tolower(DF_DiseaseGroupCategory$primaryTumorSite)
DF_DiseaseGroupCategory$Disease.Group <- tolower(DF_DiseaseGroupCategory$Disease.Group)

assign("DiseaseGroupCategory_LongFormat", DF_DiseaseGroupCategory, envir = .GlobalEnv)

# Convert to wide-format
DF_DiseaseGroupCategory_Final <- data.frame()
category.list <- sort(unique(DF_DiseaseGroupCategory$Disease.Group.category))

for (row_No in 1:length(category.list)) {
  DF_Temp <- DF_DiseaseGroupCategory[DF_DiseaseGroupCategory$Disease.Group.category == category.list[row_No],]
  primaryTumorSite.list <- sort(unique(DF_Temp$primaryTumorSite))
  Disease.Group.list <- sort(unique(DF_Temp$Disease.Group))
  
  DF_DiseaseGroupCategory_pre <- data.frame(Disease.Group.category = category.list[row_No],
                                            primaryTumorSite = paste(as.character(primaryTumorSite.list), collapse=", "),
                                            Disease.Group = paste(as.character(Disease.Group.list), collapse=", "),
                                            stringsAsFactors = FALSE)
  DF_DiseaseGroupCategory_Final <- rbind(DF_DiseaseGroupCategory_Final, DF_DiseaseGroupCategory_pre)
}

# Convert empty cells to NA
DF_DiseaseGroupCategory_Final$primaryTumorSite[DF_DiseaseGroupCategory_Final$primaryTumorSite == ""] <- NA
DF_DiseaseGroupCategory_Final$Disease.Group[DF_DiseaseGroupCategory_Final$Disease.Group == ""] <- NA

# Assign DF to global environment
assign("DiseaseGroupCategory", DF_DiseaseGroupCategory_Final, envir = .GlobalEnv)

remove(DF_DiseaseGroupCategory,DF_DiseaseGroupCategory_pre,DF_Temp,category.list,Disease.Group.list,
       primaryTumorSite.list,row_No,DF_DiseaseGroupCategory_Final)
