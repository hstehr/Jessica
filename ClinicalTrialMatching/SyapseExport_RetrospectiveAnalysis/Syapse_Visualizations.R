suppressMessages(libraries("gridExtra","gtable","grid","ggrepel","RColorBrewer"))

# Specify outdir
#----------------------------------------------
visual.outdir = paste(outdir.root,"TIFF/",sep="")
if (!dir.exists(visual.outdir)){dir.create(visual.outdir)} 

# Specify distribution of timeline
#----------------------------------------------
month.number <- c("2014-01","2014-02","2014-03","2014-04","2014-05","2014-06","2014-07","2014-08","2014-09","2014-10","2014-11","2014-12",
                  "2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08","2015-09","2015-10","2015-11","2015-12",
                  "2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08","2016-09","2016-10","2016-11","2016-12",
                  "2017-01","2017-02","2017-03","2017-04","2017-05","2017-06","2017-07","2017-08","2017-09","2017-10","2017-11","2017-12",
                  "2018-01","2018-02","2018-03","2018-04","2018-05","2018-06","2018-07","2018-08","2018-09","2018-10","2018-11","2018-12",
                  "2019-01","2019-02","2019-03","2019-04")

month.alpha <- c("Jan 2014","Feb 2014","Mar 2014","Apr 2014","May 2014","Jun 2014","Jul 2014","Aug 2014","Sep 2014","Oct 2014","Nov 2014","Dec 2014",
                 "Jan 2015","Feb 2015","Mar 2015","Apr 2015","May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
                 "Jan 2016","Feb 2016","Mar 2016","Apr 2016","May 2016","Jun 2016","Jul 2016","Aug 2016","Sep 2016","Oct 2016","Nov 2016","Dec 2016",
                 "Jan 2017","Feb 2017","Mar 2017","Apr 2017","May 2017","Jun 2017","Jul 2017","Aug 2017","Sep 2017","Oct 2017","Nov 2017","Dec 2017",
                 "Jan 2018","Feb 2018","Mar 2018","Apr 2018","May 2018","Jun 2018","Jul 2018","Aug 2018","Sep 2018","Oct 2018","Nov 2018","Dec 2018",
                 "Jan 2019","Feb 2019","Mar 2019","Apr 2019")

## STAMP database: age distribution by gender
#----------------------------------------------
col_extract <- c("PatientID","PatientGender","PatientAge")
gender_age_AbsValue_fxn (DF = rbind(STAMP_DF_plot[,col_extract],
                                        STAMP_Fusion_plot[,col_extract],
                                        STAMP_CNV_plot[!is.na(STAMP_CNV_plot$PatientAge),col_extract]),
                             cohort=cohort_id,
                             outdir = visual.outdir)

## STAMP database: pathogenicity status + variant type
#----------------------------------------------
variant_type_distribution_fxn (DF_SNVIndel = STAMP_DF_plot,
                               DF_Fusion = STAMP_Fusion_plot,
                               DF_CNV = STAMP_CNV_plot,
                               cohort=cohort_id,
                               outdir = visual.outdir)


## STAMP database: top primary tumor sites
#----------------------------------------------
col_extract <- c("PatientID","PrimaryTumorSite")
top_site_count_fxn (DF = rbind(STAMP_DF_plot[,col_extract],
                           STAMP_Fusion_plot[,col_extract],
                           STAMP_CNV_plot[!is.na(STAMP_CNV_plot$PrimaryTumorSite),col_extract]),
                cohort=cohort_id,
                outdir = visual.outdir)


## STAMP database: test volume
#----------------------------------------------
col_extract <- c("PatientID","AssayDateReceived")
test_volume_timeline_fxn (DF = rbind(STAMP_DF_plot[,col_extract],
                                     STAMP_Fusion_plot[,col_extract],
                                     STAMP_CNV_plot[!is.na(STAMP_CNV_plot$AssayDateReceived),col_extract]),
                          cohort=cohort_id, outdir = visual.outdir, 
                          PerSite = FALSE)

remove(col_extract,month.number,month.alpha)
