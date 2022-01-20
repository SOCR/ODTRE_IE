#######################################################################################
#
#                                                                                     
#   Filename    :	hyper_emergency_data_preprocessing.R    												  
#                                                                                     
#   Project     :       BiomJ article "Optimal dynamic treatment regime estimation using information extraction from unstructured clinical text"                                                             
#   Authors     :       Nina Zhou, Robert Brook, Ivo Dinov, and Lu Wang                                                                
#   Date        :       11/28/2021
#   Purpose     :       As described in BiomJ article
#																				  
#   R Version   :       4.0.2 (2020-06-22)                                                                
#   RStudio Version    :  	1.3.959
#   Command   :  	R CMD BATCH --no-save hyper_emergency_data_preprocessing.R
#
#   Input data  :    death.csv
#                    demo.csv
#                    demo_dx.csv
#                    drugs.csv
#                    firstday_lab.csv
#                    firstday_vital.csv
#                    hyper_text_ie.csv
#                    kidney_diabeties_copd_hyper.csv
#                    notes.csv
#                    vitals.csv
#   Output data files :   application_data_day1.csv, application_data_day2.csv (preprocessed data for application section)
#
#   Required R packages :  dplyr
#
########################################################################################
#Load package
library(dplyr)

#Baseline SBP
vital_firstday <-
  read.csv("firstday_vital.csv")

vital <- read.csv("vitals.csv")
vital$charttime <- as.POSIXct(vital$charttime, tz = "UTC")
vital$day <- as.Date(vital$charttime, format = "%Y-%m-%d")
vital$rownumb <- 1:nrow(vital)


count_day <-
  vital %>% group_by(icustay_id) %>% filter(!is.na(SysBP)) %>% 
  summarise(daycount = length(unique(day)))
#table(count_day$daycount) #select larger than 3 days.

secondday_vital <-
  vital %>% group_by(icustay_id) %>% filter(!is.na(SysBP)) %>%
  filter(icustay_id %in% count_day$icustay_id[count_day$daycount > 1]) %>% 
  filter(day != min(day)) %>%
  filter(day == min(day)) %>%
  summarize(
    day = max(day),
    sbp = max(SysBP, na.rm = TRUE),
    dbp = max(DiasBP, na.rm = TRUE),
    HR = mean(HeartRate, na.rm = TRUE),
    Temp = mean(TempC, na.rm = TRUE),
    SpO2 = mean(SpO2, na.rm = TRUE),
    Glucose = mean(Glucose, na.rm = TRUE)
  )

thirdday_vital <-
  vital %>% group_by(icustay_id) %>% filter(!is.na(SysBP)) %>% 
  filter(day != min(day)) %>% filter(day != min(day)) %>%
  filter(day == min(day)) %>% summarize(
    day = max(day),
    sbp = max(SysBP, na.rm = TRUE),
    dbp = max(DiasBP, na.rm = TRUE),
    HR = mean(HeartRate, na.rm = TRUE),
    Temp = mean(TempC, na.rm = TRUE),
    SpO2 = mean(SpO2, na.rm = TRUE),
    Glucose = mean(Glucose, na.rm = TRUE)
  )


#check how many patients died within 90days
demo <- read.csv("demo.csv")
demo <-
  demo %>% group_by(SUBJECT_ID, HADM_ID) %>% summarise(wt_avg = mean(WEIGHT, na.rm = TRUE),
                                                       ht_avg = mean(HEIGHT, na.rm = TRUE))

aa <- demo %>% filter(!is.na(wt_avg))

#Drugs
drug <- read.csv("drugs.csv")
drug$STARTDATE <- as.Date(drug$STARTDATE, format = "%Y-%m-%d")
drug$ENDDATE <- as.Date(drug$STARTDATE, format = "%Y-%m-%d")
names(drug)[5] <- "day"


hyper_drug_type <- function(vec) {
  out = ifelse(
    vec %in% c(
      "Atenolol",
      "Labetalol",
      "Carvedilol",
      "Propranolol",
      "Sotalol",
      "Nadolol",
      "Metoprolol",
      "Esmolol",
      "Acebutolol",
      "Betaxolol",
      "Bisoprolol Fumarate"
    ),
    "beta_blokers",
    ifelse(
      vec %in% c(
        "Bumetanide",
        "Chlorothiazide",
        "Chlorthalidone",
        "Hydrochlorothiazide",
        "Furosemide",
        "Metolazone",
        "Torsemide"
      ),
      "diuretics",
      ifelse(
        vec %in% c(
          "Captopril",
          "Lisinopril",
          "Moexipril",
          "Ramipril",
          "Trandolapril",
          "Enalapril Maleate",
          "Enalaprilat"
        ),
        "ace_inhibitors",
        ifelse(
          vec %in% c(
            "Felodipine",
            "Nifedipine",
            "Nimodipine",
            "Nicardipine",
            "Clevidipine",
            "Diltiazem"
          ),
          "calcium_channel_blockers",
          ifelse(vec %in% c("Nitroglycerin", "Nitroprusside"), "Nitrates", NA)
        )
      )
    )
  )
  return(out)
}

drug$class <- hyper_drug_type(drug$DRUG_NAME_GENERIC)

drug <-
  drug %>% filter(!(ROUTE %in% c("IV", "IV BOLUS", "IV DRIP") &
                      class != "diuretics")) %>% 
  filter(!(class == "diuretics" & !ROUTE %in% c("IV", "IV BOLUS", "IV DRIP")))

#Display drugs per patients per day.
drug <- drug %>% filter(!is.na(day))
drug_day <- data.frame()
for (i in 1:nrow(drug)) {
  if (is.na(drug[i, ]$ENDDATE) |
      drug[i, ]$ENDDATE < drug[i, ]$day)
    drug[i, ]$ENDDATE = drug[i, ]$day
  temp <-
    data.frame(
      SUBJECT_ID = drug[i, ]$SUBJECT_ID,
      day = seq(drug[i, ]$day, drug[i, ]$ENDDATE, by = "days"),
      class = drug[i, ]$class
    )
  drug_day <- rbind(drug_day, temp)
}

drug_day <-
  drug_day %>% group_by(SUBJECT_ID, day) %>% summarize(
    ace = ifelse("ace_inhibitors" %in% class, 1, 0),
    beta = ifelse("beta_blokers" %in% class, 1, 0),
    cal = ifelse("calcium_channel_blockers" %in% class, 1, 0),
    diuretics = ifelse("diuretics" %in% class, 1, 0)
  )
drug_day <- drug_day %>% filter(ace + beta + cal + diuretics > 0)

drug_count <-
  drug_day %>% group_by(SUBJECT_ID, day) %>% summarise(numb_drugs = ace +
                                                         beta + cal + diuretics)


#Text
hyper_ie <- read.csv("hyper_text_ie.csv")

hyper_ie <- hyper_ie %>% select(SUBJECT_ID, height, weight, smoke)
hyper_ie$height <- hyper_ie$height * 2.54
hyper_ie$weight <- hyper_ie$weight * 0.453592

#Here is how you come up with wt_avg investigate the missingness.
BB <- demo %>% left_join(hyper_ie, by = "SUBJECT_ID")
BB$wt_ori <- BB$wt_avg
BB$wt_avg <- ifelse(is.na(BB$wt_avg), BB$weight, BB$wt_avg)
BB$ht_avg <- ifelse(is.na(BB$ht_avg), BB$height, BB$ht_avg)
BB$smoke <- ifelse(is.na(BB$smoke), 0, BB$smoke)

temp_miss_wt = BB %>% filter(!is.na(wt_avg) &
                               is.na(wt_ori)) %>% select(SUBJECT_ID, HADM_ID) %>% distinct()

demo_ie <-
  BB %>% select(-height, -weight) %>% group_by(SUBJECT_ID, HADM_ID) %>%
  summarise(
    wt_avg = mean(wt_avg, na.rm = TRUE),
    ht_avg = mean(ht_avg, na.rm = TRUE),
    smoke = max(smoke, na.rm = TRUE)
  ) %>% distinct()

aa <- demo_ie %>% filter(!is.na(wt_avg))
ID_cmp_wt <- unique(aa$HADM_ID)#2220


#################### Merging ##############
#Merge demo_ie/vitals
vital_firstday <- vital_firstday %>% filter(hadm_id %in% ID_cmp_wt)
hadm_table <-
  vital_firstday %>% select(hadm_id, icustay_id) %>% distinct()
names(hadm_table)[1] <- "HADM_ID"

secondday_vital <-
  secondday_vital %>% filter(icustay_id %in% unique(vital_firstday$icustay_id))
thirdday_vital <-
  thirdday_vital %>% filter(icustay_id %in% unique(vital_firstday$icustay_id))

demo_ie <- demo_ie %>% left_join(hadm_table, by = "HADM_ID")

vital_demo_1 <-
  vital_firstday %>% left_join(demo_ie, by = c("icustay_id"))
vital_demo_2 <-
  secondday_vital %>% left_join(demo_ie, by = c("icustay_id"))
vital_demo_3 <-
  thirdday_vital %>% left_join(demo_ie, by = c("icustay_id"))

#Add drug info
#Day 1
day1 <-
  vital %>% group_by(SUBJECT_ID, icustay_id) %>% filter(!is.na(SysBP)) %>% 
  filter(day == min(day)) %>% summarise(day = min(day))
drug_day1 <-
  day1 %>% left_join(drug_day, by = c("SUBJECT_ID", "day")) %>% filter(!is.na(ace))


vital_demo_drug_1 <-
  vital_demo_1 %>% left_join(drug_day1, by = c("icustay_id"))
vital_demo_drug_1 <-
  vital_demo_drug_1 %>% filter(!is.na(ace) &
                                 icustay_id %in% vital_demo_2$icustay_id)
rownum <-
  unlist(sapply(vital_demo_drug_1$icustay_id, function(x)
    which(vital_demo_2 == x)))

vital_demo_drug_1$final_sbp <- vital_demo_2$sbp[rownum]
vital_demo_drug_1$final_dbp <- vital_demo_2$dbp[rownum]
#dim(vital_demo_drug_1)  #778 records

#Day 2
day2 <-
  vital %>% group_by(SUBJECT_ID, icustay_id) %>% filter(!is.na(SysBP)) %>% 
  filter(day != min(day)) %>% filter(day == min(day)) %>% summarise(day = min(day))
drug_day2 <-
  day2 %>% left_join(drug_day, by = c("SUBJECT_ID", "day")) %>% filter(!is.na(ace))

vital_demo_drug_2 <-
  vital_demo_2 %>% left_join(drug_day2, by = c("icustay_id"))
vital_demo_drug_2 <-
  vital_demo_drug_2 %>% filter(!is.na(ace) &
                                 icustay_id %in% vital_demo_3$icustay_id)

rownum <-
  unlist(sapply(vital_demo_drug_2$icustay_id, function(x)
    which(vital_demo_3 == x)))

vital_demo_drug_2$final_sbp <- vital_demo_2$sbp[rownum]
vital_demo_drug_2$final_dbp <- vital_demo_2$dbp[rownum]
#dim(vital_demo_drug_2) #748 records

###############

###Add age/death/gender
adg <- read.csv("death.csv")
adg$DOB <- as.Date(adg$DOB, format = "%Y-%m-%d")
adg$DOD <- as.Date(adg$DOD, format = "%Y-%m-%d")
adg$GENDER <- ifelse(adg$GENDER == "M", 0, 1)
vital_demo_drug_1 <-
  vital_demo_drug_1 %>% left_join(adg, by = c("subject_id" = "SUBJECT_ID"))
vital_demo_drug_1$age <-
  difftime(vital_demo_drug_1$day, vital_demo_drug_1$DOB, units = "days") /
  365
vital_demo_drug_1$age <-
  ifelse(vital_demo_drug_1$age > 300,
         as.numeric(mean(vital_demo_drug_1$age)),
         as.numeric(vital_demo_drug_1$age))
vital_demo_drug_1$ttdeath <- ifelse(
  is.na(vital_demo_drug_1$DOD),
  NA,
  difftime(vital_demo_drug_1$DOD, vital_demo_drug_1$day, units = "days")
)

vital_demo_drug_2 <-
  vital_demo_drug_2 %>% left_join(vital_demo_drug_1 %>% select(icustay_id, age, GENDER), by =
                                    "icustay_id")
vital_demo_drug_2 <- vital_demo_drug_2 %>% filter(!is.na(GENDER))

#Add kidney + diabetes + copd
kdc <-
  read.csv("kidney_diabeties_copd_hyper.csv")
kdc$kidney <- grepl("kidney|Kidney", kdc$LONG_TITLE)
kdc$diabetes <- grepl("diabetes|Diabetes", kdc$LONG_TITLE)
kdc$copd <- grepl("bronchitis", kdc$LONG_TITLE)
kdc$hyper <- grepl("hypertension", kdc$LONG_TITLE)

vital_demo_drug_1$kidney <-
  ifelse(vital_demo_drug_1$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$kidney ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_1$diabetes <-
  ifelse(vital_demo_drug_1$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$diabetes ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_1$copd <-
  ifelse(vital_demo_drug_1$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$copd ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_1$hyper <-
  ifelse(vital_demo_drug_1$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$hyper ==
                                                                    TRUE)], 1, 0)

vital_demo_drug_2$kidney <-
  ifelse(vital_demo_drug_2$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$kidney ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_2$diabetes <-
  ifelse(vital_demo_drug_2$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$diabetes ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_2$copd <-
  ifelse(vital_demo_drug_2$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$copd ==
                                                                    TRUE)], 1, 0)
vital_demo_drug_2$hyper <-
  ifelse(vital_demo_drug_2$SUBJECT_ID.x %in% kdc$SUBJECT_ID[which(kdc$hyper ==
                                                                    TRUE)], 1, 0)

#Add lab results first day
####Some how adding this shrinks the population.
lab <- read.csv("firstday_lab.csv")
lab <- lab %>% filter(icustay_id %in% vital_demo_drug_1$icustay_id)
lab$wbc_min[which(is.na(lab$wbc_min))] = mean(lab$wbc_min, na.rm = TRUE)
lab$wbc_max[which(is.na(lab$wbc_max))] = mean(lab$wbc_max, na.rm = TRUE)
lab$hemoglobin_max[which(is.na(lab$hemoglobin_max))] = mean(lab$hemoglobin_max, na.rm = TRUE)
vital_demo_drug_1 <- vital_demo_drug_1 %>%
  left_join(lab %>% select(icustay_id, hemoglobin_max, creatinine_min), by =
              "icustay_id")

demo_dx <- read.csv("demo_dx.csv")
demo_dx <- demo_dx %>% filter(HADM_ID %in% vital_demo_drug_1$hadm_id)
demo_dx$black <- ifelse(grepl("BLACK", demo_dx$ETHNICITY), 1, 0)
names(demo_dx)[2] <- "hadm_id"
vital_demo_drug_1 <-
  vital_demo_drug_1 %>% left_join(demo_dx %>% select(hadm_id, black), by =
                                    "hadm_id")
names(vital_demo_drug_2)[which(names(vital_demo_drug_2) == 'HADM_ID')] = 'hadm_id'
vital_demo_drug_2 <-
  vital_demo_drug_2 %>% left_join(demo_dx %>% select(hadm_id, black), by =
                                    "hadm_id")

write.csv(vital_demo_drug_2, "application_data_day2.csv", row.names = FALSE)
write.csv(vital_demo_drug_1, "application_data_day1.csv", row.names = FALSE)