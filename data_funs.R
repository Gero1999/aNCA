pk_aval <- function(time, Tmax, Cmax, lambda_z) {
  # Calculate the concentration at each time point
  concentration <- ifelse(
    time <= Tmax,
    # Absorption phase (linear rise to Cmax)
    Cmax * (time / Tmax),
    # Elimination phase (exponential decay after Tmax)
    Cmax * exp(-lambda_z * (time - Tmax))
  )
  concentration <- ifelse(time < 0, 0, concentration)
  # Return the results as a data frame
  concentration
}

add_sample_per_dose <- function(adnca,  nrrlt, sd = 0.05) {
  new_samples <- adnca %>%
    arrange(USUBJID, AFRLT) %>%
    group_by(USUBJID, PARAM, PCSPEC, DOSNOP) %>%
    slice(1) %>%
    group_by(USUBJID, DOSNOP) %>%
    mutate(
      nominal.dose.time = (NFRLT - NRRLT)[1],
      actual.dose.time = (AFRLT - ARRLT)[1],
      NRRLT = nrrlt,
      ARRLT = NRRLT + rnorm(n = n(), mean = 0, sd = sd),
      NFRLT = nominal.dose.time + NRRLT,
      AFRLT = actual.dose.time + ARRLT
    ) %>%
    ungroup()
  bind_rows(adnca, new_samples) %>%
    arrange(USUBJID, PCSPEC, DOSETRT, PARAM, AFRLT, NFRLT)
}

# Define subjects and number of doses
adnca <- crossing(
  STUDYID = "S1",
  USUBJID = paste0("S1-", sprintf("%02d", 1:8)),
  DOSNOP = 1:8
) %>%
  # Add other dose and subject information
  mutate(
    SEX = ifelse(USUBJID %in% paste0("S1-", sprintf("%02d", c(1,2,3,5,8))), "M", "F"),
    RACE = ifelse(USUBJID %in% paste0("S1-", sprintf("%02d", c(1,4,5,8,9))), "WHITE", "ASIAN"),
    AGE = case_when(
      USUBJID %in% paste0("S1-", sprintf("%02d", c(1,4))) ~ 25,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(2,5))) ~ 30,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(6))) ~ 35,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(3))) ~ 40,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(9))) ~ 45,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(7))) ~ 50,
      USUBJID %in% paste0("S1-", sprintf("%02d", c(8))) ~ 60,
    ),
    AGEU = "Years",
    PARAM = "DrugA",
    PCSPEC = "SERUM",
    ROUTE = "INTRAVENOUS DRIP",
    ADOSEDUR = 3,
    NDOSEDUR = 3,
    RRLTU = "Hours",
    DOSEU = "mg",
    ATPTREF = paste0("DOSE ", DOSNOP),
  ) %>%
  
  # Define different treatment doses
  mutate(
    DOSEA = ifelse(USUBJID %in% unique(USUBJID)[c(1,3,5,7,9)], 5, 10),
    TRT01A = paste0("DrugA ", DOSEA, " mg, Infusion"),
    DOSETRT = "DrugA",
    TRTRINT = 24
  ) %>%
  
  group_by(USUBJID, DOSNOP) %>%
  
  # Provide dose information
  mutate(
    nominal.dose.time = (TRTRINT * (DOSNOP - 1)),
    actual.dose.time = nominal.dose.time + rnorm(n(), mean = 0, sd = 0.1),
    actual.dose.time = nominal.dose.time + rnorm(n(), mean = 0, sd = 0.1),
    actual.dose.time = ifelse(actual.dose.time < 0, 0, actual.dose.time),
    actual.dose.time = actual.dose.time[1],
    ADOSEDUR = (NDOSEDUR[1] + rnorm(n(), mean = 0, sd = 0.05))[1]
  ) %>%
  ungroup() %>%
  group_by(USUBJID, DOSNOP) %>%
  mutate(rn = row_number()) %>%
  mutate(
    NRRLT = round(rn * 1.5, 2),
    ARRLT = (NRRLT + rnorm(n(), mean = 0, sd = 0.05))[1],
    NFRLT = nominal.dose.time + NRRLT,
    AFRLT = actual.dose.time + ARRLT
  ) %>%
  ungroup() %>%
  
  add_sample_per_dose(-5/60, sd = 0.001) %>%
  add_sample_per_dose(0.5) %>%
  add_sample_per_dose(1.5) %>%
  add_sample_per_dose(3) %>%
  add_sample_per_dose(4.5) %>%
  add_sample_per_dose(6) %>%
  add_sample_per_dose(7.5) %>%
  add_sample_per_dose(12) %>%
  
  mutate(
    AFRLT = ifelse(ATPTREF == "DOSE 1", ARRLT, AFRLT)
  ) %>%
  
  # Define concentration mainly based on dose and time
  mutate(
    AVAL = pk_aval(ARRLT, Tmax = ADOSEDUR, Cmax = 8 * (DOSEA/10), lambda_z = 0.5),
    AVALU = "ug/mL"
  ) %>%
  mutate(
    # Introduce variability per individual
    AVAL = AVAL * (ifelse(SEX == "M", 1, 1.1)),
    AVAL = AVAL * (ifelse(RACE == "WHITE", 1, 1.15)),
    AVAL = AVAL * (ifelse(USUBJID %in% paste0("S1-", sprintf("%02d", c(2,4,5,8))), 0.95, 1)),
    # Introduce sample error
    AVAL = AVAL + rnorm(n = nrow(.), mean = 0, sd = 0.05),
    AVAL = ifelse(AVAL < 0, 0, AVAL)
  ) %>%
  # Add other sample information
  mutate(
    METABFL = "",
    PCSTRESC = as.character(AVAL),
    AEFRLT = AFRLT,
    NEFRLT = NFRLT,
    AERRLT = ARRLT,
    NERRLT = NRRLT
  ) %>%
  # Make all numeric non-integer columns rounded to 3 decimals
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Add metabolite
adnca_metab <- adnca %>%
  mutate(
    PARAM = "Metab-DrugA",
    AVAL = AVAL * 0.4,
    METABFL = "Y"
  )

adnca <- bind_rows(adnca, adnca_metab)

# Introduce IV bolus doses for all doses over 5
adnca <- adnca %>%
  mutate(is.bolus = USUBJID %in% unique(USUBJID)[6:10]) %>%
  mutate(
    ADOSEDUR = ifelse(is.bolus, 0, ADOSEDUR),
    ROUTE = ifelse(is.bolus, "INTRAVENOUS BOLUS", "INTRAVENOUS DRIP"),
    TRT01A = ifelse(is.bolus, paste0("DrugA ", DOSEA, " mg, bolus"), TRT01A),
    DOSETRT = "DrugA",
    AVAL = ifelse(is.bolus, pk_aval(ARRLT, 0, 10 * (DOSEA/max(DOSEA)), 0.5), AVAL)
  ) %>%
  mutate(
    # Introduce variability per individual
    AVAL = ifelse(is.bolus, AVAL * (ifelse(SEX == "M", 1, 1.1)), AVAL),
    
    # Introduce Sample error
    AVAL = AVAL + rnorm(n = nrow(.), mean = 0, sd = 0.05),
    AVAL = ifelse(AVAL < 0, 0, AVAL),
    
    PCSTRESC = as.character(AVAL)
  )

adnca_urine <- adnca %>%
  filter(USUBJID %in% unique(USUBJID)[1:5]) %>%
  filter(METABFL == "") %>%
  mutate(
    PCSPEC = "URINE",
    VOLUME = round(rnorm(n = nrow(.), mean = 100, sd = 15)),
    VOLUMEU = "mL",
    NEFRLT = NFRLT + 2,
    NERRLT = NRRLT + 2,
    AEFRLT = round(NEFRLT + rnorm(n = nrow(.), mean = 2, sd = 0.1), 1),
    AERRLT = ARRLT + (AEFRLT - AFRLT),
    AVALU = "mg/mL",
    ATPTREF = paste0("DOSE ", DOSNOP)
  )

new_adnca <- dplyr::bind_rows(adnca, adnca_urine) %>%
  # Simplify some values
  mutate(USUBJID = gsub("111", "", USUBJID)) %>%
  
  # Order columns
  select(
    any_of(
      c("STUDYID", "USUBJID", "PCSPEC", "PARAM", "METABFL", "AFRLT", "NFRLT", "ARRLT", "NRRLT", "TRTRINT", "RRLTU", "AVAL", "AVALU", "VOLUME", "VOLUMEU", "DOSETRT", "TRT01A", "DOSEA", "DOSEU", "ROUTE", "ADOSEDUR", "ATPTREF", "AVISIT", "AGE", "AGEU", "RACE", "SEX")
    )
  ) %>%
  # IV bolus will have time duplicates, pre-remove them
  group_by(STUDYID, USUBJID, PARAM, PCSPEC, AFRLT) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(STUDYID, USUBJID, DOSETRT, PARAM, PCSPEC, AFRLT, NFRLT) %>%
  filter(!is.na(AFRLT))  %>%
  # Make all numeric non-integer columns rounded to 3 decimals
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Check dose times are ok
new_adnca  %>%
  mutate(
    actual.dose.time = AFRLT - ARRLT,
    nominal.dose.time = NFRLT - NRRLT
  ) %>%
  group_by(USUBJID, ATPTREF) %>%
  summarise(
    sd_actual = sd(actual.dose.time),
    sd_nominal = sd(nominal.dose.time),
    mean_actual = mean(actual.dose.time)
  ) %>%
  ungroup() %>%
  arrange(ATPTREF, -sd_actual, -sd_nominal)

write.csv(new_adnca, "inst/shiny/data/example-ADNCA.csv", row.names = FALSE)




create_subjects <- function(
    n = 1,
    sex_male_prob = 0.5,
    age_range = c(18, 65),
    races = list("WHITE" = 0.33, "BLACK OR AFRICAN AMERICAN" = 0.33, "ASIAN" = 0.33)
) {
    subjects <- list()
    for (i in 1:n) {
        subject_id <- sprintf("%03d", i)
        sex <- ifelse(runif(1) < sex_male_prob, "M", "F")
        age <- sample(age_range[1]:age_range[2], 1)
        race <- sample(names(races), 1, prob = unlist(races))
        subjects[[subject_id]] <- list(
            SEX = sex,
            AGE = age,
            RACE = race
        )
    }
}

create_doses <- function(
    n = 1,
    dose = 5,
    dose_unit = "mg",
    dose_treatment = "DrugA",
    route = "INTRAVENOUS DRIP",
    start_nominal_dtc = "2024-01-01T08:00:00",
    start_actual_dtc_sd = 0.1,
    dose_nominal_duration = 3,
    dose_actual_duration_sd = 0.05,
    tau = 24,
    tau_actual_sd = 0.1
) {
    doses <- list()
    for (i in 1:n) {
        dose_number <- i
        nominal_start_dtc <- as.POSIXct(start_nominal_dtc) + (dose_number - 1) * tau * 3600
        actual_start_dtc <- nominal_start_dtc + rnorm(1, mean = 0, sd = start_actual_dtc_sd * 3600)
        nominal_end_dtc <- actual_start_dtc + dose_nominal_duration * 3600
        actual_end_dtc <- actual_start_dtc + (dose_nominal_duration + rnorm(1, mean = 0, sd = dose_actual_duration_sd)) * 3600
        doses[[dose_number]] <- list(
            DOSE = dose,
            DOSEU = dose_unit,
            DOSEDUR = dose_nominal_duration,
            DOSETRT = dose_treatment,
            ROUTE = route,
            STDTC = format(nominal_start_dtc, "%Y-%m-%dT%H:%M:%S"),
            ENDTC = format(nominal_end_dtc, "%Y-%m-%dT%H:%M:%S"),
            ACTUAL_STDTM = format(actual_start_dtc, "%Y-%m-%dT%H:%M:%S"),
            ACTUAL_ENDTM = format(actual_end_dtc, "%Y-%m-%dT%H:%M:%S")
        )
    }
    return(doses)
}
    

create_samples <- function(
  nominal_start_times = c(-0.1, 0.5, 1.5, 3, 4.5, 6, 7.5, 12),
  actual_time_sd = 0.05,
  nominal_duration = 0,
  actual_duration_sd = 0,
  pcspec = "SERUM",
  analyte = "DrugA",
  is_metabolite = FALSE,
  # PK characteristics to determine concentration
  cmax = 8,
  lambda_z = 0.5,
  tmax = 3
){
  samples <- list()
  for (i in seq_along(nominal_start_times)) {
      nominal_time <- nominal_start_times[i]
      actual_time <- nominal_time + rnorm(1, mean = 0, sd = actual_time_sd)
      nominal_end_time <- nominal_time + nominal_duration
      actual_end_time <- actual_time + nominal_duration + rnorm(1, mean = 0, sd = actual_duration_sd)
      concentration <- pk_aval(actual_time, Tmax = tmax, Cmax = cmax, lambda_z = lambda_z)
      if (is_metabolite) {
          concentration <- concentration * 0.4
      }
      samples[[i]] <- list(
          NOMINAL_TIME = nominal_time,
          ACTUAL_TIME = actual_time,
          NOMINAL_END_TIME = nominal_end_time,
          ACTUAL_END_TIME = actual_end_time,
          PCSPEC = pcspec,
          PARAM = analyte,
          AVAL = round(concentration + rnorm(1, mean = 0, sd = 0.05), 3),
          AVALU = "ug/mL"
      )
  }
  return(samples)
}

    

# Create a hierarchical data structure for a study with subjects, doses, samples, and analytes
study_data <- list(
    STUDYID = "S1",
    treatment_arms = list()
)

# TRT01A = "DrugA 5 mg, Infusion".
# Population studied involved 6 WHITE (3 female, 3 male) and 6 ASIAN (3 female, 3 male) subjects across different ages. 
# Dose infusions (1h) were given every 24 hours, to WHITE and ASIAN subjects. No missing dose and no protocol deviations in dosing.
# Samples were taken at -0.1, 1.5, 3, 4.5, 6, 7.5 and 12 hours post-dose. No missing sampling, but some protocol deviations.
# Analyte Drug A was analyzed for all samples. Concentrations are additively higher in ASIAN subjects (1.1) and in women (1.15) 
# Analyte Metab-DrugA was analyzed for all samples. Concentrations are 40% of Drug A concentrations. No additional variability.

trt01a <- list()
combn_subjects <- expand.grid(c("WHITE", "ASIAN"), c("M", "F"), stringsAsFactors = FALSE)
for (i in 1:nrow(combn_subjects)) {
  trt01a[(i-1)*3 + 1] <- create_subjects(
    n = 3,
    sex_male_prob = ifelse(combn_subjects[i, 2] == "M", 1, 0),
    age_range = c(18, 65),
    races = list(combn_subjects[i, 1] = 1)
  )
  for (j in 1:3) {
    trt01a[[ (i-1)*3 + j ]]$doses <- create_doses(
      n = 8,
      dose = 5,
      dose_unit = "mg",
      dose_treatment = "DrugA",
      route = "INTRAVENOUS DRIP",
      start_nominal_dtc = "2024-01-01T08:00:00",
      start_actual_dtc_sd = 0.1,
      dose_nominal_duration = 1,
      dose_actual_duration_sd = 0.05,
      tau = 24,
      tau_actual_sd = 0.1
    )
  }
}

study_data$treatment_arms[["DrugA 5 mg, Infusion"]] <- list(
    TRT01A = "DrugA 5 mg, Infusion",
    subjects = list()
)



# Study > Treatment Arm > Subject > Dose > Sample > Analyte

study_data <- list()



data_subj <- list(
    STUDYID = "S1",
    treatment_arms = list(
        "DrugA 5 mg, Infusion" = list(
              subjects = list(
                  "001" = list(
              doses = list(
                  # Dose 1
                  list(
                      DOSE = 5,
                      DOSEU = "mg",
                      DOSEDUR = 3,
                      DOSETRT = "DrugA",
                      ROUTE = "INTRAVENOUS DRIP",
                      STDTC = "2024-01-01T08:00:00",
                      ENDTC = "2024-01-01T10:00:00",
                      samples = list(
                          list(
                              ACTUAL_TIME = -0.1,
                              NOMINAL_TIME = -0.1,
                              PCSPEC = "SERUM",
                              analytes = list(
                                  list(
                                      PARAM = "DrugA"
                                  ),
                                  list(
                                      PARAM = "Metab-DrugA"
                                  )
                              )
                          ),
                      )
                  ),
                )
              )
            )
        )
    )
)
