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



create_subjects <- function(
    n = 1,
    sex = "M",
    race = "WHITE",
    age = 21,
    subjid = NULL
) {
    subjects <- list()
    for (i in 1:n) {
        subject_id <- ifelse(is.null(subjid), sprintf("S%03d", i), subjid[i])
        subjects[[i]] <- list(
            USUBJID = if (length(subjid) == 1) subject_id else subjid[i],
            SEX = if (length(sex) == 1) sex else sex[i],
            RACE = if (length(race) == 1) race else race[i],
            AGE = if (length(age) == 1) age else age[i],
            AGEU = "Years"
        )
    }
    return(subjects)
}

create_doses <- function(
    n = 1,
    dose = 5,
    dose_unit = "mg",
    dose_treatment = "DrugA",
    route = "INTRAVENOUS DRIP",
    start_nominal_dtc = "2024-01-01 08:00:00",
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
    

# Function to generate sample time points (without analyte info)
create_samples <- function(
  nominal_start_times = c(-0.1, 0.5, 1.5, 3, 4.5, 6, 7.5, 12),
  actual_time_sd = 0.05,
  nominal_duration = 0,
  actual_duration_sd = 0,
  pcspec = "SERUM",
  volume = 0,
  volume_unit = "mL"
){
  samples <- list()
  for (i in seq_along(nominal_start_times)) {
    nominal_time <- nominal_start_times[i]
    actual_time <- nominal_time + rnorm(1, mean = 0, sd = actual_time_sd)
    nominal_end_time <- nominal_time + nominal_duration
    actual_end_time <- actual_time + nominal_duration + rnorm(1, mean = 0, sd = actual_duration_sd)
    samples[[i]] <- list(
      NOMINAL_TIME = nominal_time,
      ACTUAL_TIME = actual_time,
      NOMINAL_END_TIME = nominal_end_time,
      ACTUAL_END_TIME = actual_end_time,
      PCSPEC = pcspec,
      VOLUME = volume,
      VOLUMEU = volume_unit
    )
  }
  return(samples)
}

# Function to add analyte information and calculate concentration for a sample
create_analyte <- function(
  actual_time,
  analyte = "DrugA",
  is_metabolite = FALSE,
  cmax = 8,
  lambda_z = 0.5,
  tmax = 3,
  sd = 0.05,
  lloq = 0.05,
  avalu = "ug/mL"
){
  # Calculate concentration based on PK characteristics
  concentration <- pk_aval(actual_time, Tmax = tmax, Cmax = cmax, lambda_z = lambda_z)
  concentration <- concentration + rnorm(1, mean = 0, sd = sd)
  concentration <- ifelse(concentration <= lloq, 0, concentration)
  list(
    PARAM = analyte,
    AVAL = round(concentration, 3),
    AVALU = avalu,
    METABFL = ifelse(is_metabolite, "Y", "")
  )
}





# Create a function that generates the whole treatment arm
create_treatment_arm_epoch <- function(
  trt_epoch = "Treatment",
  # Population characteristics
  n_subjects = 12,
  subjids = paste0("S1-", sprintf("%03d", seq_along(n_subjects))),
  race_sex_proportion = list(WHITE = list(M = 0.5, F = 0.5), ASIAN = list(M = 0.5, F = 0.5)),
  age_range = c(18, 65),

  # Dosing regimen
  dose = 5,
  dose_unit = "mg",
  dose_treatment = "DrugA",
  route = "INTRAVENOUS DRIP",
  start_nominal_dtc = "2024-01-01 08:00:00",
  start_actual_dtc_sd = 0.1,
  dose_nominal_duration = 3,
  dose_actual_duration_sd = 0.05,
  tau = 24,
  tau_actual_sd = 0.1,

  # Sampling scheme
  pcspec = list("SERUM", "URINE"),
  nominal_start_times = list(
    SERUM = c(-0.1, 0.5, 1.5, 3, 4.5, 6, 7.5, 12),
    URINE = c(12, 24, 36)
  ),
  actual_time_sd = list(SERUM = 0.05, URINE = 0.1),
  nominal_duration = list(SERUM = 0, URINE = 6),
  actual_duration_sd = list(SERUM = 0, URINE = 0),
  volume = list(SERUM = 0, URINE = 0.6),
  volume_unit = "mL",

  # Analyte characteristics
  analytes = list(
    SERUM = list(
      list(
        analyte = "DrugA",
        is_metabolite = FALSE,
        cmax = 8,
        lambda_z = 0.5,
        tmax = 3,
        sd = 0.05,
        fun_effect = pk_aval,
        lloq = 0.05,
        avalu = "ug/mL"
      ),
      list(
        analyte = "Metab-DrugA",
        is_metabolite = TRUE,
        cmax = 3.2,
        lambda_z = 0.5,
        tmax = 3,
        sd = 0.075,
        fun_effect = pk_aval,
        lloq = 0.05,
        avalu = "ug/mL"
      )
    ),
    URINE = list(
      list(
        analyte = "DrugA",
        is_metabolite = FALSE,
        cmax = 80,
        lambda_z = 0.5,
        tmax = 3,
        sd = 0.1,
        fun_effect = pk_aval,
        lloq = 0.05,
        avalu = "mg/mL"
      ),
      list(
        analyte = "Metab-DrugA",
        is_metabolite = TRUE,
        cmax = 32,
        lambda_z = 0.5,
        tmax = 3,
        sd = 0.15,
        fun_effect = pk_aval,
        lloq = 0.05,
        avalu = "mg/mL"
      )
    )
  ),
  # PK variability modifiers (functions)
  sex_conc_effect = list(M = function(cmax) cmax, F = function(cmax) cmax * 1.15),
  sex_sd_effect = list(M = function(sd) sd, F = function(sd) sd),
  race_conc_effect = list(WHITE = function(cmax) cmax, ASIAN = function(cmax) cmax * 1.15),
  race_sd_effect = list(WHITE = function(sd) sd, ASIAN = function(sd) sd * 1.15),
  age_conc_effect = function(age, cmax) cmax,
  age_sd_effect = function(age, sd) sd,
  # Random seed for reproducibility
  seed = 123
) {
  set.seed(seed)
  subjects <- list()

  # Generate subject-level info
  for (i in 1:n_subjects) {
    # Determine race and sex based on specified proportions
    races <- sample(names(race_sex_proportion), 1, prob = sapply(race_sex_proportion, function(x) sum(unlist(x))))
    sexes <- sample(names(race_sex_proportion[[races]]), 1, prob = unlist(race_sex_proportion[[races]]))
    ages <- sample(age_range[1]:age_range[2], 1)
    subject_id <- ifelse(length(subjids) == 1, sprintf("S%03d", i), subjids[i])
    subj_info <- create_subjects(
      n = 1,
      sex = sexes,
      race = races,
      age = ages,
      subjid = subject_id
    )[[1]]

    # Dosing for this subject
    doses <- create_doses(
      n = 8,
      dose = dose,
      dose_unit = dose_unit,
      dose_treatment = dose_treatment,
      route = route,
      start_nominal_dtc = start_nominal_dtc,
      start_actual_dtc_sd = start_actual_dtc_sd,
      dose_nominal_duration = dose_nominal_duration,
      dose_actual_duration_sd = dose_actual_duration_sd,
      tau = tau,
      tau_actual_sd = tau_actual_sd
    )

    # For each dose, generate samples and analytes
    for (d in seq_along(doses)) {
      all_samples <- list()
      for (spec in names(nominal_start_times)) {
        # Generate samples for this specimen type
        samples <- create_samples(
          nominal_start_times = nominal_start_times[[spec]],
          actual_time_sd = actual_time_sd[[spec]],
          nominal_duration = nominal_duration[[spec]],
          actual_duration_sd = actual_duration_sd[[spec]],
          pcspec = spec,
          volume = ifelse(is.list(volume), volume[[spec]], volume),
          volume_unit = volume_unit
        )
        # For each sample, add analytes
        for (s in seq_along(samples)) {
          analyte_list <- list()
          for (a in analytes[[spec]]) {
            # Apply PK variability modifiers
            cmax_mod <- sex_conc_effect[[subj_info$SEX]](
              race_conc_effect[[subj_info$RACE]](
                age_conc_effect(subj_info$AGE, a$cmax)
              )
            )
            sd_mod <- sex_sd_effect[[subj_info$SEX]](
              race_sd_effect[[subj_info$RACE]](
                age_sd_effect(subj_info$AGE, a$sd)
              )
            )
            # Calculate analyte
            analyte_info <- create_analyte(
              actual_time = samples[[s]]$ACTUAL_TIME,
              analyte = a$analyte,
              is_metabolite = a$is_metabolite,
              cmax = cmax_mod,
              lambda_z = a$lambda_z,
              tmax = a$tmax,
              sd = sd_mod,
              lloq = a$lloq,
              avalu = a$avalu
            )
            analyte_list[[length(analyte_list) + 1]] <- analyte_info
          }
          samples[[s]]$analytes <- analyte_list
        }
        all_samples <- c(all_samples, samples)
      }
      doses[[d]]$samples <- all_samples
    }
    subj_info$doses <- doses
    subjects[[i]] <- subj_info
  }
  return(list(
    trt_epoch = trt_epoch,
    subjects = subjects
  ))
}
















    


study_data <- list(
    STUDYID = "S1",
    treatment_arms = list()
)

# The study plans to compare two treatments: TRT01A = "DrugA 5 mg, Infusion", and TRT02A = "DrugA 10 mg, Infusion". For each treatment arm, we will create a list of subjects, and for each subject, we will create a list of doses, and for each dose, we will create a list of samples with analyte concentrations.
for (trt in c("DrugA 5 mg, Infusion", "DrugA 10 mg, Infusion", "DrugA 5 mg, Bolus")) {

  # For each treatment arm, population studied involved a balanced sex, race (white and asian) composition of 12 adult volunteers across different ages.
  combn_subjects <- expand.grid(c("WHITE", "ASIAN"), c("M", "F"), stringsAsFactors = FALSE)
  combn_subjects <- combn_subjects[rep(1:nrow(combn_subjects), each = 3), ]

  study_data$treatment_arms[[trt]] <- list(
    TRT01A = trt,
    
    subjects = create_subjects(
      n = nrow(combn_subjects),
      race = combn_subjects$Var1,
      sex = combn_subjects$Var2,
      age = sample(18:65, nrow(combn_subjects), replace = TRUE),
      subjid = paste0(study_data$STUDYID, sprintf("%03d", 1:nrow(combn_subjects)))
    )

  for (subj in seq_along(study_data$treatment_arms[[trt]]$subjects)) {
    race <- combn_subjects[subj, 1]
    sex  <- combn_subjects[subj, 2]
    age <- sample(18:65, 1)

    study_data$treatment_arms[[trt]]$subjects[[subj]] <- create_subjects(
      n = 1,
      race = race,
      sex = sex,
      age = age,
      subjid = paste0(study_data$STUDYID, "-", sprintf("%03d", subj))
    )[[1]]
    # For each subject, dose infusions (2h) were given every 24 hours during 7 days. No missing dose and no protocol deviations in dosing.
    start_nominal_dtc_subj <- as.character(as.POSIXct("2024-01-01 08:00:00") + (subj - 1) * 3600 * 24 * 14)
    study_data$treatment_arms[[trt]]$subjects[[subj]]$doses <- create_doses(
      n = 8,
      dose = ifelse(trt == "DrugA 5 mg, Infusion", 5, 10),
      dose_unit = "mg",
      dose_treatment = "DrugA",
      route = "INTRAVENOUS DRIP",
      start_nominal_dtc = start_nominal_dtc_subj,
      start_actual_dtc_sd = 0.1,
      dose_nominal_duration = 2,
      dose_actual_duration_sd = 0.05,
      tau = 24,
      tau_actual_sd = 0.1
    )

    for (dose in seq_along(study_data$treatment_arms[[trt]]$subjects[[subj]]$doses)) {
      # For each dose, serum samples were taken at -0.1, 0.5, 1.5, 3, 4.5, 6, 7.5 and 12 hours post-dose. No missing sampling, but some protocol deviations.
      study_data$treatment_arms[[trt]]$subjects[[subj]]$doses[[dose]]$samples <- c(
        # For each dose, serum samples were taken at -0.1, 0.5, 1.5, 3, 4.5, 6, 7.5 and 12 hours post-dose. No missing sampling, but some protocol deviations. 
        create_samples(
          nominal_start_times = c(-0.1, 0.5, 1.5, 3, 4.5, 6, 7.5, 12),
          actual_time_sd = 0.05,
          nominal_duration = 0,
          actual_duration_sd = 0,
          pcspec = "SERUM",
          volume = 0
        ),
        # For each dose, 6h urine samples were taken at 12, 24 and 36 hours post-dose. No missing sampling, but some protocol deviations.
        create_samples(
          nominal_start_times = c(12, 24, 36),
          actual_time_sd = 0.1,
          nominal_duration = 6,
          actual_duration_sd = 0,
          pcspec = "URINE",
          volume = round(rnorm(1, mean = 0.6, sd = 0.1), 2)
        )
      )
      for (sample in seq_along(study_data$treatment_arms[[trt]]$subjects[[subj]]$doses[[dose]]$samples)) {
        pcspec <- study_data$treatment_arms[[trt]]$subjects[[subj]]$doses[[dose]]$samples[[sample]]$PCSPEC
        actual_time <- study_data$treatment_arms[[trt]]$subjects[[subj]]$doses[[dose]]$samples[[sample]]$ACTUAL_TIME
        # Drug concentrations are additively higher in women (1.15).
        cmax <- ifelse(trt == "DrugA 5 mg, Infusion", 8, 16) * (ifelse(sex == "F", 1.15, 1)) * (ifelse(pcspec == "URINE", 10, 1)) 
        # Variability in concentrations is additively higher in ASIAN subjects (1.15).
        sd <- 0.05 * (1 + ifelse(race == "ASIAN", 0.15, 0))
        # Clearance is additively lower in subjects with higher age.
        lz <- 0.5 * (1 + (study_data$treatment_arms[[trt]]$subjects[[subj]]$AGE - 18) / (65 - 18) * 0.5)
        study_data$treatment_arms[[trt]]$subjects[[subj]]$doses[[dose]]$samples[[sample]]$analytes <- c(
          # For each sample, analyte Drug A was analyzed.
          create_analyte(
            actual_time = actual_time,
            analyte = "DrugA",
            is_metabolite = FALSE,
            cmax = cmax,
            lambda_z = lz,
            sd = sd,
            tmax = 3,
            lloq = 0.05,
            avalu = ifelse(pcspec == "SERUM", "ug/mL", "mg/mL")
          ),
          # For each sample, analyte Metab-DrugA was analyzed.
          # Concentrations are 40% of Drug A concentrations. Higher individual variability.
          create_analyte(
            actual_time = actual_time,
            analyte = "Metab-DrugA",
            is_metabolite = TRUE,
            cmax = cmax * 0.4,
            lambda_z = lz,
            sd = sd * 1.5,
            tmax = 3,
            lloq = 0.05,
            avalu = ifelse(pcspec == "SERUM", "ug/mL", "mg/mL")
          )
        )
      }
    }
  }
}



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
