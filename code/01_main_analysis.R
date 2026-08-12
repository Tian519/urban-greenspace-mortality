# ============================================================
# Main analysis
# Urban greenspace exposure / inequality and mortality
# ============================================================

rm(list = ls())
gc()

# Packages ---------------------------------------------------
library(dplyr)
library(purrr)
library(readxl)
library(glmmTMB)
library(broom.mixed)
library(splines)


# ============================================================
# 1. File paths
# ============================================================

# Use relative paths for reproducibility.
# The mortality dataset is not publicly available due to
# data-sharing restrictions.

data_path  <- "data/demo_data.xlsx"
output_dir <- "results"

# Create results folder if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 2. Read and prepare data
# ============================================================

global_alldata <- read_xlsx(
  data_path,
  sheet = 2
)

global_alldata <- global_alldata %>%
  mutate(
    cityID = as.factor(cityID),
    tmean  = as.numeric(tmean),
    all    = as.numeric(all),
    cvd    = as.numeric(cvd),
    resp   = as.numeric(resp)
  ) %>%
  # Restrict to locations with >=50% urban population
  filter(Per_URBANPOP >= 50)


# ============================================================
# 3. Function for country-specific models
# ============================================================

run_country_models <- function(data,
                               outcome,
                               exposure,
                               output_name) {
  
  # ----------------------------------------------------------
  # Prepare outcome-specific dataset
  # ----------------------------------------------------------
  
  analysis_data <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      .data[[outcome]] > 0
    )
  
  # Number of cities included in each country
  city_number <- analysis_data %>%
    group_by(COUNTRY) %>%
    summarise(
      n_CITIES = n_distinct(cityname),
      .groups = "drop"
    )
  
  # Split data by country
  country_data <- split(
    analysis_data,
    analysis_data$COUNTRY
  )
  
  
  # ----------------------------------------------------------
  # Fit country-specific negative binomial models
  # ----------------------------------------------------------
  
  model_list <- lapply(
    country_data,
    function(df) {
      
      tryCatch({
        
        # Country-specific degrees of freedom for temporal trend
        df_value <- unique(na.omit(df$df))
        
        if (length(df_value) != 1) {
          stop("Expected exactly one temporal-trend df value per country.")
        }
        
        # Model formula
        formula_str <- paste0(
          outcome,
          " ~ offset(log(POP)) + ",
          exposure,
          " + ns(tmean, df = 2)",
          " + PM25",
          " + log(GDP)",
          " + per_older",
          " + ns(year_id, df = ",
          df_value,
          ")"
        )
        
        dynamic_formula <- as.formula(formula_str)
        
        glmmTMB(
          formula = dynamic_formula,
          family = nbinom2,
          data = df
        )
        
      }, error = function(e) {
        
        message(
          "Error in ",
          unique(df$COUNTRY),
          " [", outcome, ", ", exposure, "]: ",
          e$message
        )
        
        return(NULL)
      })
    }
  )
  
  
  # ----------------------------------------------------------
  # Extract coefficients and calculate RR per 0.1-unit increase
  # ----------------------------------------------------------
  
  rr_table <- imap_dfr(
    model_list,
    function(model, country_name) {
      
      if (is.null(model)) {
        
        return(
          tibble(
            COUNTRY = country_name,
            est = NA_real_,
            se = NA_real_,
            RR = NA_real_,
            RR_low = NA_real_,
            RR_high = NA_real_,
            p.value = NA_real_,
            AIC = NA_real_
          )
        )
      }
      
      model_result <- broom.mixed::tidy(model)
      
      exposure_result <- model_result %>%
        filter(term == exposure)
      
      if (nrow(exposure_result) != 1) {
        
        return(
          tibble(
            COUNTRY = country_name,
            est = NA_real_,
            se = NA_real_,
            RR = NA_real_,
            RR_low = NA_real_,
            RR_high = NA_real_,
            p.value = NA_real_,
            AIC = AIC(model)
          )
        )
      }
      
      est <- exposure_result$estimate
      se  <- exposure_result$std.error
      
      tibble(
        COUNTRY = country_name,
        est = est,
        se = se,
        
        # RR per 0.1-unit increase
        RR = exp(est * 0.1),
        
        RR_low = exp(
          (est - 1.96 * se) * 0.1
        ),
        
        RR_high = exp(
          (est + 1.96 * se) * 0.1
        ),
        
        p.value = exposure_result$p.value,
        AIC = AIC(model)
      )
    }
  )
  
  
  # Add number of cities
  rr_table <- rr_table %>%
    left_join(
      city_number,
      by = "COUNTRY"
    )
  
  
  # ----------------------------------------------------------
  # Save results
  # ----------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(output_name, "_", outcome, ".csv")
  )
  
  write.csv(
    rr_table,
    output_file,
    row.names = FALSE
  )
  
  return(rr_table)
}


# ============================================================
# 4. Run analyses for all mortality outcomes
# ============================================================

outcomes <- c(
  "all",
  "cvd",
  "resp"
)


for (outcome in outcomes) {
  
  # ----------------------------------------------------------
  # Greenspace exposure
  # ----------------------------------------------------------
  
  run_country_models(
    data = global_alldata,
    outcome = outcome,
    exposure = "weighted_evi",
    output_name = "EVI_RR"
  )
  
  
  # ----------------------------------------------------------
  # Greenspace inequality
  # ----------------------------------------------------------
  
  run_country_models(
    data = global_alldata,
    outcome = outcome,
    exposure = "EVI_gini",
    output_name = "EVIgini_RR"
  )
}