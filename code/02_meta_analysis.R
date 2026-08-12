rm(list=ls());
gc() #清理内存

library(readxl)
library(dplyr)
library(purrr)
library(mixmeta)
library(ggplot2)
library(writexl)



# ------------------------------------------------------------
# 1. File and sheets
# ------------------------------------------------------------

file <- "~data/RR.xlsx"

sheets <- c(
  "evi_all",
  "evigini_all",
  "evi_cvd",
  "evigini_cvd",
  "evi_resp",
  "evigini_resp"
)


# ------------------------------------------------------------
# 2. Function: random-effects pooling
# ------------------------------------------------------------

pool_meta <- function(dat) {
  
  dat <- dat %>%
    filter(
      !is.na(est),
      !is.na(se),
      is.finite(est),
      is.finite(se),
      se > 0
    )
  
  # Random-effects meta-analysis
  # est = log(RR) coefficient per 1-unit increase
  # S = within-country variance
  fit <- mixmeta(
    est,
    S = se^2,
    data = dat,
    method = "reml"
  )
  
  sm <- summary(fit)
  
  beta <- as.numeric(coef(fit)[1])
  beta_se <- sqrt(as.numeric(vcov(fit)[1, 1]))
  
  beta_low  <- beta - 1.96 * beta_se
  beta_high <- beta + 1.96 * beta_se
  
  p_value <- 2 * pnorm(
    abs(beta / beta_se),
    lower.tail = FALSE
  )
  
  # Transform to RR per 0.1-unit increase
  RR      <- exp(beta * 0.1)
  RR_low  <- exp(beta_low * 0.1)
  RR_high <- exp(beta_high * 0.1)
  
  # Between-country variance
  tau2 <- as.numeric(fit$Psi)[1]
  
  # I-squared
  I2 <- as.numeric(sm$i2stat)[1]
  
  tibble(
    beta = beta,
    beta_se = beta_se,
    RR = RR,
    RR_low = RR_low,
    RR_high = RR_high,
    p_value = p_value,
    tau2 = tau2,
    I2 = I2,
    n_country = nrow(dat),
    n_cities = sum(dat$n_CITIES, na.rm = TRUE)
  )
}


# ------------------------------------------------------------
# 3. Function: leave one country out
# ------------------------------------------------------------

leave_one_country_out <- function(dat, sheet_name) {
  
  # Make sure country is character
  dat <- dat %>%
    mutate(COUNTRY = as.character(COUNTRY))
  
  
  # -------------------------
  # Full pooled model
  # -------------------------
  
  full <- pool_meta(dat) %>%
    mutate(
      sheet = sheet_name,
      omitted_country = "None (full model)",
      omitted_n_CITIES = 0
    )
  
  
  # -------------------------
  # Leave-one-country-out
  # -------------------------
  
  loo <- map_dfr(
    unique(dat$COUNTRY),
    
    function(country_i) {
      
      dat_sub <- dat %>%
        filter(COUNTRY != country_i)
      
      omitted_n <- dat %>%
        filter(COUNTRY == country_i) %>%
        pull(n_CITIES)
      
      omitted_n <- sum(omitted_n, na.rm = TRUE)
      
      result <- pool_meta(dat_sub)
      
      result %>%
        mutate(
          sheet = sheet_name,
          omitted_country = country_i,
          omitted_n_CITIES = omitted_n
        )
    }
  )
  
  
  # -------------------------
  # Compare with full model
  # -------------------------
  
  loo <- loo %>%
    mutate(
      full_RR = full$RR,
      
      # Absolute difference from original pooled RR
      RR_difference = RR - full$RR,
      
      # Relative percentage change
      percent_change =
        100 * (RR / full$RR - 1),
      
      abs_percent_change =
        abs(percent_change),
      
      # Does the leave-one-out CI contain null?
      CI_crosses_1 =
        RR_low <= 1 & RR_high >= 1
    ) %>%
    arrange(desc(omitted_n_CITIES))
  
  
  list(
    full = full,
    loo = loo
  )
}


# ------------------------------------------------------------
# 4. Run all six sheets
# ------------------------------------------------------------

results <- map(
  sheets,
  function(sheet_i) {
    
    dat <- read_excel(
      file,
      sheet = sheet_i
    )
    
    leave_one_country_out(
      dat = dat,
      sheet_name = sheet_i
    )
  }
)

names(results) <- sheets


# ------------------------------------------------------------
# 5. Combine full pooled results
# ------------------------------------------------------------

full_results <- map_dfr(
  results,
  "full"
) %>%
  select(
    sheet,
    RR,
    RR_low,
    RR_high,
    p_value,
    tau2,
    I2,
    n_country,
    n_cities
  )

print(full_results)


# ------------------------------------------------------------
# 6. Combine all leave-one-out results
# ------------------------------------------------------------

loo_results <- map_dfr(
  results,
  "loo"
) %>%
  select(
    sheet,
    omitted_country,
    omitted_n_CITIES,
    RR,
    RR_low,
    RR_high,
    p_value,
    tau2,
    I2,
    n_country,
    n_cities,
    full_RR,
    RR_difference,
    percent_change,
    abs_percent_change,
    CI_crosses_1
  )

print(loo_results)




# ------------------------------------------------------------
# 8. Countries with greatest influence on pooled RR
# based on absolute percentage change
# ------------------------------------------------------------

most_influential <- loo_results %>%
  group_by(sheet) %>%
  arrange(desc(abs_percent_change), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

print(most_influential)


# ------------------------------------------------------------
# 9. Summary of leave-one-out RR range
# Useful for manuscript/reviewer response
# ------------------------------------------------------------

loo_summary <- loo_results %>%
  group_by(sheet) %>%
  summarise(
    full_RR = first(full_RR),
    
    min_LOO_RR = min(RR, na.rm = TRUE),
    max_LOO_RR = max(RR, na.rm = TRUE),
    
    max_abs_percent_change =
      max(abs_percent_change, na.rm = TRUE),
    
    most_influential_country =
      omitted_country[
        which.max(abs_percent_change)
      ],
    
    .groups = "drop"
  )

print(loo_summary)


# ------------------------------------------------------------
# 10. Export results to Excel
# ------------------------------------------------------------

output_list <- list(
  
  full_pooled = full_results,
  
  all_leave_one_out = loo_results,
  
  largest_contributors = largest_contributors,
  
  most_influential = most_influential,
  
  loo_summary = loo_summary
)


# Also save each analysis separately
for (sheet_i in sheets) {
  
  output_list[[paste0("LOO_", sheet_i)]] <-
    results[[sheet_i]]$loo
}


write_xlsx(
  output_list,
  "pooled_results.xlsx"
)




