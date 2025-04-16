# Load necessary libraries
library(rms)
library(survival)
library(timeROC)

# Set seed for reproducibility
set.seed(123)

# Generate some example data
n <- 200
age <- rnorm(n, mean = 60, sd = 10)
gender <- sample(c("Male", "Female"), n, replace = TRUE)
treatment <- sample(c("Treatment A", "Treatment B"), n, replace = TRUE)
time <- rexp(n, rate = 0.1)
status <- sample(c(0, 1), n, replace = TRUE)

# Create a data frame
data <- data.frame(age, gender, treatment, time, status)

# Convert categorical variables to factors
data$gender <- as.factor(data$gender)
data$treatment <- as.factor(data$treatment)

# Show the first few rows of the data
head(data)


# Specify the data distribution
dd <- datadist(data)
options(datadist = 'dd')

# Fit a Cox proportional hazards model
cox_model <- cph(Surv(time, status) ~ age + gender + treatment, data = data, x = TRUE, y = TRUE, surv = TRUE)


# Create a nomogram based on the Cox model
nomogram <- nomogram(cox_model, fun=list(function(x) Survival(cox_model, x, time=12)),
                     lp=TRUE, fun.at=c(0.8, 0.9), funlabel="12-month Survival Probability")

nomogram <- nomogram(cox_model, fun = list(surv1 = function(x) Survival(cox_fit)(12, x)),
                     lp=TRUE,
                     funlabel = c("1-year surv prob"))

# Plot the nomogram
plot(nomogram) 


# Calculate the linear predictor from the Cox model
lp <- predict(cox_model, type = "lp")

# Generate the time-dependent ROC curves for 12 and 24 months
roc_12 <- timeROC(T = data$time, delta = data$status, marker = lp, cause = 1, times = 12, iid = TRUE)
roc_24 <- timeROC(T = data$time, delta = data$status, marker = lp, cause = 1, times = 24, iid = TRUE)

# Plot the ROC curves
plot(roc_12, time = 12, col = "blue", title = "Time-dependent ROC at 12 months")
plot(roc_24, time = 24, col = "red", add = TRUE, lty = 2)
legend("bottomright", legend = c("12 months", "24 months"), col = c("blue", "red"), lty = c(1, 2))



# ------ OUR EXAMPLE ----- #

# --- Proteomics --- #
dat <- df_pat_sum2 %>% 
  rename(npx_bin = RS_binr)

dat <- dat %>% mutate(npx_bin = factor(npx_bin, levels = c("low", "high"))) %>% 
  drop_na()

dd <- datadist(dat)
options(datadist = "dd") 

data <- dat





time_point <- 12

# Fit a Cox proportional hazards model
cox_model_1 <- cph(Surv(OS, death) ~ psa + ldh + alk_ph + alb, data = dat, x = TRUE, y = TRUE, surv = TRUE)
cox_model_2 <- cph(Surv(OS, death) ~ psa + ldh + alk_ph + alb + npx_bin, data = dat, x = TRUE, y = TRUE, surv = TRUE)

# Calculate the linear predictor from the Cox model
lp_1 <- predict(cox_model_1, type = "lp")
lp_2 <- predict(cox_model_2, type = "lp")

# Generate the time-dependent ROC curves for 12 and 24 months
roc_12_1 <- timeROC(T = data$OS, delta = data$death, marker = lp_1, cause = 1, times = time_point, iid = TRUE)
roc_12_2 <- timeROC(T = data$OS, delta = data$death, marker = lp_2, cause = 1, times = time_point, iid = TRUE)

au_1 <- round(roc_12_1$AUC[2],2)
au_2 <- round(roc_12_2$AUC[2],2)
# Plot the ROC curves
pdf(file.path(wd$outCurr, paste0("AUC_Prost_", time_point, ".pdf")),width = 6, height = 6)
plot(roc_12_1, time = time_point, col = "#FFCEAC", lwd=3, title="")
title(paste0(time_point, " months"))
plot(roc_12_2, time = time_point, col = "#BFDBAB", lwd=3, add = TRUE)
legend("bottomright", 
       legend = c(paste0("Clin_", au_1),
                  paste0("Clin + Proteome_", au_2)), 
       col = c("#FFCEAC", "#BFDBAB"), lty = c(1, 1))
dev.off()



# ---- -Methylation ---- #
dat_l <- read_xlsx("../ZW_025_DoD_grant/Book1.xlsx") %>% 
  mutate(death = case_when(death_y_n == "N" ~ 0,
                           death_y_n == "Y" ~ 1),
         OS = s1_to_death_or_last_fu_mo,
         alk_ph = as.numeric(alk_phos_at_s1),
         #alb = Sample1_ALB,
         #ldh = as.numeric(ldh_at_s1),
         HCI_cID = Samples
  ) %>%  
  rename(npx_bin = 'Methylome Risk') %>% 
  mutate(npx_bin = factor(npx_bin, levels=c("High", "Low")))


# Liang is missing PSA, so we add it manually
df_hci <- df_meta_f %>%  
  filter(!HCI_cID %in% serial_sam)

d1 <- df_hci %>% filter(isPsom == "Y") %>% mutate(sample_id = paste0("Ps_", sample_id_psom)) %>% 
  select(HCI_cID, sample_id) %>% distinct()
d2 <- df_hci %>% filter(isQ13356 == "Y") %>% mutate(sample_id = paste0("Q1335Q_", sample_id_Q13356)) %>% 
  select(HCI_cID, sample_id) %>% distinct()
d3 <- df_hci %>% filter(isQ15806 == "Y") %>% mutate(sample_id = paste0("Q1580Q_", sample_id_Q15806)) %>% 
  select(HCI_cID, sample_id) %>% distinct()

d_all <- rbind(d1,d2,d3)

dat_psa <- dat %>% select(sample_id, death,psa:alb)
d_all <- d_all %>% left_join(dat_psa)

dat_l_sub <- dat_l %>% select(HCI_cID, npx_bin)

dat_l2 <- dat_l_sub %>% left_join(d_all) %>% 
  filter(!is.na(death)) %>% 
  select(-sample_id) %>% distinct()


dat <- dat_l2 %>% 
  drop_na()

dd <- datadist(dat)
options(datadist = "dd") 

data <- dat

time_point <- 36

# Fit a Cox proportional hazards model
cox_model_1 <- cph(Surv(OS, death) ~ psa + ldh + alk_ph + alb, data = dat, x = TRUE, y = TRUE, surv = TRUE)
cox_model_2 <- cph(Surv(OS, death) ~ psa + ldh + alk_ph + alb + npx_bin, data = dat, x = TRUE, y = TRUE, surv = TRUE)

# Calculate the linear predictor from the Cox model
lp_1 <- predict(cox_model_1, type = "lp")
lp_2 <- predict(cox_model_2, type = "lp")

# Generate the time-dependent ROC curves for 12 and 24 months
roc_12_1 <- timeROC(T = data$OS, delta = data$death, marker = lp_1, cause = 1, times = time_point, iid = TRUE)
roc_12_2 <- timeROC(T = data$OS, delta = data$death, marker = lp_2, cause = 1, times = time_point, iid = TRUE)

au_1 <- round(roc_12_1$AUC[2],2)
au_2 <- round(roc_12_2$AUC[2],2)
# Plot the ROC curves
pdf(file.path(wd$outCurr, paste0("AUC_Methy_", time_point, ".pdf")),width = 6, height = 6)
plot(roc_12_1, time = time_point, col = "#FFCEAC", lwd=3, title="")
title(paste0(time_point, " months"))
plot(roc_12_2, time = time_point, col = "#A7BFD0", lwd=3, add = TRUE)
legend("bottomright", 
       legend = c(paste0("Clin_", au_1),
                  paste0("Clin + Methy_", au_2)), 
       col = c("#FFCEAC", "#A7BFD0"), lty = c(1, 1))
dev.off()







# If we jus tread liangs data without using updated clinical data
dat <- read_xlsx("../ZW_025_DoD_grant/Book1.xlsx") %>% 
  mutate(death = case_when(death_y_n == "N" ~ 0,
                           death_y_n == "Y" ~ 1),
         OS = s1_to_death_or_last_fu_mo,
         alk_ph = as.numeric(alk_phos_at_s1),
         #alb = Sample1_ALB,
         ldh = as.numeric(ldh_at_s1),
         HCI_cID = Samples
  ) %>%  
  rename(npx_bin = 'Methylome Risk') %>% 
  mutate(npx_bin = factor(npx_bin, levels=c("High", "Low")))

data <- dat %>% 
  drop_na()

dd <- datadist(data)
options(datadist = "dd") 

time_point <- 12

# Fit a Cox proportional hazards model
cox_model_1 <- cph(Surv(OS, death) ~   ldh + alk_ph, data = data, x = TRUE, y = TRUE, surv = TRUE)
cox_model_2 <- cph(Surv(OS, death) ~   ldh + alk_ph + npx_bin, data = data, x = TRUE, y = TRUE, surv = TRUE)

# Calculate the linear predictor from the Cox model
lp_1 <- predict(cox_model_1, type = "lp")
lp_2 <- predict(cox_model_2, type = "lp")

# Generate the time-dependent ROC curves for 12 and 24 months
roc_12_1 <- timeROC(T = data$OS, delta = data$death, marker = lp_1, cause = 1, times = time_point, iid = TRUE)
roc_12_2 <- timeROC(T = data$OS, delta = data$death, marker = lp_2, cause = 1, times = time_point, iid = TRUE)

au_1 <- round(roc_12_1$AUC[2],2)
au_2 <- round(roc_12_2$AUC[2],2)
# Plot the ROC curves
pdf(file.path(wd$outCurr, paste0("AUC_Meth_", time_point, ".pdf")),width = 6, height = 6)
plot(roc_12_1, time = time_point, col = "#FFCEAC", lwd=3, title="")
title(paste0(time_point, " months"))
plot(roc_12_2, time = time_point, col = "#BFDBAB", lwd=3, add = TRUE)
legend("bottomright", 
       legend = c(paste0("Clin_", au_1),
                  paste0("Clin + Methy", au_2)), 
       col = c("#FFCEAC", "#BFDBAB"), lty = c(1, 1))
dev.off()

