library(dplyr)
library(mgcv)
library(pammtools)
library(ggplot2)
theme_set(theme_bw())
library(survival)
library(pec)

patient <-patient

## split data into train and test data
n_train   <- 1400
train_idx <- sample(seq_len(nrow(patient)), n_train)
test_idx  <- setdiff(seq_len(nrow(patient)), train_idx)
## data transformation
patient_ped <- as_ped(patient[train_idx, ], Surv(survhosp, PatientDied)~.)

# some simple models for comparison

pam1 <- pamm( formula = ped_status ~ s(tend)  + Age + BMI, data = patient_ped)
pam2 <- pamm(formula = ped_status ~ s(tend)  + Age + BMI + Gender, data = patient_ped)
pam3 <- pamm( formula = ped_status ~s(tend, by = Gender) + BMI, data = patient_ped)

# calculate prediction error curves (on test data)

pec <- pec( list(pam1 = pam1, pam2 = pam2, pam3 = pam3),
  Surv(Survdays, PatientDied) ~ 1, # formula for IPCW
  data = patient[test_idx, ], # new data not used for model fit
  times = seq(.01, 1200, by = 10),
  start = .01,
  exact = FALSE)

plot(pec)

crps(pec, times = quantile(patient_ped$Survdays[patient_ped$PatientDied == 1], c(.25, .5, .75)))

```{r}
plot(cindex(
  list(pam1 = pam1, pam2 = pam2, pam3 = pam3),
  Surv(Survdays, PatientDied) ~ 1,
  data = patient[test_idx, ],
  eval.times = quantile(patient$days[patient$PatientDied == 1], c(.25, .5, .75))))
```
summary(patient)
















$$Se^I = P(X_i>c|T_i=t) \tag{3}$$
$$Sp^D = P(X_i\|T_i>t) \tag{4}$$
$$\mathrm{IAUC}=\operatorname{Pr}\left\{z\left(\mathbf{X}_{i}\right)>z\left(\mathbf{X}_{j}\right) \mid D_{i}=1  \&  D_{j}=0\right\} \tag{7} $$
