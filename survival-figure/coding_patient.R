patient <-patient

## split data into train and test data
n_train   <- 1400
train_idx <- sample(seq_len(nrow(patient)), n_train)
test_idx  <- setdiff(seq_len(nrow(patient)), train_idx)
## data transformation
patient_ped <- as_ped(patient[train_idx, ], Surv(survhosp, PatientDied)~.)

# some simple models for comparison

pam1 <- pamm( formula = Survdays ~ s(tend)  + Age + BMI, data = patient_ped)
pam2 <- pamm(formula = Survdays ~ s(tend)  + Age + BMI + Gender, data = patient_ped)
pam3 <- pamm( formula = Survdays ~s(tend, by = Gender) + BMI, data = patient_ped)

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
