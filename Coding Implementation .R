## See Andreas Bender Survival Analysis Model Evaluation
library(dplyr)
library(mgcv)
library(pammtools)
library(ggplot2)
theme_set(theme_bw())
library(survival)
library(pec)

patient <- patient
tumor <- tumor


## split data into train and test data
n_train   <- 400
train_idx <- sample(seq_len(nrow(tumor)), n_train)
test_idx  <- setdiff(seq_len(nrow(tumor)), train_idx)
## data transformation
tumor_ped <- as_ped(tumor[train_idx, ], Surv(days, status)~.)

# some simple models for comparison
pam1 <- pamm(
  formula = ped_status ~ s(tend) + charlson_score + age,
  data = tumor_ped)
pam2 <- pamm(
  formula = ped_status ~ s(tend) + charlson_score + age + metastases + complications,
  data = tumor_ped)
pam3 <- pamm(
  formula = ped_status ~s(tend, by = complications) + charlson_score + age +
    metastases,
  data = tumor_ped)
# calculate prediction error curves (on test data)
pec <- pec(
  list(pam1 = pam1, pam2 = pam2, pam3 = pam3),
  Surv(days, status) ~ 1, # formula for IPCW
  data = tumor[test_idx, ], # new data not used for model fit
  times = seq(.01, 1200, by = 10),
  start = .01,
  exact = FALSE
)

plot(pec)

```{r}
plot(cindex(
  list(pam1 = pam1, pam2 = pam2, pam3 = pam3),
  Surv(days, status) ~ 1,
  data = tumor[test_idx, ],
  eval.times = quantile(tumor$days[tumor$status == 1], c(.25, .5, .75))))
```


```{r}

# For further information: https://mlr3book.mlr-org.com/survival.html
library("mlr3learners")
library("mlr3")
library("mlr3proba")
library("survival")
library("ranger")
library("mlr3viz")


TaskSurv$new(id = "interval_censored", backend = survival::bladder2[,-c(1, 7)],
             time = "start", time2 = "stop", type = "interval2")

task = tsk("rats")

# some integrated learners
learners = lrns(c("surv.coxph", "surv.kaplan", "surv.ranger"))
print(learners)

??msr()
# Harrell's C-Index for survival
measure = msr("surv.cindex")
print(measure)

set.seed(1)
bmr = benchmark(benchmark_grid(task, learners, rsmp("cv", folds = 3)))
bmr$aggregate(measure)
autoplot(bmr, measure = measure)

```


```{r}
# for reference same with brier score
# Plot after
measure = msr("surv.brier")
print(measure)

set.seed(1)
bmr = benchmark(benchmark_grid(task, learners, rsmp("cv", folds = 3)))
bmr$aggregate(measure)
autoplot(bmr, measure = measure)

```
