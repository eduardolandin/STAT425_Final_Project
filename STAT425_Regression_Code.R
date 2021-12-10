library(BAS)
library(devtools)
library(dplyr)
library(ggplot2)
library(grid)
library(gridBase)
library(mvtnorm)
library(pracma)
library(stringr)

## Data joining ----------------------------------------------------------------
data <- read.csv("hierarchical_df.csv")
data <- data %>% mutate(Cases.Per.1k = 1000 * Cases / Population)

# mask data
mask_data <- read.csv("mask_use_by_county.csv")
data_v2 <- merge(
  x = data,
  y = mask_data,
  by.x = "Fips",
  by.y = "COUNTYFP",
  all.x = TRUE,
  all.y = FALSE
)
colnames(data_v2) <- c("Fips", "X", "County", "State", "Date", "Cases", "Deaths", "Population", "Cases.Per.1k", "Mask.Never", "Mask.Rarely", "Mask.Sometimes", "Mask.Frequently", "Mask.Always")

# election data
data_election <- read.csv("county_pres_data_2000_2020.csv")
election_2020 <- data_election %>% filter(year == 2020)
election_2020 <- election_2020 %>% select(-office, - candidate)
election_2020 <- election_2020 %>% filter(mode == "TOTAL")
election_2020 <- election_2020 %>% mutate(Pct.Vote = candidatevotes / totalvotes)

election_2020_totals <- election_2020 %>% filter(party == "DEMOCRAT")
election_2020_totals <- election_2020_totals %>% select(-year, -state, -state_po, -county_name, -party, -candidatevotes, -version, -mode, -Pct.Vote)
colnames(election_2020_totals) <- c("County.Fips", "Total.Votes")

election_2020_democrat <- election_2020 %>% filter(party == "DEMOCRAT")
election_2020_democrat <- election_2020_democrat %>% select(-year, -state, -state_po, -county_name, -party, -candidatevotes, -totalvotes, -version, -mode)
colnames(election_2020_democrat) <- c("County.Fips", "Pct.Democrat")

election_2020_republican <- election_2020 %>% filter(party == "REPUBLICAN")
election_2020_republican <- election_2020_republican %>% select(-year, -state, -state_po, -county_name, -party, -candidatevotes, -totalvotes, -version, -mode)
colnames(election_2020_republican) <- c("County.Fips", "Pct.Republican")

data_v3 <- merge(
  x = data_v2,
  y = election_2020_totals,
  by.x = "Fips",
  by.y = "County.Fips",
  all.x = TRUE
)

data_v3 <- merge(
  x = data_v3,
  y = election_2020_democrat,
  by.x = "Fips",
  by.y = "County.Fips",
  all.x = TRUE
)

data_v3 <- merge(
  x = data_v3,
  y = election_2020_republican,
  by.x = "Fips",
  by.y = "County.Fips",
  all.x = TRUE
)

data_v3 <- data_v3 %>% select(-X)

## Data Cleaning ---------------------------------------------------------------
data <- read.csv("regression_df.csv")
data <- data %>% filter(!is.na(Pct.Democrat))
data <- data %>% mutate(Pct.Other = 1 - Pct.Democrat - Pct.Republican)
data <- data %>% mutate(Votes.Per.Capita = Total.Votes / Population)
data <- data %>% mutate(Mask.Positivity = (Mask.Frequently + Mask.Always) - (Mask.Rarely + Mask.Never))
data <- data %>% mutate(Intercept = rep(1, nrow(data)))

# EDA plots:
plot(x = data$Votes.Per.Capita, y = data$Cases.Per.1k)
plot(x = data$Pct.Republican, y = data$Cases.Per.1k)
plot(x = data$Pct.Democrat, y = data$Cases.Per.1k)
plot(x = data$Population, y = data$Cases.Per.1k)
plot(x = data$Mask.Never, y = data$Cases.Per.1k)
plot(x = data$Mask.Rarely, y = data$Cases.Per.1k)
plot(x = data$Mask.Sometimes, y = data$Cases.Per.1k)
plot(x = data$Mask.Frequently, y = data$Cases.Per.1k)
plot(x = data$Mask.Always, y = data$Cases.Per.1k)
plot(x = data$Pct.Republican, y = data$Pct.Other)
plot(x = data$Mask.Never, y = data$Mask.Always)
ggplot(data) +
  geom_point(mapping = aes(x = Mask.Never, y = Mask.Always), col = 'firebrick3') +
  xlab("Mask Never") +
  ylab("Mask Always") +
  labs(title = paste0(
    "Correlation: ",
    round(cor(data$Mask.Never, data$Mask.Always), 2))
  )
ggplot(data) +
  geom_point(mapping = aes(x = Mask.Always, y = Mask.Positivity), col = 'slateblue4') +
  xlab("Mask Always") +
  ylab("Mask Mask Positivity") +
  labs(title = paste0(
    "Correlation: ",
    round(cor(data$Mask.Always, data$Mask.Positivity), 2))
  )
ggplot(data) +
  geom_point(mapping = aes(x = Mask.Always, y = Cases.Per.1k), col = 'slateblue4') +
  xlab("Mask Always") +
  ylab("Cases Per 1K")

ggplot(data) +
  geom_point(mapping = aes(x = Pct.Other, y = Cases.Per.1k), col = 'purple4') +
  xlab("Percentage Other") +
  ylab("Cases Per 1K")

ggplot(data) +
  geom_point(mapping = aes(x = Votes.Per.Capita, y = Cases.Per.1k), col = 'aquamarine4') +
  xlab("Votes Per Capita") +
  ylab("Cases Per 1K")

reg_func <- function(num_obs, num_preds, var_mat, beta0, g, v0, nu_0, sigma2_0, data_Xmat, data_y) {
  Identity_n <- diag(rep(1, num_obs))
  sigma_post_param1 <- (nu_0 + num_obs) / 2
  S2_g <- t(data_y) %*% (Identity_n - ((g / (g + 1)) * data_Xmat %*% var_mat %*% t(data_Xmat))) %*% data_y
  sigma_post_param2 <- 0.5 * ((nu_0 * sigma2_0) + S2_g)
  
  sigma2 <- 1 / rgamma(10000, sigma_post_param1, sigma_post_param2)
  
  beta_post_mean <- (g / (g + 1)) * ((beta0 / g) + (var_mat %*% t(data_Xmat) %*% data_y))
  beta_est <- matrix(rep(0, length(sigma2) * num_preds), nrow = length(sigma2))
  for (i in 1:length(sigma2)) {
    beta_post_var <- (g / (g + 1)) * sigma2[i] * var_mat
    beta_est[i, ] <- mvtnorm::rmvnorm(1, mean = beta_post_mean, sigma = beta_post_var)
  }
  
  beta_matrix <- matrix(rep(0, num_preds * 3), nrow = num_preds, ncol = 3)
  for (i in 1:num_preds) {
    beta_matrix[i, 1] <- quantile(beta_est[, i], probs = 0.025)
    beta_matrix[i, 2] <- beta_post_mean[i]
    beta_matrix[i, 3] <- quantile(beta_est[, i], probs = 0.975)
  }
  beta_matrix <- as.data.frame(
    beta_matrix,
    row.names = colnames(data_Xmat)
  )
  colnames(beta_matrix) <- c("2.5%", "Mean", "97.5%")
  return(beta_matrix)
}

## Reg v1: all numerical -------------------------------------------------------
data_regv1 <- data %>% 
  select(-X, -County, -State, -Date, -Cases, -Deaths, -Total.Votes, -Pct.Democrat)
## Reg v1: scale and center ----------------------------------------------------
data_regv1 <- data_regv1 %>% 
  mutate(
    Population = (Population - mean(data_regv1$Population)) / 
      sd(data_regv1$Population)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Never = (Mask.Never - mean(data_regv1$Mask.Never)) / 
      sd(data_regv1$Mask.Never)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Rarely = (Mask.Rarely - mean(data_regv1$Mask.Rarely)) / 
      sd(data_regv1$Mask.Rarely)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Sometimes = (Mask.Sometimes - mean(data_regv1$Mask.Sometimes)) /
      sd(data_regv1$Mask.Sometimes)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Frequently = (Mask.Frequently - mean(data_regv1$Mask.Frequently)) /
      sd(data_regv1$Mask.Frequently)
    )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Always = (Mask.Always - mean(data_regv1$Mask.Always)) /
      sd(data_regv1$Mask.Always)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Pct.Republican = (Pct.Republican - mean(data_regv1$Pct.Republican)) /
      sd(data_regv1$Pct.Republican)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Votes.Per.Capita = (Votes.Per.Capita - mean(data_regv1$Votes.Per.Capita)) /
      sd(data_regv1$Votes.Per.Capita)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Pct.Other = (Pct.Other - mean(data_regv1$Pct.Other)) /
      sd(data_regv1$Pct.Other)
  )
data_regv1 <- data_regv1 %>% 
  mutate(
    Mask.Positivity = (Mask.Positivity - mean(data_regv1$Mask.Positivity)) /
      sd(data_regv1$Mask.Positivity)
  )

## Reg v1: data analysis -------------------------------------------------------
write.csv(data_regv1, "regression_dfv1.csv")
data_regv1 <- read.csv("regression_dfv1.csv")
data_Xmat <- as.matrix(
  data_regv1 %>% 
    select(-X, -Fips, -Cases.Per.1k)
)
data_y <- as.matrix(data_regv1 %>% select(Cases.Per.1k))

# Initialize prior parameters
num_obs <- as.numeric(nrow(data_Xmat))
num_preds <- as.numeric(ncol(data_Xmat))
var_mat <- pinv(t(data_Xmat) %*% data_Xmat)
beta0 <- rep(0, num_preds)
g <- num_obs
v0 <- g * var_mat
nu_0 <- 2
sigma2_0 <- 1

# Perform regression analysis
beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0,
                     sigma2_0, data_Xmat, data_y)
beta_mat <- round(beta_mat, 3)
rownames(beta_mat) <- c("Population", "Mask Never", "Mask Rarely",
                        "Mask Sometimes", "Mask Frequently", "Mask Always",
                        "Pct Republican", "Pct Other", "Votes Per Capita",
                        "Mask Positivity", "Intercept")
knitr::kable(beta_mat)

grid.newpage()
vp1 <- viewport(x = 0, y = 0, width = 0.5, height = 1,
                just = c('left', 'bottom'))
vp2 <- viewport(x = 0.5, y = 0, width = 0.5, height = 1,
                just = c('left', 'bottom'))

a <- ggplot(data_regv1) +
  geom_point(mapping = aes(x = Mask.Never, y = Mask.Always), col = 'slateblue4') +
  xlab("Mask Always") +
  ylab("Percentage Other") +
  labs(title = paste0(
    "Correlation: ",
    round(cor(data_regv1$Mask.Never, data_regv1$Mask.Always), 2))
  )

b <- ggplot(data_regv1) +
  geom_point(mapping = aes(x = Mask.Always, y = Mask.Positivity),
             col = 'aquamarine4') +
  xlab("Mask Always") +
  ylab("Mask Positivity") +
  labs(title = paste0(
    "Correlation: ",
    round(cor(data_regv1$Mask.Always, data_regv1$Mask.Positivity), 2))
  )

print(a, vp = vp1)
print(b, vp = vp2)
## Variable selection ----------------------------------------------------------
#https://cran.r-project.org/web/packages/BAS/BAS.pdf
data_regvs <- data_regv1 %>% select(-Intercept)
num_obs <- as.numeric(nrow(data_regvs))
covid_bayes <- bas.lm(
  Cases.Per.1k ~ Population + Mask.Never + Mask.Rarely + Mask.Sometimes +
    Mask.Frequently + Mask.Always + Pct.Republican + Pct.Other +
    Votes.Per.Capita + Mask.Positivity,
  data = data_regvs,
  prior = "g-prior",
  alpha = num_obs,
  modelprior = uniform(),
  method = "deterministic"
)
beta <- coef(covid_bayes)
prob_matrix <- matrix(rep(0, 2 * length(beta$probne0)), ncol = 2)
prob_matrix[, 1] <- beta$namesx
prob_matrix[, 2] <- round(beta$probne0, 3)
prob_matrix <- as.data.frame(prob_matrix)
colnames(prob_matrix) <- c(
  "Variable",
  "Posterior Probability of Non-Zero Coefficient"
)
knitr::kable(prob_matrix)

num_sig_vars <- sum(prob_matrix[, 2] > 0.95)
mat_coeff_est <- matrix(rep(0, 4 * num_sig_vars),
                        nrow = num_sig_vars, ncol = 4)
confint_beta <- confint(beta)
index <- 1
for (i in 1:nrow(prob_matrix)) {
  if (prob_matrix[i, 2] > 0.95) {
    mat_coeff_est[index, 1] <- prob_matrix[i, 1]
    mat_coeff_est[index, 2] <- round(confint_beta[i, 1], 3)
    mat_coeff_est[index, 3] <- round(confint_beta[i, 3], 3)
    mat_coeff_est[index, 4] <- round(confint_beta[i, 2], 3)
    index <- index + 1
  }
}
colnames(mat_coeff_est) <- c("Variable", "2.5% Estimate",
                             "Mean Estimate", "97.5% Estimate")
knitr::kable(mat_coeff_est)

alpha_vec <- seq(from = 1, to = 2 * num_obs, by = 2)
inclusion_mat <- matrix(rep(0, 11 * (length(alpha_vec) + 1)),
                        nrow = 11, ncol = length(alpha_vec) + 1)
for (index in 1:length(alpha_vec)) {
  covid_bayes <- bas.lm(
    Cases.Per.1k ~ Population + Mask.Never + Mask.Rarely +
      Mask.Sometimes + Mask.Frequently + Mask.Always + Pct.Republican +
      Pct.Other + Votes.Per.Capita + Mask.Positivity,
    data = data_regvs,
    prior = "g-prior",
    alpha = alpha_vec[index],
    modelprior = uniform(),
    method = "deterministic"
  )
  beta <- coef(covid_bayes)
  if (index == 1) {
    inclusion_mat[, index] <- beta$namesx
  }
  inclusion_mat[, index + 1] <- round(beta$probne0, 3)
}
par(mfrow = c(1, 3))
plot(x = alpha_vec, y = inclusion_mat[7, 2:ncol(inclusion_mat)], type = 'l',
     xlab = "prior g-value", ylab = "Inclusion probability",
     main = "Mask Always")
plot(x = alpha_vec, y = inclusion_mat[9, 2:ncol(inclusion_mat)], type = 'l',
     xlab = "prior g-value", ylab = "Inclusion probability",
     main = " % Other")
plot(x = alpha_vec, y = inclusion_mat[10, 2:ncol(inclusion_mat)], type = 'l',
     xlab = "prior g-value", ylab = "Inclusion probability",
     main = "Votes Per Capita")


## Fitting Smaller Model + Diagnostics + Cross validation ------------------------------------------------------------
data_Xmat <- as.matrix(
  data_regv1 %>% 
    select(Mask.Always, Pct.Other, Votes.Per.Capita, Intercept)
)
data_y <- as.matrix(data_regv1 %>% select(Cases.Per.1k))
num_obs <- as.numeric(nrow(data_Xmat))
num_preds <- as.numeric(ncol(data_Xmat))
var_mat <- solve(t(data_Xmat) %*% data_Xmat)
beta0 <- rep(0, num_preds)
g <- num_obs
v0 <- g * var_mat
nu_0 <- 2
sigma2_0 <- 1

beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0, sigma2_0,
                     data_Xmat, data_y)
rownames(beta_mat) <- c("Mask Always", "Pct Other", "Votes Per Capita", "Intercept")
knitr::kable(round(beta_mat, 3))

## Autocorrelation between regressors
grid.newpage()
vp1 <- viewport(x = 0, y = 0, width = 0.333, height = 1,
                just = c('left', 'bottom'))
vp2 <- viewport(x = 0.333, y = 0, width = 0.333, height = 1,
                just = c('left', 'bottom'))
vp3 <- viewport(x = 0.666, y = 0, width = 0.333, height = 1,
                just = c('left', 'bottom'))

a <- ggplot(data_regv1) +
  geom_point(mapping = aes(x = Mask.Always, y = Pct.Other), col = 'slateblue4') +
  xlab("Mask Always") +
  ylab("Percentage Other") +
  labs(title = paste0(
    "Cor: ",
    round(cor(data_regv1$Mask.Always, data_regv1$Pct.Other), 2))
  )

b <- ggplot(data_regv1) +
  geom_point(mapping = aes(x = Mask.Always, y = Votes.Per.Capita),
             col = 'aquamarine4') +
  xlab("Mask Always") +
  ylab("Votes Per Capita") +
  labs(title = paste0(
    "Cor: ",
    round(cor(data_regv1$Mask.Always, data_regv1$Votes.Per.Capita), 2))
  )

c <- ggplot(data_regv1) +
  geom_point(mapping = aes(x = Votes.Per.Capita, y = Pct.Other),
             col = 'purple4') +
  xlab("Votes Per Capita") +
  ylab("Percentage Other") +
  labs(title = paste0(
    "Cor: ",
    round(cor(data_regv1$Votes.Per.Capita, data_regv1$Pct.Other), 2))
  )
print(a, vp = vp1)
print(b, vp = vp2)
print(c, vp = vp3)

## Residuals
beta_means <- beta_mat[, 2]
pred <- data_Xmat %*% beta_means
residuals <- data_y - pred
MS_res <- mean(residuals^2)
residuals <- residuals / sqrt(MS_res)
residuals_df <- matrix(
  c(seq(from = 1, to = length(residuals), by = 1), residuals),
  ncol = 2,
  byrow = FALSE
)
residuals_df <- as.data.frame(residuals_df)
colnames(residuals_df) <- c("observation", "residual")

ggplot(data = residuals_df, mapping = aes(x = residual)) +
  geom_histogram(bins = 30, color="black", fill="grey")
ggplot(data = residuals_df) +
  geom_point(mapping = aes(x = observation, y = residual), col = "#F05B5B") +
  xlab("Observation") +
  ylab("Residual Value")

par(mfrow=c(1, 2))
hist(residuals, main = "Residuals Histogram")
qqnorm(residuals, main = "Normal QQ Plot of Residuals")
abline(a = 0, b = 1, col = "red")

## Sensitivity analysis
reg_func <- function(num_obs, num_preds, var_mat, beta0, g, v0, nu_0,
                     sigma2_0, data_Xmat, data_y) {
  # compute posterior parameters
  sigma_post_param1 <- (nu_0 + num_obs) / 2
  identity_mat <- diag(rep(1, num_obs))
  intermediate_mat <- identity_mat - 
    ((g / (g + 1)) * data_Xmat %*% var_mat %*% t(data_Xmat))
  S2_g <- t(data_y) %*% intermediate_mat %*% data_y
  sigma_post_param2 <- 0.5 * ((nu_0 * sigma2_0) + S2_g)
  
  # compute posterior mean
  beta_post_mean <- (g / (g + 1)) * 
    ((beta0 / g) + (var_mat %*% t(data_Xmat) %*% data_y))
  return(beta_post_mean)
}
g_vec <- c(seq(from = 0.1, to = 75, by = 1))
num_obs <- as.numeric(nrow(data_Xmat))
num_preds <- as.numeric(ncol(data_Xmat))
var_mat <- solve(t(data_Xmat) %*% data_Xmat)
beta0 <- rep(0, num_preds)
beta_sens_mat <- matrix(rep(0, 4 * length(g_vec)), nrow = 4, ncol = length(g_vec))
for (i in 1:length(g_vec)) {
  print(i)
  g <- g_vec[i]
  v0 <- g * var_mat
  nu_0 <- 2
  sigma2_0 <- 1
  beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0, sigma2_0,
                       data_Xmat, data_y)
  beta_sens_mat[, i] <- beta_mat
}
par(mfrow = c(2, 2))
plot(x = g_vec, y = beta_sens_mat[1,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Mask Always")
plot(x = g_vec, y = beta_sens_mat[2,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Pct Other")
plot(x = g_vec, y = beta_sens_mat[3,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Votes Per Capita")
plot(x = g_vec, y = beta_sens_mat[4,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Intercept")

## Cross validation
bigk <- 10
test_size <- ceiling(2188 / bigk)
indices <- seq(from = 1, to = 2188, by = 1)
remaining_indices <- indices
mse_per_fold <- rep(0, bigk)
beta_avg <- rep(0, 4)
for (k in 1:bigk) {
  print(k)
  if (k == bigk) {
    test_indices <- remaining_indices
  } else {
    test_indices <- sample(remaining_indices, test_size, replace = FALSE)
  }
  
  remaining_indices <- setdiff(remaining_indices, test_indices)
  train_indices <- setdiff(indices, test_indices)
  
  data_Xmat <- as.matrix(
    data_regv1 %>% 
      select(Mask.Always, Pct.Other, Votes.Per.Capita, Intercept)
  )
  data_Xmat_train <- data_Xmat[train_indices, ]
  data_Xmat_test <- data_Xmat[test_indices, ]
  data_y <- as.matrix(data_regv1 %>% select(Cases.Per.1k))
  data_y_train <- data_y[train_indices, ]
  data_y_test <- data_y[test_indices, ]
  
  num_obs <- as.numeric(nrow(data_Xmat_train))
  num_preds <- as.numeric(ncol(data_Xmat_train))
  var_mat <- solve(t(data_Xmat_train) %*% data_Xmat_train)
  beta0 <- rep(0, num_preds)
  g <- num_obs
  v0 <- g * var_mat
  nu_0 <- 2
  sigma2_0 <- 1
  
  beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0, sigma2_0, data_Xmat_train, data_y_train)
  beta_means <- beta_mat[, 2]
  beta_avg <- beta_avg + beta_means
  #print("beta means")
  #print(beta_means)
  
  pred <- data_Xmat_test %*% beta_means
  #print("predicted values")
  #hist(pred)
  
  residuals <- data_y_test - pred
  mse_per_fold[k] <- mean(residuals^2)
  hist(residuals, xlab = "Residuals", main = "Histogram of Residuals")
  plot(x = seq(from = 1, to = length(residuals), by = 1), y = residuals, xlab = "Observation", ylab = "Residuals", main = "Residuals Plot")
  #print("error")
  #print(mse_per_fold[k])
  #print("y values")
  #hist(data_y_test)
}
mean(mse_per_fold)
print(beta_avg / 10)

## Vaccine model ---------------------------------------------------------------
new_cases <- read.csv("us-counties-2021.csv")
nc_filt <- filter(new_cases, date == "2021-06-12")
cvax <- read.csv("countyvax.csv")
cv_filt <- filter(cvax, Date == "12/07/2021")
cv_non <- filter(cv_filt, Recip_County %in% nc_filt$county)
cv_yes <- filter(cv_filt, !(Recip_County %in% nc_filt$county))
cv_yes %>% mutate(Recip_County = word(Recip_County , 1  , -2)) -> cv_yes
nrow(cv_non)
nrow(cv_yes)
nrow(cv_filt)
cv_yes_filt <- filter(cv_yes, Recip_County %in% nc_filt$county)
combined_data <- rbind(cv_non, cv_yes_filt)
nc_filt$state_abb <- setNames(
  state.abb, state.name)[as.vector(nc_filt$state)]
joined_data <- inner_join(combined_data, nc_filt, by=c("Recip_County"="county", "Recip_State"="state_abb"))
joined_data <- joined_data %>% mutate(FIPS = as.numeric(FIPS))
joined_data <- joined_data %>% select(FIPS, Series_Complete_Pop_Pct, cases_avg_per_100k)
reg_data <- read.csv("regression_dfv1.csv")
reg_vac_data <- merge(
  x = reg_data,
  y = joined_data,
  by.x = "Fips",
  by.y = "FIPS",
  all.x = TRUE,
  all.y = FALSE
) 
reg_vac_data <- reg_vac_data %>% 
  select(Fips, Mask.Always, Pct.Other, Votes.Per.Capita, Intercept, Series_Complete_Pop_Pct, cases_avg_per_100k)
reg_vac_data <- reg_vac_data %>% mutate(new_cases_avg_per_100k = cases_avg_per_100k)
reg_vac_data <- reg_vac_data %>% select(-cases_avg_per_100k)
reg_vac_data <- reg_vac_data %>% filter(!is.na(Series_Complete_Pop_Pct))
reg_vac_data <- reg_vac_data %>% filter(!is.na(new_cases_avg_per_100k))
reg_vac_data <- reg_vac_data %>% 
  mutate(Series_Complete_Pop_Pct = (Series_Complete_Pop_Pct - mean(reg_vac_data$Series_Complete_Pop_Pct))/sd(reg_vac_data$Series_Complete_Pop_Pct))
write.csv(reg_vac_data, "regression_df_vac.csv")

reg_vac_data <- read.csv("regression_df_vac.csv")
data_Xmat <- as.matrix(
  reg_vac_data %>% 
    select(-X, -Fips, -new_cases_avg_per_100k)
)
data_y <- as.matrix(reg_vac_data %>% select(new_cases_avg_per_100k))

num_obs <- as.numeric(nrow(data_Xmat))
num_preds <- as.numeric(ncol(data_Xmat))
var_mat <- solve(t(data_Xmat) %*% data_Xmat)
beta0 <- rep(0, num_preds)
g <- num_obs
v0 <- g * var_mat
nu_0 <- 2
sigma2_0 <- 1
reg_func <- function(num_obs, num_preds, var_mat, beta0, g, v0, nu_0,
                     sigma2_0, data_Xmat, data_y) {
  # compute posterior parameters
  sigma_post_param1 <- (nu_0 + num_obs) / 2
  identity_mat <- diag(rep(1, num_obs))
  intermediate_mat <- identity_mat - 
    ((g / (g + 1)) * data_Xmat %*% var_mat %*% t(data_Xmat))
  S2_g <- t(data_y) %*% intermediate_mat %*% data_y
  sigma_post_param2 <- 0.5 * ((nu_0 * sigma2_0) + S2_g)
  
  # sample variance values from posterior distribution
  sigma2 <- 1 / rgamma(10000, sigma_post_param1, sigma_post_param2)
  
  # compute posterior mean
  beta_post_mean <- (g / (g + 1)) * 
    ((beta0 / g) + (var_mat %*% t(data_Xmat) %*% data_y))
  
  # compute 95% CI for the coefficients
  beta_est <- matrix(rep(0, length(sigma2) * num_preds),
                     nrow = length(sigma2))
  for (i in 1:length(sigma2)) {
    beta_post_var <- (g / (g + 1)) * sigma2[i] * var_mat
    beta_est[i, ] <- mvtnorm::rmvnorm(1, mean = beta_post_mean,
                                      sigma = beta_post_var)
  }
  beta_matrix <- matrix(rep(0, num_preds * 3),
                        nrow = num_preds, ncol = 3)
  for (i in 1:num_preds) {
    beta_matrix[i, 1] <- quantile(beta_est[, i], probs = 0.025)
    beta_matrix[i, 2] <- beta_post_mean[i]
    beta_matrix[i, 3] <- quantile(beta_est[, i], probs = 0.975)
  }
  beta_matrix <- as.data.frame(beta_matrix,
                               row.names = colnames(data_Xmat))
  colnames(beta_matrix) <- c("2.5%", "Mean", "97.5%")
  return(beta_matrix)
}
beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0,
                     sigma2_0, data_Xmat, data_y)
knitr::kable(beta_mat)
plot(x = reg_vac_data$Series_Complete_Pop_Pct, y = reg_vac_data$new_cases_avg_per_100k)

reg_func <- function(num_obs, num_preds, var_mat, beta0, g, v0, nu_0,
                     sigma2_0, data_Xmat, data_y) {
  # compute posterior parameters
  sigma_post_param1 <- (nu_0 + num_obs) / 2
  identity_mat <- diag(rep(1, num_obs))
  intermediate_mat <- identity_mat - 
    ((g / (g + 1)) * data_Xmat %*% var_mat %*% t(data_Xmat))
  S2_g <- t(data_y) %*% intermediate_mat %*% data_y
  sigma_post_param2 <- 0.5 * ((nu_0 * sigma2_0) + S2_g)
  
  # compute posterior mean
  beta_post_mean <- (g / (g + 1)) * 
    ((beta0 / g) + (var_mat %*% t(data_Xmat) %*% data_y))
  return(beta_post_mean)
}
num_obs <- as.numeric(nrow(data_Xmat))
num_preds <- as.numeric(ncol(data_Xmat))
g_vec <- c(seq(from = 0.1, to = 100, by = 1))
var_mat <- solve(t(data_Xmat) %*% data_Xmat)
beta0 <- rep(0, num_preds)
beta_sens_mat <- matrix(rep(0, num_preds * length(g_vec)),
                        nrow = num_preds, ncol = length(g_vec))
for (i in 1:length(g_vec)) {
  g <- g_vec[i]
  v0 <- g * var_mat
  nu_0 <- 2
  sigma2_0 <- 1
  beta_mat <- reg_func(num_obs, num_preds, var_mat, beta0, g, v0, nu_0, sigma2_0,
                       data_Xmat, data_y)
  beta_sens_mat[, i] <- beta_mat
}
par(mfrow = c(2, 2))
plot(x = g_vec, y = beta_sens_mat[1,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Mask Always")
plot(x = g_vec, y = beta_sens_mat[2,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Pct Other")
plot(x = g_vec, y = beta_sens_mat[3,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Votes Per Capita")
plot(x = g_vec, y = beta_sens_mat[5,], type = "l",
     xlab = "g prior value", ylab = "Mean Coeff. Est.",
     main = "Vaccination Rate")
## Regression Prior Plot -------------------------------------------------------
x_vals <- seq(from = 0.1, to = 25, by = 0.01)
y_vals <- dinvgamma(x_vals, shape = 1, rate = 1)
plot(x_vals, y_vals, type = 'l', xlab = "x", ylab = "Prior Density")