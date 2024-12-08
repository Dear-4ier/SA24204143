## ----include=FALSE------------------------------------------------------------
library(dplyr)
library(Rcpp)
library(microbenchmark)
library(corrplot)
library(ggplot2)
library(boot)
library(DAAG)
library(lpSolve)
library(bootstrap)
library(coda)
library(stats)
library(bench)

## -----------------------------------------------------------------------------
x = rnorm(100000)
hist(x)

## -----------------------------------------------------------------------------
estimate_pi <- function(n) {
  I <- 0
  for (i in 1:n) {
    # Generate random x and y coordinates between -1 and 1
    x <- runif(1, -1, 1)
    y <- runif(1, -1, 1)
    # Check if the point is inside the unit circle
    if (x^2 + y^2 <= 1) {
      I <- I + 1
    }
  }
  # Estimate pi
  return (4 * I / n)
}

# Generate different numbers of points
points <- c(10^1, 10^2, 10^3, 10^4, 10^5)
pi_estimates <- sapply(points, estimate_pi)

# Display results
data.frame(Number_of_Points = points, Estimated_Pi = pi_estimates)

## -----------------------------------------------------------------------------
(R = cor(state.x77))
corrplot(R,diag=F)

## -----------------------------------------------------------------------------
# generate random samples from a Rayleigh $(\sigma)$ distribution
set.seed(123456)
rayleigh <- function(sigma) {
  n <- 1000
  U <- runif(n)
  return(sigma * sqrt(-2 * log(U)))
}
# Generate Rayleigh $(\sigma)$ samples for several choices of $\sigma>0$
sigma_values <- c(2, 3, 4)

for (sigma in sigma_values) {
  samples <- rayleigh(sigma)
  
  # check the histogram
  hist(
    samples,
    main = paste("Histogram for Rayleigh(", sigma, ") Distribution", sep = ""),
    xlab = "Value",
    ylab = "Frequency"
  )
  
  # Theoretical mode is sigma
  cat("Theoretical mode for sigma =", sigma, "is", sigma, "\n")
  
  # Estimated mode from histogram
  estimated_mode <- density(samples)$x[which.max(density(samples)$y)]
  cat("mode of the generated samples for sigma =",
      sigma,
      "is",
      estimated_mode,
      "\n")
}

## -----------------------------------------------------------------------------
set.seed(123456)
normal_mixture <- function(p_1) {
  p_2 <- 1 - p_1
  n <- 1000
  sample_with_mean_0 <- rnorm(n * p_1, 0, 1)
  sample_with_mean_3 <- rnorm(n * p_2, 3, 1)
  mixture_sample = c(sample_with_mean_0, sample_with_mean_3)
  
  # Graph the histogram of the sample with density superimposed
  p = ggplot(data = data.frame(x = mixture_sample), aes(x = x, label = "Density")) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = 0.5,
      color = "black",
      fill = "gray"
    ) +
    stat_density(geom = "line", color = "red") +
    ggtitle(paste0("Histogram of Mixture Distribution with p1 = ", p_1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  return(p)
}

## -----------------------------------------------------------------------------
for (p_1 in c(0.75, 0.5, 0.25, 0.1)) {
  normal_mixture(p_1)
}

## -----------------------------------------------------------------------------
set.seed(123456)
compound_poisson_gamma <- function(lambda, shape, scale, t) {
  N <- rpois(1, lambda * t)
  if (N == 0) {
    return(0)
  } else {
    Y <- rgamma(N, shape = shape, scale = scale)
    return(sum(Y))
  }
}

## -----------------------------------------------------------------------------
lambda <- c(1, 2, 3)
shape <- c(2, 3, 4)
scale <- c(1, 2, 3)
num_sims <- 1000
t <- 10
results <- data.frame(
  lambda = rep(lambda, each = length(shape) * length(scale)),
  shape = rep(rep(shape, each = length(scale)), length(lambda)),
  scale = rep(scale, length(lambda) * length(shape)),
  simulated_mean = NA,
  theoretical_mean = NA,
  mean_error = NA,
  simulated_var = NA,
  theoretical_var = NA,
  var_error = NA
)
for (i in 1:nrow(results)) {
  sims <- replicate(
    num_sims,
    compound_poisson_gamma(results$lambda[i], results$shape[i], results$scale[i], t)
  )
  results$simulated_mean[i] <- mean(sims)
  results$simulated_var[i] <- var(sims)
  results$theoretical_mean[i] <- results$lambda[i] * t * (results$shape[i] * results$scale[i])
  results$theoretical_var[i] <- results$lambda[i] * t * (results$shape[i] * (results$shape[i] +
                                                                               1) * results$scale[i] ^ 2)
  results$mean_error <- abs(results$simulated_mean[i] - results$theoretical_mean[i]) /
    results$theoretical_mean[i]
  results$var_error <- abs(results$simulated_var[i] - results$theoretical_var[i]) /
    results$theoretical_var[i]
}

## -----------------------------------------------------------------------------
print(results, width = 1000)

## -----------------------------------------------------------------------------
ggplot(results, aes(x = lambda, y = simulated_mean, color = "Simulation")) +
  geom_point() +
  geom_line(aes(y = theoretical_mean, color = "Theory")) +
  facet_grid(shape ~ scale) +
  labs(title = "Comparison of Mean", x = "lambda", y = "Mean") +
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
ggplot(results, aes(x = lambda, y = simulated_var, color = "Simulation")) +
  geom_point() +
  geom_line(aes(y = theoretical_var, color = "Theory")) +
  facet_grid(shape ~ scale) +
  labs(title = "Comparison of Variance", x = "lambda", y = "Variance") +
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
set.seed(12345)
monte_carlo_beta_cdf_mean <- function(x) {
  u <- rbeta(1e4, shape1 = 3, shape2 = 3)
  return(mean(u <= x))
}

## -----------------------------------------------------------------------------
x_values <- seq(0.1, 0.9, by = 0.1)
monte_carlo_estimates <- numeric(length(x_values)) # Monte Carlo估计值
pbeta_values <- numeric(length(x_values)) # pbeta 函数值

for (i in seq_along(x_values)) {
  monte_carlo_estimates[i] <- format(round(monte_carlo_beta_cdf_mean(x_values[i]), 4), nsmall = 4)
  pbeta_values[i] <- format(round(pbeta(
    x_values[i], shape1 = 3, shape2 = 3
  ), 4), nsmall = 4)
}

# 输出结果进行比较
df <- data.frame(x = x_values,
                 monte_carlo = monte_carlo_estimates,
                 pbeta = pbeta_values)
df

## -----------------------------------------------------------------------------
ggplot(df, aes(x = x)) +
  geom_line(aes(y = as.numeric(monte_carlo), color = "蒙特卡洛估计")) +
  geom_line(aes(y = as.numeric(pbeta), color = "pbeta 值")) +
  geom_point(aes(y = as.numeric(monte_carlo), color = "蒙特卡洛估计")) +
  geom_point(aes(y = as.numeric(pbeta), color = "pbeta 值")) +
  labs(x = "x 值", y = "CDF 值", color = "曲线类型") +
  scale_color_manual(values = c("蒙特卡洛估计" = "blue", "pbeta 值" = "red"))

## -----------------------------------------------------------------------------
sigma <- 1
n <- 1e4

## -----------------------------------------------------------------------------
# Implement a function to generate samples from a Rayleigh($\sigma)$ distribution, using antithetic variables.
antithetic_rayleigh <- function(sigma, n) {
  u <- runif(n)
  X <- sigma * sqrt(-2 * log(u))
  X_prime <- sigma * sqrt(-2 * log(1 - u))
  return((X + X_prime) / 2)
}

# 对偶变量法采样
Sample_antithetic <- antithetic_rayleigh(sigma, n)
var_antithetic <- var(Sample_antithetic)

## -----------------------------------------------------------------------------
# 生成瑞利分布样本的函数(上次作业)
rayleigh <- function(sigma, n) {
  U <- runif(n)
  return(sigma * sqrt(-2 * log(U)))
}

# 独立采样
X_1 <- rayleigh(sigma, n)
X_2 <- rayleigh(sigma, n)
Sample_independent <- (X_1 + X_2) / 2
var_independent <- var(Sample_independent)

## -----------------------------------------------------------------------------
reduced_percentage <- (var_independent - var_antithetic) / var_independent * 100
cat("对于参数 sigma =", sigma, "和样本数量 n =", n, "\n")
cat("独立采样的样本方差为", var_independent, "\n")
cat("对偶变量法采样的样本方差为", var_antithetic, "\n")
cat("方差减少的百分比为：", reduced_percentage, "%\n")

## -----------------------------------------------------------------------------
g <- function(x) {
  ifelse(x > 1, x ^ 2 / sqrt(2 * pi) * exp(-x ^ 2 / 2), 0)
}

## -----------------------------------------------------------------------------
f1 <- function(x) {
  ifelse(x > 1, exp(-(x - 1)), 0)
}
f2 <- function(x) {
  ifelse(x > 1, 2 * exp(-2 * (x - 1)), 0)
}

## -----------------------------------------------------------------------------
# 生成数据
set.seed(12345)
x <- seq(1.01, 5, by = 0.1)
data <- data.frame(
  x = x,
  g = g(x),
  f1 = f1(x),
  f2 = f2(x)
)

# 绘制图像
ggplot(data, aes(x = x)) +
  geom_line(aes(y = g, color = "g(x)")) +
  geom_line(aes(y = f1, color = "f1")) +
  geom_line(aes(y = f2, color = "f2")) +
  labs(x = "x", y = "函数值", color = "函数") +
  scale_color_manual(values = c(
    "g(x)" = "red",
    "f1" = "blue",
    "f2" = "black"
  ))

## -----------------------------------------------------------------------------
m <- 1e4
# Importance sampling for f1
samples1 <- rexp(m, 1) + 1
gf1 <- g(samples1) / f1(samples1)
var1 <- var(gf1)

# Importance sampling for f2
samples2 <- rexp(m, 2) + 1
gf2 <- g(samples2) / f1(samples2)
var2 <- var(gf2)

# 打印结果
cat("f1在估计积分时产生的方差为：", var1, "\n")
cat("f2在估计积分时产生的方差为：", var2, "\n")

# Determine which variance is smaller
if (var1 < var2) {
  cat("f1在估计积分时产生的方差更小")
} else {
  cat("f2在估计积分时产生的方差更小")
}

## -----------------------------------------------------------------------------
qsort <- function(x) {
  if (length(x) <= 1)
    return(x)
  pivot <- x[sample(length(x), 1)]
  left <- x[x < pivot]
  right <- x[x > pivot]
  return(c(qsort(left), pivot, qsort(right)))
}

## ----eval=FALSE---------------------------------------------------------------
# set.seed(12345)
# n_values <- c(1 * 10 ^ 4, 2 * 10 ^ 4, 4 * 10 ^ 4, 6 * 10 ^ 4, 8 * 10 ^ 4)
# m <- 100 #  100 simulations
# times <- numeric(length(n_values))
# 
# for (i in seq_along(n_values)) {
#   n <- n_values[i]
#   total_time <- 0
#   for (j in 1:m) {
#     x <- 1:n
#     start_time <- Sys.time()
#     qsort(x)
#     end_time <- Sys.time()
#     total_time <- total_time + as.numeric(difftime(end_time, start_time, units = "secs"))
#   }
#   times[i] <- total_time / m
# }

## ----eval=FALSE---------------------------------------------------------------
# t_values <- n_values * log(n_values)
# time_fit <- lm(times ~ t_values)
# summary(time_fit)

## ----eval=FALSE---------------------------------------------------------------
# ggplot(data = data.frame(n = n_values, times = times, t = t_values), aes(x = t, y = times)) + geom_point() +
#   geom_smooth(method = "lm", se = FALSE) + labs(x = "nlog(n)", y = "平均计算时间")

## -----------------------------------------------------------------------------
set.seed(12345)
skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  std_dev <- sd(x)
  m3 <- sum((x - mean_x) ^ 3) / n
  return(m3 / std_dev ^ 3)
}

## -----------------------------------------------------------------------------
n <- 1000
m <- 1e4

skewness_values <- numeric(m)
for (i in 1:m) {
  x <- rnorm(n)
  skewness_values[i] <- skewness(x)
}

# 计算分位数和标准误差
quantiles <- c(0.025, 0.05, 0.95, 0.975)
estimated_quantiles <- quantile(skewness_values, quantiles)

density <- dnorm(estimated_quantiles)
var <- quantiles * (1 - quantiles) / (n * density ^ 2)
se <- sqrt(var)

## -----------------------------------------------------------------------------
large_sample_quantiles <- qnorm(quantiles, mean = 0, sd = sqrt(6 / n))

## -----------------------------------------------------------------------------
results <- data.frame(
  分位数 = quantiles,
  标准差 = se,
  估计分位数 = estimated_quantiles,
  大样本分位数 = large_sample_quantiles
)

print(results)

## -----------------------------------------------------------------------------
bivariate_normal_test <- function(n, m) {
  p_values_Pearson <- rep(NA, m)
  p_values_spearman <- rep(NA, m)
  p_values_kendall <- rep(NA, m)
  for (i in 1:m) {
    x <- rnorm(n)
    y <- rnorm(n)
    cor_res <- cor.test(x, y)
    p_values_Pearson[i] <- cor_res$p.value
    spearman_res <- cor.test(x, y, method = 'spearman')
    p_values_spearman[i] <- spearman_res$p.value
    kendall_res <- cor.test(x, y, method = 'kendall')
    p_values_kendall[i] <- kendall_res$p.value
  }
  list(
    p_values_Pearson = p_values_Pearson,
    p_values_spearman = p_values_spearman,
    p_values_kendall = p_values_kendall
  )
}

## -----------------------------------------------------------------------------
compute_power <- function(p_values) {
  mean(p_values < 0.05)
}

## -----------------------------------------------------------------------------
print_power_bivariate_normal <- function(n, m) {
  test_results <- bivariate_normal_test(n, m)
  power_Pearson <- compute_power(test_results$p_values_Pearson)
  power_spearman <- compute_power(test_results$p_values_spearman)
  power_kendall <- compute_power(test_results$p_values_kendall)
  cat('在二元正态分布下，样本量为', n, '，重复次数为', m, '时：\n')
  cat('Pearson相关检验的功效：', power_Pearson, '\n')
  cat('Spearman检验的功效：', power_spearman, '\n')
  cat('Kendall检验的功效：', power_kendall, '\n')
}

## -----------------------------------------------------------------------------
# 生成基于指数分布的数据并进行检验的函数
exponential_test <- function(n, m) {
  p_values_Pearson <- rep(NA, m)
  p_values_spearman <- rep(NA, m)
  p_values_kendall <- rep(NA, m)
  for (i in 1:m) {
    x <- rexp(n)
    y <- rexp(n)
    cor_res <- cor.test(x, y)
    p_values_Pearson[i] <- cor_res$p.value
    spearman_res <- cor.test(x, y, method = 'spearman')
    p_values_spearman[i] <- spearman_res$p.value
    kendall_res <- cor.test(x, y, method = 'kendall')
    p_values_kendall[i] <- kendall_res$p.value
  }
  list(
    p_values_Pearson = p_values_Pearson,
    p_values_spearman = p_values_spearman,
    p_values_kendall = p_values_kendall
  )
}


## -----------------------------------------------------------------------------
print_power_exponential <- function(n, m) {
  test_results <- exponential_test(n, m)
  power_Pearson <- compute_power(test_results$p_values_Pearson)
  power_spearman <- compute_power(test_results$p_values_spearman)
  power_kendall <- compute_power(test_results$p_values_kendall)
  cat('在非正态分布下，样本量为', n, '，重复次数为', m, '时：\n')
  cat('Pearson相关检验的功效：', power_Pearson, '\n')
  cat('Spearman检验的功效：', power_spearman, '\n')
  cat('Kendall检验的功效：', power_kendall, '\n')
}

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 100
m <- 1000
print_power_bivariate_normal(n, m)# 二元正态分布
cat('\n')
print_power_exponential(n, m)# 非正态分布

## -----------------------------------------------------------------------------
# 计算 Z 统计量
p1 <- 0.651
p2 <- 0.676
n <- 10000
se1 <- sqrt(p1 * (1 - p1) / n)
se2 <- sqrt(p2 * (1 - p2) / n)
se <- sqrt(se1 ^ 2 + se2 ^ 2)
z <- (p2 - p1) / se
# 计算 P 值
p_value <- 2 * pnorm(-abs(z))

# 判断是否拒绝零假设
cat('p=', p_value,'\n')
if (p_value < 0.05) {
  cat("在 0.05 水平上，两种方法的功效有显著差异。\n")
} else {
  cat("在 0.05 水平上，两种方法的功效没有显著差异。\n")
}

## ----eval = FALSE-------------------------------------------------------------
# set.seed(123)
# N <- 1000
# n0 <- 950
# n1 <- 50
# m <- 10000
# alpha <- 0.1
# 
# # 循环生成模拟的 p 值
# p_values <- matrix(NA, nrow = m, ncol = N)
# for (i in 1:m) {
#   for (j in 1:N) {
#     if (j <= n0) {
#       p_values[i, j] <- runif(1)
#     } else {
#       p_values[i, j] <- rbeta(1, 0.1, 1)
#     }
#   }
# }
# 
# # 调整后的 p 值
# p_bonferroni <- p_values * N
# p_bh <- matrix(NA, nrow = m, ncol = N)
# for (i in 1:m) {
#   sorted_p <- sort(p_values[i, ])
#   for (j in 1:N) {
#     p_bh[i, j] <- sorted_p[j] * N / (j)
#   }
# }
# 
# # 计算
# FWER_bonferroni <- mean(rowSums(p_bonferroni <= alpha) > 0)
# FDR_bonferroni <- mean(apply(p_bonferroni <= alpha, 1, function(x)
#   mean(p_values[x, ] <= alpha / N)))
# TPR_bonferroni <- mean(apply(p_bonferroni <= alpha, 1, function(x)
#   mean(p_values[x, (n0 + 1):N] <= alpha / N)))
# 
# FWER_bh <- mean(rowSums(p_bh <= alpha) > 0)
# FDR_bh <- mean(apply(p_bh <= alpha, 1, function(x)
#   mean(p_values[x, ] <= alpha / N)))
# TPR_bh <- mean(apply(p_bh <= alpha, 1, function(x)
#   mean(p_values[x, (n0 + 1):N] <= alpha / N)))
# 
# # 结果
# result <- matrix(
#   c(FWER_bonferroni, FWER_bh, FDR_bonferroni, FDR_bh, TPR_bonferroni, TPR_bh),
#   nrow = 3,
#   ncol = 2,
#   byrow = FALSE
# )
# colnames(result) <- c("Bonferroni correction", "B-H correction")
# rownames(result) <- c("FWER", "FDR", "TPR")
# result

## -----------------------------------------------------------------------------
times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
# 危险率lambda的极大似然估计值
n <- length(times)
lambda_hat <- n/sum(times)
cat("MLE of the hazard rate lambda：", lambda_hat, "\n")

## -----------------------------------------------------------------------------
set.seed(123)
# Bootstrap
B <- 1000 # 设置自助抽样次数
lambda_boot <- numeric(B)
for (i in 1:B) {
  boot_sample <- sample(times, replace = TRUE)
  lambda_boot[i] <- length(boot_sample)/sum(boot_sample)
}

bias <- mean(lambda_boot) - lambda_hat
se <- sd(lambda_boot)

cat("Bias:", bias, "\n")
cat("Standard error:", se)

## -----------------------------------------------------------------------------
set.seed(123)
# 计算平均故障间隔时间
mean_time <- 1/lambda_hat

# 计算 95%置信区间
## 标准正态方法
z <- qnorm(0.975)
normal_interval <- mean_time + c(-z * se, z * se) * (1/lambda_hat^2)

## 基本方法
basic_interval <- quantile(lambda_boot, c(0.025, 0.975))
basic_interval <- 1/basic_interval

## 百分位数方法
percentile_interval <- quantile(1/lambda_boot, c(0.025, 0.975))

## BCa方法
boot_obj <- boot(data = times, statistic = function(x, i) {
  n <- length(x[i])
  return(n/sum(x[i]))
}, R = B)
bca_interval <- boot.ci(boot_obj, type = "bca")$bca[4:5]
bca_interval <- 1/bca_interval

cat("标准正态法 95%置信区间：[", normal_interval[1], ",", normal_interval[2], "]\n")
cat("基本方法 95%置信区间：[", basic_interval[2], ",", basic_interval[1], "]\n")
cat("分位数法 95%置信区间：[", percentile_interval[1], ",", percentile_interval[2], "]\n")
cat("BCa法 95%置信区间：[", bca_interval[2], ",", bca_interval[1], "]\n")

## -----------------------------------------------------------------------------
sc<-scor
set.seed(123)

# 计算theta的函数
estimate_theta <- function(x) {
    cov(x) %>%
    prcomp() %>%
    {.$sdev[1]/sum(.$sdev)}
}

n<-nrow(sc)
theta.j<- numeric(n)
for (i in 1:n){    
  theta.j[i]<-estimate_theta(sc[-i,])
}
theta.hat<-estimate_theta(sc)

bias<-(n-1)*(mean(theta.j)-theta.hat) #BIAS
se<-sqrt((n-1)*var(theta.j)) #SE

cat("Bias:", round(bias,3), "\n")
cat("standard error:", round(se,3),"\n")

## -----------------------------------------------------------------------------
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1

  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2

  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3

  J4 <- lm(y ~ x + I(x^2)+I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +J4$coef[3] * chemical[k]^2 +J4$coef[3] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}

round(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)),3)

## -----------------------------------------------------------------------------
# 函数：拟合模型，返回模型列表
fit_models <- function(y, x) {
  L1 <- lm(y ~ x)
  L2 <- lm(y ~ x + I(x ^ 2))
  L3 <- lm(log(y) ~ x)
  L4 <- lm(y ~ x + I(x ^ 2) + I(x ^ 3))
  return(list(L1, L2, L3, L4))
}

# maximum adjusted R^2 选择模型
fit_models(magnetic, chemical) %>%
  sapply(., function(model)
    round(summary(model)$adj.r.squared, 3)) %>%
  cat("Adjusted-R^2:\n", .)

## -----------------------------------------------------------------------------
# 计算 Cramer-von Mises 统计量
cvm_stat <- function(x, y) {
  n <- length(x)
  m <- length(y)
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  term1 <- sum((Fn(x) - Gm(x)) ^ 2)
  term2 <- sum((Fn(y) - Gm(y)) ^ 2)
  return((m * n) / (m + n) ^ 2 * (term1 + term2))
}

# 置换检验函数
permutation_test <- function(x, y, B = 1000) {
  observed_stat <- cvm_stat(x, y)
  z <- c(x, y)
  n <- length(x)
  m <- length(y)
  stats <- numeric(B)
  for (i in 1:B) {
    permuted <- sample(z)
    x_perm <- permuted[1:n]
    y_perm <- permuted[(n + 1):(n + m)]
    stats[i] <- cvm_stat(x_perm, y_perm)
  }
  p_value <- mean(stats >= observed_stat)
  return(list(observed_statistic = observed_stat, p_value = p_value))
}

# 结果打印函数
print_permutation_test_result <- function(x, y) {
  set.seed(123)
  permutation_test(x, y)$observed_statistic %>%
    cat("Observed Cramer-von Mises statistic:", ., "\n")
  p_value <- permutation_test(x, y)$p_value
  cat("P-value:", p_value, "\n")
  if (p_value < 0.05)
    cat("p < 0.05，拒绝原假设。\n")
  else
    cat("p > 0.05，不能拒绝原假设\n")
}

## -----------------------------------------------------------------------------
# 数据
attach(chickwts)
eg_x <- sort(as.vector(weight[feed == "soybean"]))
eg_y <- sort(as.vector(weight[feed == "linseed"]))

print_permutation_test_result(eg_x, eg_y)

## -----------------------------------------------------------------------------
# 置换检验函数
permutation_test_spearman <- function(x, y, B = 1000) {
  observed_stat <- cor(x, y, method = "spearman")
  n <- length(x)
  stats <- numeric(B)
  for (i in 1:B) {
    permuted_y <- sample(y)
    stats[i] <- cor(x, permuted_y, method = "spearman")
  }
  p_value <- mean(abs(stats) >= abs(observed_stat))
  return(list(observed_statistic = observed_stat, p_value = p_value))
}

# 生成随机数据
set.seed(123)
x <- rnorm(1e3)
y <- rnorm(1e3)

perm_result <- permutation_test_spearman(x, y)
cor_test_result <- cor.test(x, y, method = "spearman")

# 结果打印
round(perm_result$observed_statistic, 3) %>%
  cat("Observed Spearman correlation statistic from permutation test:",
      .,
      "\n")
round(perm_result$p_value, 3) %>%
  cat("P-value from permutation test:", ., "\n")
round(cor_test_result$p.value, 3) %>%
  cat("P-value from cor.test:", ., "\n")

## -----------------------------------------------------------------------------
set.seed(123)

# 目标密度函数（标准柯西分布）
target_density <- function(x) {
  return(1 / (pi * (1 + x^2)))
}

# Metropolis-Hastings
metropolis_hastings <- function(start_value, iterations) {
  samples <- numeric(iterations)
  current <- start_value
  for (i in 1:iterations) {
    proposal <- rnorm(1, current, 1)
    acceptance_ratio <- target_density(proposal) / target_density(current)
    if (runif(1) < acceptance_ratio) {
      current <- proposal
    }
    samples[i] <- current
  }
  return(mcmc(samples))
}

# Metropolis-Hastings并监测收敛
run_and_monitor <- function() {
  n_chains <- 2
  chain_length <- 10000
  chains <- list()
  for (i in 1:n_chains) {
    chains[[i]] <- metropolis_hastings(rnorm(1), chain_length)
  }
  combined_chains <- mcmc.list(chains)
  gelman_rubin <- gelman.diag(combined_chains, multivariate = FALSE)$psrf[1]
  while (gelman_rubin > 1.2) {
    new_chains <- list()
    for (i in 1:n_chains) {
      new_chains[[i]] <- c(chains[[i]], metropolis_hastings(chains[[i]][length(chains[[i]])], chain_length))
    }
    chains <- new_chains
    combined_chains <- mcmc.list(chains)
    gelman_rubin <- gelman.diag(combined_chains, multivariate = FALSE)$psrf[1]
  }
  return(chains)
}

# 结果输出
chains <- run_and_monitor()
all_samples <- unlist(chains)
trimmed_samples <- all_samples[-(1:1000)]
generated_deciles <- quantile(trimmed_samples, probs = seq(0.1, 0.9, by = 0.1))
true_deciles <- qcauchy(p = seq(0.1, 0.9, by = 0.1))
cat("生成的观测值十分位数：", generated_deciles, "\n")
cat("标准柯西分布十分位数：", true_deciles, "\n")

## -----------------------------------------------------------------------------
df <- function(x, y) {
  # 一般二项式系数
  gamma(n + 1) / (gamma(x + 1) * gamma(n - x + 1)) * y ^ (x + a - 1) * (1 - y) ^
    (n - x + b - 1)
}

# Gibbs 采样函数
gibbs_sampler <- function(niter = 10000) {
  d <- 2
  x1 <- matrix(0, nrow = niter, ncol = d)
  x2 <- matrix(0.1, nrow = niter, ncol = d)
  
  for (i in 2:niter) {
    x1[i, ] <- x1[i - 1, ] %>% {
     .[1] <- rbinom(1, n,.[2])
     .[2] <- rbeta(1,.[1] + a, n -.[1] + b)
     .
    }
    x2[i, ] <- x2[i - 1, ] %>% {
     .[1] <- rbinom(1, n,.[2])
     .[2] <- rbeta(1,.[1] + a, n -.[1] + b)
     .
    }
  }
  chain1 <- as.mcmc(x1)
  chain2 <- as.mcmc(x2)
  return(mcmc.list(chain1, chain2))
}

# 结果
set.seed(123)
n <- 100
a <- 2
b <- 3
converged <- FALSE
while (!converged) {
  combined_chain <- gibbs_sampler() %>% gelman.diag()
  if (combined_chain$mpsrf[1] < 1.2) {
    converged <- TRUE
  } else {
    combined_chain <- gibbs_sampler()
  }
}
cat("Gibbs 采样收敛情况：", converged, "\n")
cat("目标联合密度函数中的参数 n =", n, ", a =", a, ", b =", b, "\n")
cat("生成的马尔可夫链中前几个样本点（x,y）：\n")
head(combined_chain[[1]])

## -----------------------------------------------------------------------------
# function to compute the kth term
compute_kth_term <- function(k, a) {
  d <- length(a)
  norm_a <- sqrt(sum(a^2))  # Euclidean norm
  term <- (-1)^k / (factorial(k) * 2^k) * (norm_a^(2*k+2)) / ((2*k+1)*(2*k+2)) *
    gamma((d+1)/2) * gamma(k + 3/2) / gamma(k + d/2 + 1)
  return(term)
}

# function to compute sum
compute_sum <- function(a, tol = 1e-10) {
  norm_a <- sqrt(sum(a^2))  # Euclidean norm
  sum_series <- 0
  k <- 0
  term <- compute_kth_term(k, a)
  while (abs(term) > tol) {
    sum_series <- sum_series + term
    k <- k + 1
    term <- compute_kth_term(k, a)
  }
  return(sum_series)
}

# Evaluate the sum when a=(1,2)
a <- c(1, 2)
result <- compute_sum(a)
cat('sum when a = (1,2) :', result)

## -----------------------------------------------------------------------------
# 定义方程的函数
define_equation <- function(k, a){
  int_func <- function(u, n){(1 + u^2/(n - 1))^(-n/2)}
  get_c <- function(n, a){sqrt(a^2 * n / (n + 1 - a^2))}
  expr <- function(n, a) {
    this_int_func <- function(u){
      int_func(u, n)}
    c_val <- get_c(n - 1, a)
    2/sqrt(pi*(n - 1)) * exp(lgamma(n/2)-lgamma((n - 1)/2)) * 
      integrate(this_int_func, lower = 0, upper = c_val)$value}
  
  lhs <- expr(k, a)
  rhs <- expr(k + 1, a)
  return(lhs - rhs)
}

# 解方程的函数
solve_eq <- function(k) {
  f <- function(a){return(define_equation(k, a))}
  eps <- 1e-2
    if (f(eps) < 0 && f(sqrt(k) - eps) > 0 || f(eps) > 0 && f(sqrt(k) - eps) < 0) {
    return(uniroot(f, interval = c(eps, sqrt(k) - eps))$root)
  } else {
    return(NA)
  }
}

# 结果输出
r11.5 <- sapply(c(4:25, 100, 500, 1000), function(k) {
  solve_eq(k) %>% 
    return()
})
r11.5

## -----------------------------------------------------------------------------
# EX 11.4
findIntersection = function (k) {
  s_k.minus.one = function (a) {1 - pt(sqrt(a ^ 2 * (k - 1) / (k - a ^ 2)), df = k - 1)}
  s_k = function (a) {1 - pt(sqrt(a ^ 2 * k / (k + 1 - a ^ 2)), df = k)}
  f = function (a) {s_k(a) - s_k.minus.one(a)}
  
  eps = 1e-2
  return(uniroot(f, interval = c(eps, sqrt(k) - eps))$root)
}

r11.4 <- sapply(c(4:25, 100, 500, 1000), function (k) {
  findIntersection(k)
})

r11.4

## -----------------------------------------------------------------------------
# E步函数
E_step <- function(Y, tau, lambda0) {
  n <- length(Y)
  loglikelihood <- ifelse(Y < tau, log(lambda0) - lambda0 * Y, log(lambda0) - lambda0 * (tau + 1 / lambda0))
  return(sum(loglikelihood))
}

# M步函数
M_step <- function(Y, tau, lambda0) {
  n <- length(Y)
  expected_Ti <- ifelse(Y < tau, Y, tau + 1 / lambda0)
  return(n / sum(expected_Ti))
}

# EM算法
EM_algorithm <- function(Y, tau, lambda0) {
  while (TRUE) {
    Q_old <- E_step(Y, tau, lambda0)
    lambda1 <- M_step(Y, tau, lambda0)
    Q_new <- E_step(Y, tau, lambda1)
    # 同时考虑lambda的变化和对数似然期望的变化
    if (abs(lambda1 - lambda0) < 1e-6 && abs(Q_new - Q_old) < 1e-6) {
      break
    }
    lambda0 <- lambda1
  }
  return(lambda1)
}

# 观测数据的似然函数
loglikelihood <- function(Y, tau, lambda) {
  loglik <- ifelse(Y < tau, log((1 - exp(-lambda * tau)) * lambda * exp(-lambda * Y)), log(exp(-lambda * tau)))
  return(-sum(loglik))
}

# 结果输出
Y <- c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)
tau <- 1
lambda0 <- 1

EM_algorithm(Y, tau, lambda0) %>%
  round(., 3) %>%
  cat('EM算法估计的lambda值为：', ., '\n')

result_optim <- optim(
  par = lambda0, fn = loglikelihood, method = 'BFGS', Y = Y, tau = tau)
  lambda_mle <- round(result_optim$par, 3) %>%
    cat('观测数据MLE估计的lambda值为：',., '\n')

## -----------------------------------------------------------------------------
obj_coeffs <- c(4, 2, 9)
constraint_matrix <- matrix(c(2, 1, 1, 1, -1, 3), nrow = 2, byrow = TRUE)
rhs_constants <- c(2, 3)
constraint_directions <- c('<', '<')

solution <- lp(direction = 'min',
                objective.in = obj_coeffs,
                const.mat = constraint_matrix,
                const.dir = constraint_directions,
                const.rhs = rhs_constants)

# 输出结果
cat('目标函数的最优值为：', solution$objval,'\n')
cat('x的值为：', solution$solution[1],'\n')
cat('y的值为：', solution$solution[2],'\n')
cat('z的值为：', solution$solution[3],'\n')

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
# for loops
fl_1 <- list()
for (i in 1:length(formulas)) {
  # 在每次循环中，使用当前公式拟合线性模型，并将结果存储到列表中
  fl_1[[i]] <- lm(formulas[[i]], data = mtcars)
}
fl_1

## -----------------------------------------------------------------------------
# lapply
la_1 <- lapply(formulas, lm, data = mtcars)
la_1

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

## -----------------------------------------------------------------------------
# for loops
fl_2 <- list()
for (i in seq_along(bootstraps)) {
  fl_2[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
fl_2

## -----------------------------------------------------------------------------
# lapply
la_2 <- lapply(bootstraps, lm, formula = mpg ~ disp)
la_2

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
# Exercise 3
sapply(fl_1, rsq) %>%
  cat('For loops model:\n\t', ., '\n')
sapply(la_1, rsq) %>%
  cat('lapply model:\n\t', ., '\n')

## -----------------------------------------------------------------------------
# Exercise 4
sapply(fl_2, rsq) %>%
  cat('For loops model:\n\t', ., '\n')
sapply(la_2, rsq) %>%
  cat('lapply model:\n\t', ., '\n')

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

## -----------------------------------------------------------------------------
# anonymous function:
sapply(trials, function(x) x$p.value)

## -----------------------------------------------------------------------------
# without anonymous function:
sapply(trials, '[[', 'p.value')

## -----------------------------------------------------------------------------
# 自写lapply()
my_lapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  if (simplify == TRUE) return(simplify2array(out))
  else return(Map(function(x) vapply(x, FUN, FUN.VALUE), X))
}

# 做一个简单测试
list_test <- list(cars, mtcars)
my_lapply(list_test, mean, numeric(1))

## -----------------------------------------------------------------------------
# 自写chisq.test()
my_chisq.test <- function(x, y) {
  # 检查输入是否为两个数值向量且无缺失值
  if (!is.numeric(x) ||!is.numeric(y) || any(is.na(x)) || any(is.na(y))) {
    stop('应该输入两个无缺失值的数值向量。')
  }

  m <- rbind(x, y)
  margin1 <- rowSums(m)
  margin2 <- colSums(m)
  n <- sum(m)
  me <- tcrossprod(margin1, margin2) / n

  # 计算卡方检验统计量
  chi_square_stat <- sum((m - me)^2 / me)
  # 计算自由度
  degrees_of_freedom <- (length(margin1) - 1) * (length(margin2) - 1)
  # 计算p值
  p_value <- pchisq(chi_square_stat, df = degrees_of_freedom, lower.tail = FALSE)

  # 返回卡方检验统计量、自由度和p值组成的列表
  return(list(chi_square_statistic = chi_square_stat,
              degrees_of_freedom = degrees_of_freedom,
              p_value = p_value))
}

## ----results = ''-------------------------------------------------------------
# 测试
set.seed(123)
x <- sample(1:10, 100, replace = TRUE)
y <- sample(1:10, 100, replace = TRUE)
m <- cbind(x, y)

# chisq.test()结果
chisq.test(m)
# 自写chisq.test()结果
my_chisq.test(x, y)
# 性能评估和比较
suppressWarnings(
  bench::mark(
    chisq.test(m), 
    my_chisq.test(x, y), 
    check = FALSE
  )
)

## -----------------------------------------------------------------------------
# 自写table()
my_table <- function(x, y) {
  # 检查输入是否为两个整数型向量且无缺失值
  if (!is.integer(x) ||!is.integer(y) || any(is.na(x)) || any(is.na(y))) {
    stop('应该输入两个无缺失值的整数型向量。')
  }
  
  m <- cbind(x, y)
  unique_x <- sort(unique(x))
  unique_y <- sort(unique(y))
  
  num_unique_x <- length(unique_x)
  num_unique_y <- length(unique_y)
  
  table_dimensions  <- c(num_unique_x, num_unique_y)
  combinations_total <- num_unique_x * num_unique_y
  dimension_names <- list(x = unique_x, y = unique_y)
  
  bin <- fastmatch::fmatch(x, unique_x) +
    num_unique_x * fastmatch::fmatch(y, unique_y) - num_unique_x
  count_vector <- tabulate(bin, combinations_total)
  
  count_array  <- array(count_vector, dim = table_dimensions, dimnames = dimension_names)
  class(count_array) <- 'table'
  
  return(count_array)
}

## ----results = ''-------------------------------------------------------------
# 测试
set.seed(123)
x <- sample(1:10, 1000, replace = TRUE)
y <- sample(1:10, 1000, replace = TRUE)

# 比较二者是否一致
identical(table(x, y), my_table(x, y))

## -----------------------------------------------------------------------------
# 性能评估和比较
bench::mark(table(x, y), my_table(x, y), check = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# #include <Rcpp.h>
# using namespace Rcpp;
# 
# // 计算二项式系数
# int binomialCoefficient(int n, int k) {
#     if (k == 0 || k == n) return 1;
#     return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
# }
# 
# // Gibbs
# // [[Rcpp::export]]
# List gibbsSampler(int n, int a, int b, int m) {
#     // 存储生成的样本
#     NumericMatrix xSamples(m, 1);
#     NumericMatrix ySamples(m, 1);
#     // 初始化 x 和 y
#     int x = 0;
#     double y = 0.5;
#     for (int i = 0; i < m; ++i) {
#         x = R::rbinom(n, y);
#         y = R::rbeta(x + a, n - x + b);
#         xSamples(i, 0) = x;
#         ySamples(i, 0) = y;
#     }
#     return List::create(Named("xSamples") = xSamples,
#                         Named("ySamples") = ySamples);
# }
# 

## ----eval = FALSE-------------------------------------------------------------
# sourceCpp('Gibbs.cpp')
# 
# ## 设置参数
# n <- 10
# a <- 5
# b <- 8
# m <- 1000
# 
# ## Cpp结果
# results <- gibbsSampler(n, a, b, m)
# x_samples <- results$xSamples
# y_samples <- results$ySamples
# head(x_samples)
# head(y_samples)

## ----eval = FALSE-------------------------------------------------------------
# set.seed(123)
# m <- 1000
# 
# ## 使用R自带函数生成随机数
# Gibbs_R <- function(n, a, b, m) {
#   x_samples_r <- numeric(m)
#   y_samples_r <- numeric(m)
# 
#   # 初始化y
#   y_temp <- runif(1)
# 
#   for (i in 1:m) {
#     x_samples_r[i] <- rbinom(1, n, y_temp)
#     y_temp <- rbeta(1, x_samples_r[i] + a, n - x_samples_r[i] + b)
#     y_samples_r[i] <- y_temp
#   }
# 
#   return(list(x_samples = x_samples_r, y_samples = y_samples_r))
# }
# 
# ## 调用Gibbs_R生成随机数
# result_R <- Gibbs_R(n, a, b, m)
# x_r_samples <- result_R$x_samples
# y_r_samples <- result_R$y_samples
# 
# ## 画图
# qqplot(x_r_samples, x_samples,
#        xlab = 'Random x generated by R built-in function',
#        ylab = 'Random x generated by Rcpp function',
#        main = 'QQ plot of x random numbers')
# 
# qqplot(y_r_samples, y_samples,
#        xlab = 'Random y generated by R built-in function',
#        ylab = 'Random y generated by Rcpp function',
#        main = 'QQ plot of y random numbers')

## ----eval = FALSE-------------------------------------------------------------
# ## 比较
# microbenchmark(
#     Rcpp_function = gibbsSampler(n, a, b, m),
#     R_function = Gibbs_R(n, a, b, m)
# ) %>%
#   summary()

