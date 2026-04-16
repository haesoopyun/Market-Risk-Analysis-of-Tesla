# Market Risk Analysis of Tesla
# Notes:
# - All VaR/ES reported as positive loss magnitudes in simple return terms.
# - 99% one-day horizon
# - Mean daily simple return assumed 0 in VaR/ES models.

req <- c("data.table","xts","ggplot2","timeDate","lubridate",
         "evir","scales","moments","knitr","kableExtra")
missing <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Install them and re-run.")
}

library(data.table); library(xts); library(ggplot2)
library(timeDate);  library(lubridate); library(evir)
library(scales);    library(moments);   library(knitr); library(kableExtra)
library(gridExtra)
theme_rm <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title  = element_text(face = "bold", hjust = 0),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(colour = "black")
  )

alpha <- 0.01              # 99% VaR/ES
lambda_ewma <- 0.94        # EWMA decay
rf_1m <- 0.04364           # 1m T-bill simple annual rate
rf_cont <- log(1 + rf_1m)  # continuous comp
trading_days <- 252
as_of <- as.Date("2025-08-04")
set.seed(50)

fig_dir <- file.path("..", "R_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

dt <- data.table::fread(file = "TSLA_Data.csv")

print(head(dt))
print(names(dt))

data.table::setnames(dt, old = "tsla", new = "adj_close", skip_absent = TRUE)
dt[, date := as.Date(
  as.character(date),
  tryFormats = c("%d/%m/%Y", "%Y-%m-%d")
)]

dt_hist <- dt[date <= as_of & adj_close > 0]

S_t <- dt_hist[.N, adj_close]

# Historical daily returns
dt_use <- data.table::copy(dt_hist)
dt_use[, r_simple := adj_close / shift(adj_close) - 1]
dt_use[, r_log    := log(adj_close) - log(shift(adj_close))]
dt_use <- dt_use[!is.na(r_simple)]
dt_use[, date := as.Date(date)]
str(dt_use$date)

desc <- dt_use[, .(
  n = .N,
  mean_simple = mean(r_simple),
  sd_simple   = sd(r_simple),
  skew_simple = moments::skewness(r_simple),
  kurt_simple = moments::kurtosis(r_simple),
  min_simple  = min(r_simple),
  max_simple  = max(r_simple)
)]

plt_price <- ggplot(dt_use, aes(x = date, y = adj_close)) +
  geom_line(linewidth = 0.4) +
  labs(
    title = "TSLA Adjusted Close Price",
    x = "Date",
    y = "Price (USD)"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  theme_rm

print(plt_price)
ggsave(file.path(fig_dir, "fig_price.png"),   plt_price, width = 7, height = 4, dpi = 300)

plt_ret <- ggplot(dt_use, aes(x = date, y = r_simple)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(colour = "grey40", linewidth = 0.25) +
  labs(
    title = "TSLA Daily Simple Returns",
    x = "Date",
    y = "Return"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_rm

print(plt_ret)
ggsave(file.path(fig_dir, "fig_return.png"),   plt_ret, width = 7, height = 4, dpi = 300)

# EWMA variance recursion on simple returns
ewma_var <- function(r, lambda, v0 = var(r, na.rm = TRUE)){
  out <- numeric(length(r)); out[1] <- v0
  for(i in 2:length(r)) out[i] <- lambda*out[i-1] + (1-lambda)*r[i-1]^2
  out
}

# Black–Scholes call
pn <- pnorm
bs_call <- function(S, K, r_cont, sigma_a, T_years){
  if(T_years <= 0) return(pmax(S - K, 0))
  d1 <- (log(S/K) + (r_cont + 0.5*sigma_a^2)*T_years) / (sigma_a*sqrt(T_years))
  d2 <- d1 - sigma_a*sqrt(T_years)
  S*pn(d1) - K*exp(-r_cont*T_years)*pn(d2)
}

# Volatility (unconditional & EWMA)
# Unconditional sigma
sigma_uncond <- sd(dt_use$r_simple)

# Conditional sigma_t via EWMA and 1-step forecast
var_ewma <- ewma_var(dt_use$r_simple, lambda_ewma, v0 = sigma_uncond^2)
sigma_ewma <- sqrt(var_ewma)
sigma_cond_forecast <- tail(sigma_ewma, 1)

# 99% one day normal VaR and ES series on return scale
z_99 <- qnorm(alpha)
VaR_cond_line <- z_99 * sigma_ewma
ES_cond_mag   <- sigma_ewma * dnorm(z_99) / alpha
ES_cond_line  <- -ES_cond_mag

plt_cond <- ggplot(
  data = data.table(
    date          = dt_use$date,
    r_simple      = dt_use$r_simple,
    VaR_cond_line = VaR_cond_line,
    ES_cond_line  = ES_cond_line
  ),
  aes(x = date)
) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(aes(y = r_simple), colour = "grey40", linewidth = 0.25) +
  geom_line(aes(y = VaR_cond_line), colour = "red", linewidth = 0.4) +
  geom_line(aes(y = ES_cond_line),  colour = "blue", linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "TSLA Returns with 1-day 99% Conditional VaR and ES (VCV, EWMA σ)",
    x = "Date",
    y = "Return"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_rm

print(plt_cond)

ggsave(file.path(fig_dir, "fig_condVaR.png"), plt_cond,  width = 7, height = 4, dpi = 300)

# VCV VaR & ES (uncond vs cond, Normal)
qz <- qnorm(alpha)
VaR_vcv_uncond <- - sigma_uncond * qz
ES_vcv_uncond  <-   sigma_uncond * dnorm(qz) / alpha
VaR_vcv_cond   <- - sigma_cond_forecast * qz
ES_vcv_cond    <-   sigma_cond_forecast * dnorm(qz) / alpha

dt_use[, r_std := (r_simple - mean(r_simple)) / sd(r_simple)]
plt_tail <- ggplot(dt_use, aes(sample = r_std)) +
  stat_qq(distribution = qnorm, linewidth = 0.3) +
  stat_qq_line(distribution = qnorm, colour = "red", linewidth = 0.4) +
  labs(
    title = "Q-Q plot: TSLA Returns vs Standard Normal",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_rm

print(plt_tail)
ggsave(file.path(fig_dir, "fig_QQ.png"),plt_tail,  width = 7, height = 4, dpi = 300)

# Standardised Historical Simulation
mu_sample <- mean(dt_use$r_simple)
z <- (dt_use$r_simple - mu_sample) / sigma_ewma
z_valid <- z[max(1, 251):length(z)]

q_std <- as.numeric(quantile(z_valid, probs = alpha, na.rm = TRUE))
ES_std <- mean(z_valid[z_valid <= q_std], na.rm = TRUE)

VaR_hs_uncond <- - q_std * sigma_uncond
ES_hs_uncond  <- - ES_std * sigma_uncond
VaR_hs_cond   <- - q_std * sigma_cond_forecast
ES_hs_cond    <- - ES_std * sigma_cond_forecast

z_df <- data.table(date = dt_use$date, z = z)
z_df <- z_df[251:.N]
plt_z <- ggplot(z_df, aes(x = z)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = "grey85",
    colour = "white"
  ) +
  stat_function(
    fun  = dnorm,
    args = list(mean = 0, sd = 1),
    linewidth = 0.5,
    colour = "red"
  ) +
  labs(
    title = "Standardised Returns z_t after EWMA (first 250 dropped)",
    x = "z_t",
    y = "Density"
  ) +
  theme_rm

print(plt_z)
ggsave(file.path(fig_dir, "fig_std_hist.png"),plt_z,     width = 7, height = 4, dpi = 300)

# Monte Carlo (GBM) for stock
# Choose volatility for GBM
sigma_choice <- "cond"
sigma_daily <- switch(sigma_choice,
                      "uncond" = sigma_uncond,
                      "cond"   = sigma_cond_forecast,
                      "implied"= 0.46 / sqrt(trading_days))

nsims <- 200000
z_mc <- rnorm(nsims)
r1_mc <- (rf_cont/trading_days - 0.5*sigma_daily^2) + sigma_daily*z_mc
R1_mc <- exp(r1_mc) - 1

VaR_mc_stock <- - as.numeric(quantile(R1_mc, alpha))
ES_mc_stock  <- - mean(R1_mc[R1_mc <= -VaR_mc_stock])

mc_df <- data.table(R1_mc = R1_mc)
mu_mc <- mean(mc_df$R1_mc)
sd_mc <- sd(mc_df$R1_mc)

plt_mc <- ggplot(mc_df, aes(x = R1_mc)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins  = 60,
    fill  = "grey",
    colour = "white",
    alpha = 0.8
  ) +
  stat_function(
    fun  = function(x) dnorm(x, mean = mu_mc, sd = sd_mc),
    linewidth = 0.7
  ) +
  # VaR marker in the left tail
  geom_vline(
    xintercept = -VaR_mc_stock,
    colour     = "red",
    linewidth  = 0.8,
    linetype   = "dashed"
  ) +
  labs(
    title = "Monte Carlo (GBM): 1-day Return Distribution for TSLA",
    x     = "1-day Simple Return",
    y     = "Density"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_rm

print(plt_mc)
ggsave(file.path(fig_dir, "fig_mc.png"),      plt_mc,    width = 7, height = 4, dpi = 300)

# POT (GPD) for stock
# Losses L = -r. Threshold at 95th percentile of losses.
L <- as.numeric(-dt_use$r_simple)
L <- L[is.finite(L)]
u <- as.numeric(quantile(L, 0.95, na.rm = TRUE))

# Fit GPD on losses above u using evir.
fit <- evir::gpd(L, threshold = u)
xi   <- as.numeric(fit$par.ests["xi"])
beta <- as.numeric(fit$par.ests["beta"])
p_u  <- mean(L > u)

# VaR and ES in loss terms (positive), with xi -> 0 guard
eps <- 1e-6
if (abs(xi) < eps) {
  # Limit xi -> 0: q_alpha = u + beta * log(p_u/alpha)
  VaR_pot_loss <- u + beta * log(p_u/alpha)
  ES_pot_loss  <- VaR_pot_loss + beta
} else {
  VaR_pot_loss <- u + (beta/xi) * ((alpha/p_u)^(-xi) - 1)
  ES_pot_loss  <- if (xi < 1) (VaR_pot_loss + (beta - xi*u)) / (1 - xi) else NA_real_
}
VaR_pot <- VaR_pot_loss
ES_pot  <- ES_pot_loss


png(file.path(fig_dir, "fig_pot_diagnostics.png"), width = 1400, height = 600, res = 150)

par(mfrow = c(1, 2)) 

# 1. Mean excess plot
evir::meplot(L)

# 2. QQ plot for GPD fit
excesses <- L[L > u] - u
n_ex <- length(excesses)
p_seq <- (1:n_ex) / (n_ex + 1)
if (abs(xi) < 1e-6) {
  q_theoretical <- -beta * log(1 - p_seq)
} else {
  q_theoretical <- (beta/xi) * ((1 - p_seq)^(-xi) - 1)
}
qqplot(
  q_theoretical, sort(excesses),
  main = "GPD QQ Plot (Empirical vs Theoretical)",
  xlab = "Theoretical Quantiles",
  ylab = "Empirical Quantiles",
  pch = 19,
  col = "blue"
)
abline(0, 1, col = "red", lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# Option VaR & ES (1-day MC then BS reprice)
opt_grid <- data.table(
  strike = c(305, 310, 315, 320),
  bid    = c(19.65, 17.05, 14.75, 12.65),
  ask    = c(19.80, 17.20, 14.85, 12.80),
  iv     = 0.46
)
k_row <- opt_grid[which.min(abs(strike - S_t))]
X <- k_row$strike
C_mid <- (k_row$bid + k_row$ask)/2
iv_a <- k_row$iv

# One day ahead date and time to expiry
t1 <- as.Date("2025-08-05"); t_exp <- as.Date("2025-09-05")
T_years_toX <- as.numeric(t_exp - t1) / 365

# Simulate S_{t+1} with chosen stock sigma (same as Task 4)
S1 <- S_t * exp( (rf_cont/trading_days - 0.5*sigma_daily^2) + sigma_daily*rnorm(nsims) )

# Reprice call with Black–Scholes using implied vol
C1 <- bs_call(S = S1, K = X, r_cont = rf_cont, sigma_a = iv_a, T_years = T_years_toX)
R_call <- (C1 - C_mid)/C_mid
VaR_call_mc <- - as.numeric(quantile(R_call, alpha))
ES_call_mc  <- - mean(R_call[R_call <= -VaR_call_mc])

call_df <- data.table(R_call = R_call)
plt_call <- ggplot(call_df, aes(x = R_call)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = "grey85",
    colour = "white"
  ) +
  geom_vline(
    xintercept = -VaR_call_mc,
    colour = "red",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  labs(
    title = "ATM Call Option: 1-day Return Distribution",
    x = "Option Simple Return",
    y = "Density"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_rm

print(plt_call)
ggsave(file.path(fig_dir, "fig_call.png"),    plt_call,  width = 7, height = 4, dpi = 300)

# Presentation
tab_stock <- data.table(
  method    = c("VCV uncond","VCV cond","HS uncond","HS cond","MC (GBM)","POT (GPD)"),
  VaR_99_1d = c(VaR_vcv_uncond, VaR_vcv_cond, VaR_hs_uncond, VaR_hs_cond, VaR_mc_stock, VaR_pot),
  ES_99_1d  = c(ES_vcv_uncond, ES_vcv_cond, ES_hs_uncond, ES_hs_cond, ES_mc_stock, ES_pot)
)

tab_option <- data.table(
  S_t = S_t, strike = X, mid = C_mid, iv_annual = iv_a,
  VaR_99_1d = VaR_call_mc, ES_99_1d = ES_call_mc
)

# tables
desc_show <- copy(desc)
desc_show[, c("mean_simple","sd_simple","skew_simple","kurt_simple","min_simple","max_simple")
          := lapply(.SD, function(v) round(v, 4)),
          .SDcols = c("mean_simple","sd_simple","skew_simple","kurt_simple","min_simple","max_simple")]

kable(desc_show, caption = "Descriptive statistics of daily simple returns") %>%
  kable_styling(full_width = FALSE)

tab_stock_show <- copy(tab_stock)
tab_stock_show[, c("VaR_99_1d","ES_99_1d") := lapply(.SD, function(v) round(100*v, 2)),
               .SDcols = c("VaR_99_1d","ES_99_1d")]
setnames(tab_stock_show, c("VaR_99_1d","ES_99_1d"), c("VaR 99% (%, 1d)","ES 99% (%, 1d)"))

kable(tab_stock_show, caption = "Tesla stock: 99% 1-day VaR and ES across methods") %>%
  kable_styling(full_width = FALSE)

tab_option_show <- copy(tab_option)
tab_option_show[, c("VaR_99_1d","ES_99_1d") := lapply(.SD, function(v) round(100*v, 1)),
                .SDcols = c("VaR_99_1d","ES_99_1d")]
setnames(tab_option_show, c("VaR_99_1d","ES_99_1d"), c("VaR 99% (%, 1d)","ES 99% (%, 1d)"))

kable(tab_option_show, caption = "ATM call option (repriced after 1 day): 99% VaR and ES") %>%
  kable_styling(full_width = FALSE)

# Excesses above threshold
  excesses <- L[L > u] - u
  n_ex     <- length(excesses)
  excess_sorted <- sort(excesses)
  
# Mean excess plot data
  mean_excess <- sapply(excess_sorted, function(t) mean(excesses[excesses > t]))
  
  me_df <- data.table(
    threshold   = excess_sorted,
    mean_excess = mean_excess
  )
  
  plt_me <- ggplot(me_df, aes(x = threshold, y = mean_excess)) +
    geom_line(linewidth = 0.4) +
    labs(
      title = "Mean Excess Plot for Losses Above u",
      x = "Excess Over Threshold u",
      y = "Mean Excess"
    ) +
    theme_rm
print(plt_me)
ggsave(file.path(fig_dir, "fig_Meanexcess.png"),    plt_me,  width = 7, height = 4, dpi = 300)

# GPD QQ plot data
# Fit GPD
fit <- evir::gpd(L, threshold = u)
xi   <- as.numeric(fit$par.ests["xi"])
beta <- as.numeric(fit$par.ests["beta"])

# Excesses above threshold
excesses <- L[L > u] - u
excess_sorted <- sort(excesses)
n_ex <- length(excess_sorted)

# Theoretical quantiles
p_seq <- (1:n_ex) / (n_ex + 1)
if (abs(xi) < 1e-6) {
  q_theoretical <- -beta * log(1 - p_seq)
} else {
  q_theoretical <- (beta/xi) * ((1 - p_seq)^(-xi) - 1)
}

# Data for QQ plot
qq_df <- data.table(
  theoretical = q_theoretical,
  empirical   = excess_sorted
)

# Plot
plt_pot_qq <- ggplot(qq_df, aes(x = theoretical, y = empirical)) +
  geom_point(size = 0.6, alpha = 0.7, colour = "blue") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 0.5) +
  labs(
    title = "GPD Q-Q Plot for Excess Losses",
    x = "Theoretical GPD Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_rm

print(plt_pot_qq)
ggsave(file.path(fig_dir, "fig_GPDQQ.png"),plt_pot_qq,  width = 7, height = 4, dpi = 300)
  
# Export tables
data.table::fwrite(desc,      file.path(fig_dir, "table_descriptives.csv"))
data.table::fwrite(tab_stock, file.path(fig_dir, "table_stock_VaR_ES.csv"))
data.table::fwrite(tab_option,file.path(fig_dir, "table_option_VaR_ES.csv"))

# Sensitivity analysis (MC volatility and EWMA lambda)

# Monte Carlo VaR/ES sensitivity to volatility input
sigma_daily_iv <- 0.46 / sqrt(trading_days)

z_mc_iv  <- rnorm(nsims)
r1_mc_iv <- (rf_cont/trading_days - 0.5 * sigma_daily_iv^2) + sigma_daily_iv * z_mc_iv
R1_mc_iv <- exp(r1_mc_iv) - 1

VaR_mc_stock_iv <- -as.numeric(quantile(R1_mc_iv, alpha))
ES_mc_stock_iv  <- -mean(R1_mc_iv[R1_mc_iv <= -VaR_mc_stock_iv])

sens_mc <- data.table(
  vol_type    = c("Conditional EWMA sigma", "Implied volatility 46%"),
  sigma_daily = c(sigma_cond_forecast,      sigma_daily_iv),
  VaR_99_1d   = c(VaR_mc_stock,             VaR_mc_stock_iv),
  ES_99_1d    = c(ES_mc_stock,              ES_mc_stock_iv)
)

sens_mc[, `:=`(
  sigma_daily = round(100 * sigma_daily, 3),
  VaR_99_1d   = round(100 * VaR_99_1d, 2),
  ES_99_1d    = round(100 * ES_99_1d, 2)
)]

print(sens_mc)

# VCV VaR sensitivity to EWMA decay factor lambda

lambda_grid <- c(0.90, 0.94, 0.97)
VaR_lambda  <- numeric(length(lambda_grid))
ES_lambda   <- numeric(length(lambda_grid))

for (i in seq_along(lambda_grid)) {
  var_ewma_tmp   <- ewma_var(dt_use$r_simple, lambda = lambda_grid[i], v0 = sigma_uncond^2)
  sigma_cond_tmp <- sqrt(tail(var_ewma_tmp, 1))
  
  VaR_lambda[i] <- -sigma_cond_tmp * qz
  ES_lambda[i]  <-  sigma_cond_tmp * dnorm(qz) / alpha
}

sens_lambda <- data.table(
  lambda    = lambda_grid,
  VaR_99_1d = VaR_lambda,
  ES_99_1d  = ES_lambda
)

sens_lambda[, `:=`(
  VaR_99_1d = round(100 * VaR_99_1d, 2),
  ES_99_1d  = round(100 * ES_99_1d, 2)
)]

print(sens_lambda)

if (exists("SAVE_PLOTS") && isTRUE(SAVE_PLOTS)) {
  
  sens_table <- rbind(
    data.table(
      Component     = "Monte Carlo (GBM)",
      Specification = "Conditional EWMA volatility",
      sigma_daily   = round(100 * sigma_cond_forecast, 3),
      VaR_99_1d     = round(100 * VaR_mc_stock, 2),
      ES_99_1d      = round(100 * ES_mc_stock, 2)
    ),
    data.table(
      Component     = "Monte Carlo (GBM)",
      Specification = "Implied volatility 46%",
      sigma_daily   = round(100 * sigma_daily_iv, 3),
      VaR_99_1d     = round(100 * VaR_mc_stock_iv, 2),
      ES_99_1d      = round(100 * ES_mc_stock_iv, 2)
    ),
    data.table(
      Component     = "VCV (EWMA λ sensitivity)",
      Specification = "λ = 0.90",
      sigma_daily   = NA_real_,
      VaR_99_1d     = sens_lambda$VaR_99_1d[1],
      ES_99_1d      = sens_lambda$ES_99_1d[1]
    ),
    data.table(
      Component     = "VCV (EWMA λ sensitivity)",
      Specification = "λ = 0.94 (baseline)",
      sigma_daily   = NA_real_,
      VaR_99_1d     = sens_lambda$VaR_99_1d[2],
      ES_99_1d      = sens_lambda$ES_99_1d[2]
    ),
    data.table(
      Component     = "VCV (EWMA λ sensitivity)",
      Specification = "λ = 0.97",
      sigma_daily   = NA_real_,
      VaR_99_1d     = sens_lambda$VaR_99_1d[3],
      ES_99_1d      = sens_lambda$ES_99_1d[3]
    )
  )
  
  print(sens_table)
  
  fwrite(sens_table, file.path(fig_dir, "table_sensitivity_analysis.csv"))
}

