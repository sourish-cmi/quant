library(quantmod)
library(tseries)
library(evd)  # for rgpd
library(mvtnorm)
library(ggplot2)
library(zoo)
library(dplyr)
library(xts)


source("functions_v1.R")

## Download data from yahoo

NSEI = get.hist.quote("^NSEI")
head(NSEI)
n=nrow(NSEI)
NSEI$rt[2:n] = diff(log(NSEI$Close))*100

NSEI$emp_vol <- NA
alpha <- 0.9
k=21
for(i in k:n){
  sigma_sqr <- var(NSEI$rt[(i-k+1):(i-1)],na.rm=TRUE)
  NSEI$emp_vol[i]<-alpha*sigma_sqr+(1-alpha)*NSEI$rt[i]^2
  NSEI$emp_vol[i]<-sqrt(NSEI$emp_vol[i])*sqrt(250)
}
tail(NSEI)

NSEI$gk_vol <- mapply(gk_volatility,
                 open = NSEI$Open,
                 high = NSEI$High,
                 low = NSEI$Low,
                 close = NSEI$Close)

NSEI$gk_vol<- sqrt(NSEI$gk_vol)*sqrt(250)*100


# Assuming NSEI is a data frame with columns emp_vol and gk_vol
# Step 1: Convert zoo or xts object to data frame with proper Date column
NSEI_df <- data.frame(Date = index(NSEI), coredata(NSEI))
save(NSEI_df,file="NSEI.RData")

##load data
load("NSEI.RData",verbose = T)

## Download data from yahoo

GSPC = get.hist.quote("^GSPC")

GOLD = get.hist.quote("GC=F")
# Convert GSPC to data frame with Date
GSPC_df <- GSPC
GSPC_df <- data.frame(Date = index(GSPC), coredata(GSPC))

GOLD_df <- GOLD
GOLD_df <- data.frame(Date = index(GOLD), coredata(GOLD))


# Step 1: Calculate log returns
GSPC_df <- GSPC_df %>%
  arrange(Date) %>%
  mutate(rt_snp500 = 100 * c(NA, diff(log(Close))))

GOLD_df <- GOLD_df %>%
  arrange(Date) %>%
  mutate(rt_gold = 100 * c(NA, diff(log(Close))))


# Step 2: Empirical volatility using EWMA
alpha <- 0.9
k <- 21
n <- nrow(GSPC_df)
GSPC_df$emp_vol_snp500 <- NA
for (i in k:n) {
  sigma_sqr <- var(GSPC_df$rt_snp500[(i - k + 1):(i - 1)], na.rm = TRUE)
  GSPC_df$emp_vol_snp500[i] <- sqrt(alpha * sigma_sqr + (1 - alpha) * GSPC_df$rt_snp500[i]^2) * sqrt(250)
}

alpha <- 0.9
k <- 21
n <- nrow(GOLD_df)
GOLD_df$emp_vol_gold <- NA
for (i in k:n) {
  sigma_sqr <- var(GOLD_df$rt_gold[(i - k + 1):(i - 1)], na.rm = TRUE)
  GOLD_df$emp_vol_gold[i] <- sqrt(alpha * sigma_sqr + (1 - alpha) * GOLD_df$rt_gold[i]^2) * sqrt(250)
}


# Step 3: Garman-Klass volatility
GSPC_df$gk_vol_snp500 <- mapply(gk_volatility, GSPC_df$Open, GSPC_df$High, GSPC_df$Low, GSPC_df$Close)
GSPC_df$gk_vol_snp500 <- sqrt(GSPC_df$gk_vol_snp500)*sqrt(250)*100

GOLD_df$gk_vol_gold <- mapply(gk_volatility, GOLD_df$Open, GOLD_df$High, GOLD_df$Low, GOLD_df$Close)
GOLD_df$gk_vol_gold <- sqrt(GOLD_df$gk_vol_gold)*sqrt(250)*100

# Step 4: Rename Close as snp500, drop Open, High, Low
GSPC_df <- GSPC_df %>%
  rename(snp500 = Close) %>%
  select(Date, snp500, rt_snp500, emp_vol_snp500, gk_vol_snp500)

GOLD_df <- GOLD_df %>%
  rename(gold = Close) %>%
  select(Date, gold, rt_gold, emp_vol_gold, gk_vol_gold)

# Step 5: Merge with NSEI_df by Date

merged_df <- merge(NSEI_df, GSPC_df, by = "Date")
merged_df <- merge(merged_df, GOLD_df, by = "Date")

# View merged data
head(merged_df)

save(merged_df,file="NSEI.RData")
