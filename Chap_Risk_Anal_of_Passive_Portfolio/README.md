# Risk Analysis of Passive Portfolio

## Sourish Das

### **Abstract**

In this work, we present an alternative passive investment strategy. The passive investment philosophy comes from the Efficient Market Hypothesis (EMH), and its adoption is widespread. If EMH is true, one cannot outperform market by actively managing their portfolio for a long time. Also, it requires little to no intervention. People can buy an exchange-traded fund (ETF) with a long-term perspective. As the economy grows over time, one expects the ETF to grow. For example, in India, one can invest in NETF (see,\cite{NETF}), which suppose to mimic the Nifty50 return. However, the weights of the Nifty 50 index are based on market capitalisation. These weights are not necessarily optimal for the investor. In this work, we present that volatility risk and extreme risk measures of the Nifty50 portfolio are uniformly larger than Markowitz's optimal portfolio. However, common people can't create an optimised portfolio. So we proposed an alternative passive investment strategy of an equal-weight portfolio. We show that if one pushes the maximum weight of the portfolio towards equal weight, the idiosyncratic risk of the portfolio would be minimal. The empirical evidence indicates that the risk profile of an equal-weight portfolio is similar to that of Markowitz's optimal portfolio. Hence instead of buying Nifty50 ETFs, one should equally invest in the stocks of Nifty50 to achieve a uniformly better risk profile than the Nifty 50 ETF portfolio. We also presengt an analysis about how portfolios perform to idiosyncratic events like Russian invasion of Ukraine. We found that the equal weight portfolio has uniformly lower risk that Nifty 50 portfolio, before and during the Russia-Uknraine war.

## R-Code:

### Download data from Yahoo

```R
rm(list=ls())
library(quantmod)
library(tseries)

getSymbols("^NSEI",src = "yahoo")
Nifty50 = NSEI$NSEI.Adjusted
plot(Nifty50)

### log-return
log_return = diff(log(Nifty50))*100
n = length(log_return)

plot(log_return)
```
