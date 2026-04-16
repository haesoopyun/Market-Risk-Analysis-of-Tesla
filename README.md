# Market Risk Analysis of Tesla (TSLA)

1-day 99% Value at Risk (VaR) and Expected Shortfall (ES) for Tesla stock and an ATM call option, implemented in both **R** and **Python**.

## Methods

- **Variance-Covariance (VCV)** with unconditional and EWMA(0.94) conditional volatility
- **Historical Simulation** with standardised residuals (unconditional and conditional)
- **Monte Carlo Simulation** (200,000 GBM paths) with conditional volatility
- **Extreme Value Theory** (Peaks Over Threshold, GPD fit)
- **Option full re-pricing** via Black-Scholes under Monte Carlo scenarios

## Repository Structure

```
Market-Risk-Analysis-of-Tesla/
├── src/
│   ├── market_risk_analysis.r        # R implementation
│   ├── market_risk_analysis.ipynb    # Python (Jupyter) implementation
│   └── TSLA_Data.csv                 # Historical price data
├── R_figures/                        # Figures and tables from R
├── Python_figures/                   # Figures and tables from Python
├── report/
│   └── Market Risk Analysis of Tesla-Report.pdf
├── .gitignore
└── README.md
```

## Running

### R

```r
# From the src/ directory
setwd("src")
source("market_risk_analysis.r")
```

Requires: `data.table`, `xts`, `ggplot2`, `timeDate`, `lubridate`, `evir`, `scales`, `moments`, `knitr`, `kableExtra`, `gridExtra`

### Python

```bash
# From the src/ directory
jupyter notebook market_risk_analysis.ipynb
```

Requires: `numpy`, `pandas`, `yfinance`, `scipy`, `matplotlib`, `seaborn`

## Author

Haesoo Pyun
