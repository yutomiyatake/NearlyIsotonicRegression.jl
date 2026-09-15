```@meta
CurrentModule = NearlyIsotonicRegression
```

# NearlyIsotonicRegression.jl

NearlyIsotonicRegression.jl implements isotonic and nearly isotonic regression:

- standard and weighted isotonic regression
- generalized isotonic regression
- standard and weighted nearly isotonic regression
- generalized nearly isotonic regression

## Installation

From the Julia package prompt, run:

```julia
pkg> add NearlyIsotonicRegression
```

## Simple examples

```julia
using NearlyIsotonicRegression

n = 10
x = (1:n)/n + randn(n)
y = isotonic_regression(x)
```

The original short names `iso`, `iso!`, `neariso`, and `neariso_path` remain
available for compatibility.


## References

- R. E. Barlow, D. J. Bartholomew, J. M. Bremner and H. D. Brunk: *Statistical Inference Under Order Restrictions* (1972)
- P. Groeneboom and G. Jongbloed: *Nonparametric Estimation Under Shape Constraints* (2014)
- T. Matsuda and Y. Miyatake: “Piecewise monotone estimation in one-parameter exponential families” (2025)
- T. Robertson, F. T. Wright and R. L. Dykstra: *Order Restricted Statistical Inference* (1988)
- R. J. Tibshirani, H. Hoefling and R. Tibshirani: “Nearly-isotonic regression” (2011)
- C. van Eeden: *Restricted Parameter Space Estimation Problems* (2006)
