# NearlyIsotonicRegression.jl

[![CI](https://github.com/yutomiyatake/NearlyIsotonicRegression.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/yutomiyatake/NearlyIsotonicRegression.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://yutomiyatake.github.io/NearlyIsotonicRegression.jl/dev/)

NearlyIsotonicRegression.jl implements isotonic and nearly isotonic regression,
including weighted estimators and helpers for normal, binomial, Poisson, and
chi-squared models. See the
[documentation](https://yutomiyatake.github.io/NearlyIsotonicRegression.jl/dev/)
for the model definitions and complete API.

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

The descriptive API consists of `isotonic_regression`,
`isotonic_regression!`, `nearly_isotonic_regression`, and
`nearly_isotonic_path`. The original short names `iso`, `iso!`, `neariso`, and
`neariso_path` remain available. Mutating functions require floating-point
input so that pooled estimates can be written back to the supplied vector.


## References

- R. E. Barlow, D. J. Bartholomew, J. M. Bremner and H. D. Brunk: *Statistical Inference Under Order Restrictions* (1972)
- P. Groeneboom and G. Jongbloed: *Nonparametric Estimation Under Shape Constraints* (2014)
- T. Matsuda and Y. Miyatake: “Piecewise monotone estimation in one-parameter exponential families” (2025)
- T. Robertson, F. T. Wright and R. L. Dykstra: *Order Restricted Statistical Inference* (1988)
- R. J. Tibshirani, H. Hoefling and R. Tibshirani: “Nearly-isotonic regression” (2011)
- C. van Eeden: *Restricted Parameter Space Estimation Problems* (2006)

## Development note

The v1.0.0 registration-preparation work used OpenAI Codex for code review,
test expansion, and maintenance-file updates. The maintainer is responsible for
reviewing these changes before release.
