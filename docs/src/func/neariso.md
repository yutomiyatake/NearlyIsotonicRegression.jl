# Nearly isotonic regression

Suppose that we have ``n`` independent normal observations ``X_i \sim \mathrm{N}(\mu_i,1)`` for ``i=1,\dots,n``, where ``\mu_1,\mu_2,\dots,\mu_n`` is a **piecewize** monotone sequence of normal means.
The **nearly isotonic regression problem** is formulated as
```math
\begin{aligned}
\hat{\mu} &=\argmin_{\mu_1, \dots, \mu_n}\left(\sum_{i=1}^n -\log p(x_i\mid \mu_i) +\lambda\sum_{i=1}^{n-1} (\mu_i-\mu_{i+1})_+\right)\\
&= \argmin_{\mu_1,\dots , \mu_n}\left(\sum_{i=1}^n\frac{1}{2} (x_i-\mu_i)^2+\lambda\sum_{i=1}^{n-1} (\mu_i-\mu_{i+1})_+\right),
\end{aligned}
```
where ``(a)_+=\mathrm{max}(a,0)`` and ``\lambda>0`` is a regularization parameter.

The weighted formulation reads
```math
\hat{\mu} 
=\argmin_{\mu_1, \dots, \mu_n} \left( \sum_{i=1}^n \frac{1}{2} w_i (x_i-\mu_i)^2+\lambda\sum_{i=1}^{n-1} (\mu_i-\mu_{i+1})_+\right)
```
where ``w_i > 0`` ``(i=1,\dots,n)``.


The function ```neariso``` solves these problems.
Note that the result ``\hat{\mu} = (\hat{\mu}_1,\dots, \hat{\mu}_n)`` consists of several clusters. 
This function outputs ``\hat{\mu}`` and the number of clusters.

The descriptive public names are `nearly_isotonic_regression` and
`nearly_isotonic_path`; the shorter names remain available for compatibility.

```@docs
nearly_isotonic_regression
nearly_isotonic_path
```

```@docs
neariso
```

The number of clusters decreases if one chooses a bigger ``\lambda``.
The function ```neariso_path``` outputs the ``\lambda`` values at which several clusters merge, and the corresponding number of clusters.

```@docs
neariso_path
```

# Generalized nearly isotonic regression


The above formulation can be generalized to one parameter exponential families, as is done for the isotonic regression.
Such a formulation is called the **generalized nearly isotonic regression**.

When ``X_i\sim p_i(x_i\mid \theta_i)`` for ``i=1,\dots,n``, where ``p_i(x_i\mid\theta_i )`` is a one parameter exponential family defined by
 ```math
p_i(x_i\mid \theta_i) = h_i(x_i) \exp (\theta_i x_i -w_i\psi(\theta_i)),
```
the problem 
```math
\hat{\theta} = \argmin_{\theta_1, \dots , \theta_n} \left(\sum_{i=1}^n -\log p_i(x_i\mid \theta_i)+\lambda\sum_{i=1}^{n-1} (\theta_i-\theta_{i+1})_+\right)
```
can also be solved by the function ```neariso```, as the formulation is equivalent to
```math
\hat{\eta} = \argmin_{\eta_1,\dots , \eta_n}\left(\sum_{i=1}^n\frac{1}{2} w_i (x_i-\eta_i)^2+\lambda\sum_{i=1}^{n-1} (\eta_i-\eta_{i+1})_+\right) \quad \text{with} \quad \eta = \mathrm{E}_\theta [x] = \psi ^\prime (\theta).
```

Further, apart from ```neariso_DistibutionName``` and ```neariso_path_DistributionName```, the function ```neariso_AIC_value_DistributionName``` provides the AIC value at ``\lambda``.
Additionally, ```neariso_AIC_DistributionName``` outouts the collection of the AIC values at each checkpoint of the relaxation parameter.
Here, the Akaike Information Criterion (AIC) is defined by
```math
\mathrm{AIC}(\lambda) = -2 \sum_{i=1}^n \log p_i(x_i\mid (\hat{\theta}_\lambda)_i) + 2K_\lambda,
```
where ``\hat{\theta}_\lambda = (\hat{\theta}_1,\dots,\hat{\theta}_n)`` is the solution for ``\lambda``, and ``K_\lambda`` signifies the corresponding number of clusters.


---
```@docs
neariso_Normal
```

```@docs
neariso_path_Normal
```

```@docs
neariso_AIC_Normal
```

```@docs
neariso_AIC_value_Normal
```


---
```@docs
neariso_Binomial
```

```@docs
neariso_path_Binomial
```

```@docs
neariso_AIC_Binomial
```

```@docs
neariso_AIC_value_Binomial
```

---
```@docs
neariso_Poisson
```

```@docs
neariso_path_Poisson
```

```@docs
neariso_AIC_Poisson
```

```@docs
neariso_AIC_value_Poisson
```

---
```@docs
neariso_Chisq
```

```@docs
neariso_path_Chisq
```

```@docs
neariso_AIC_Chisq
```

```@docs
neariso_AIC_value_Chisq
```
---
