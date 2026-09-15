# Isotonic regression

Suppose that we have ``n`` independent normal observations ``X_i \sim \mathrm{N}(\mu_i,1)`` for ``i=1,\dots,n``, where
```math
\mu_1\leq \mu_2 \leq \dots \leq \mu_n
```
is a monotone sequence of means of a normal distribution.
The maximum likelihood estimate (MLE) of ``\mu`` is the solution to the constrained optimization problem:
```math
\hat{\mu} = \argmax_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n\log p(x_i\mid \mu_i) = \argmin_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n \frac{1}{2} (x_i-\mu_i)^2.
```
This formulation is 
the **isotonic regression** of ``x_1,\dots,x_n`` with uniform weights (``w_i=1`` in the formulation below).

The weighted formulation reads
```math
\hat{\mu} = \argmin_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n \frac{1}{2} w_i (x_i-\mu_i)^2,
```
where ``w_i > 0`` ``(i=1,\dots,n)``.

The function ```iso``` or ```iso!``` solves these problems using the algorithm called PAVA (Pool-Adjacent-Violators algorithm).

The descriptive public names are `isotonic_regression` and
`isotonic_regression!`; the shorter names remain available for compatibility.

```@docs
isotonic_regression
isotonic_regression!
```

```@docs
iso
```

```@docs
iso!
```

# Generalized isotonic regression

The above formulation can be generalized to accommodate distributions other than the normal distribution with a uniform variance.
Such a formulation is often called the **generalized isotonic regression**.

When ``X_i\sim p_i(x_i\mid \theta_i)`` for ``i=1,\dots,n``, where ``p_i(x_i\mid\theta_i )`` is a one parameter exponential family defined by
 ```math
p_i(x_i\mid \theta_i) = h_i(x_i) \exp (\theta_i x_i -w_i \psi(\theta_i)),
```
the problem 
```math
\hat{\theta} = \argmax_{\theta_1\leq \dots \leq \theta_n} \sum_{i=1}^n\log p_i(x_i\mid \theta_i)
```
can also be solved by the function ```iso``` or ```iso!```, as the formulation is equivalent to
```math
\hat{\eta} = \argmin_{\eta_1\leq \dots \leq \eta_n} \sum_{i=1}^n \frac{w_i}{2} \Big(\frac{x_i}{w_i}-\eta_i\Big)^2 \quad \text{with} \quad \eta = \mathrm{E}_\theta [x] = \psi ^\prime (\theta).
```

Below, we list four examples of one parameter exponential families.
For each case, the relation between the natural parameter ``\theta_i`` and the expectation parameter ``\eta_i = \mathrm{E}_\theta [x_i] = \psi ^\prime (\theta_i)`` is specified. 

---
### Normal

The probability distribution function (PDF) of the **Normal Distribution** with mean ``\mu_i`` and standard deviation ``\sigma_i\geq 0`` is given by

```math
p_i(x_i\mid \mu_i) = \frac{1}{\sqrt{2\pi\sigma_i^2}} \exp \bigg( -\frac{(x_i-\mu_i)^2}{2\sigma_i^2} \bigg).
```

Since
```math
\begin{aligned}
\log p_i(x_i\mid \mu_i) &= \mu_i \frac{x_i}{\sigma_i^2} - \frac{\mu_i^2}{2\sigma_i^2} - \frac{x_i^2}{2\sigma_i^2} - \frac{1}{2}\log(2\pi\sigma_i^2) \\
&= \theta_i \tilde{x}_i - \sigma_i^{-2}\frac{\theta_i^2}{2} - 
\frac{\sigma_i^2 \tilde{x}_i^2}{2} - \frac{1}{2}\log(2\pi\sigma_i^2) \quad \bigg(\theta_i := \frac{\mu_i}{\sigma_i^2},\quad  \tilde{x}_i=\frac{x_i}{\sigma_i^2}\bigg),
\end{aligned}
```
we can infer the following:
+ natural parameter: ``\displaystyle \theta = \frac{\mu}{\sigma^2}``
+ ``\displaystyle \psi(\theta) = \frac{\sigma^2}{2}\theta^2``
+ weight: ``\displaystyle w_i = \sigma_i^{-2}``
+ expectation parameter: ``\eta_i = \psi ^\prime (\theta_i) = \sigma^2 \theta = \mu_i``

Suppose we have ``n`` observations ``X_i \sim \mathrm{N}(\mu_i,\sigma_i^2)``.
The functions  ```iso_Normal``` and ```iso_Normal!``` solve the corresponding problem. Their input and output are as follows.

Input 
+ ``x = (x_1,\dots,x_n)``
+ ``\mathrm{variance} = (\sigma_1^2,\dots,\sigma_n^2)``
Output
+ ``\hat{\eta} = (\hat{\eta}_1,\dots,\hat{\eta}_n) = (\hat{\mu}_1,\dots,\hat{\mu}_n)``, which is monotone

```@docs
iso_Normal
```

```@docs
iso_Normal!
```

---
### Binomial

The probability mass function (PMF) of the **Binomial Distribution** with the number of trials ``N_i`` and the probability of success ``r_i`` in an individual trial is given by


```math
p_i(x_i\mid r_i) = \begin{pmatrix} N_i \\ x_i \end{pmatrix} r_i^{x_i} (1-r_i)^{N_i-x_i}.
```

Since
```math
\begin{aligned}
\log p_i(x_i\mid r_i) &= \big( \log r_i - \log(1-r_i) \big) x_i + N_i \log (1-r_i) + \log \begin{pmatrix} N_i \\ x_i \end{pmatrix} \\
&= \theta_i x_i - N_i \log(1+\mathrm{e}^\theta_i) + \log \begin{pmatrix} N_i \\ x_i \end{pmatrix},
\end{aligned}
```
we can infer the following:
+ natural parameter: ``\displaystyle \theta_i = \log r_i - \log(1-r_i) \quad \bigg( r_i = \frac{\mathrm{e}^\theta_i}{1+\mathrm{e}^\theta_i} \bigg)``
+ ``\displaystyle \psi(\theta_i) =  \log(1+\mathrm{e}^\theta_i)``
+ weight: ``\displaystyle w_i = N_i``
+ expectation parameter: ``\displaystyle \eta_i  = \psi ^\prime (\theta) =  \frac{\mathrm{e}^\theta_i}{1+\mathrm{e}^\theta_i} = r_i``

Suppose we have ``n`` observations ``X_i \sim \mathrm{Bi}(N_i,r_i)``.
The functions  ```iso_Binomial``` and ```iso_Binomial!``` solve the corresponding problem. Their input and output are as follows.

Input 
+ ``x = (x_1,\dots,x_n)``
+ ``w = (N_1,\dots,N_2)``
Output
+ ``\hat{\eta} = (\hat{\eta}_1,\dots,\hat{\eta}_n) = (\hat{r}_1,\dots,\hat{r}_n)``, which is monotone

```@docs
iso_Binomial
```

```@docs
iso_Binomial!
```

---
### Poisson

The PMF of the **Poisson Distribution** with the average rate of occurrence ``\lambda`` is gicen by

```math
p(x_i\mid \lambda_i) = \frac{\lambda_i^{x_i} \mathrm{e}^{-\lambda_i}}{x_i!}.
```

Since
```math
\begin{aligned}
\log p(x_i\mid \lambda_i) &= x_i \log \lambda_i-\lambda_i-\log x_i! \\
&= \theta_i x_i-\mathrm{e}^\theta_i-\log x_i!,
\end{aligned}
```
we can infer the following:
+ natural parameter: ``\displaystyle \theta_i = \log\lambda_i``
+ ``\displaystyle \psi(\theta_i) = \mathrm{e}^\theta_i``
+ weight: ``\displaystyle w_i = 1``
+ expectation parameter: ``\displaystyle \eta_i  = \psi ^\prime (\theta_i) = \mathrm{e}^\theta_i = \lambda_i``

Suppose we have ``n`` observations ``X_i \sim \mathrm{Po}(\lambda_i)``.
The functions  ```iso_Poisson``` and ```iso_Poisson!``` solve the corresponding problem. Their input and output are as follows.

Input 
+ ``x = (x_1,\dots,x_n)``
Output
+ ``\hat{\eta} = (\hat{\eta}_1,\dots,\hat{\eta}_n) = (\hat{\lambda}_1,\dots,\hat{\lambda}_n)``, which is monotone

```@docs
iso_Poisson
```

```@docs
iso_Poisson!
```
---
### Chisq

The PDF of the **Chi Squared Distribution** (typically written ``\chi^2``) with ``d_i`` degrees of freedom and scale parameter ``s_i`` is given by

```math
p_i(x_i\mid s_i) = \frac{1}{\Gamma(d_i/2) (2s_i)^{d_i/2}} x^{d_i/2-1} \mathrm{e}^{-\frac{x_i}{2s_i}} .
```

Since 
```math
\begin{aligned}
\log p_i(x_i\mid s_i) &= -\frac{1}{2s_i} x_i - \frac{d_i}{2} \log (2s_i) + \bigg( \frac{d_i}{2}-1\bigg)\log x_i - \log \Gamma (d_i/2)  \\
&= \theta_i x_i - \frac{d_i}{2} (-\log (-\theta_i)) + \bigg( \frac{d_i}{2}-1\bigg)\log x_i - \log \Gamma (d_i/2) ,
\end{aligned}
```
we can infer the following:
+ natural parameter: ``\displaystyle \theta_i = -\frac{1}{2s_i}``
+ ``\displaystyle \psi(\theta_i) = - \log (\theta_i)``
+ weight: ``\displaystyle w_i = \frac{d_i}{2}``
+ expectation parameter: ``\displaystyle \eta_i  = \psi ^\prime (\theta_i) = -\frac{1}{\theta_i} = 2s_i``

Suppose we have ``n`` observations ``X_i \sim s_i \chi^2(d_i)``.
The functions  ```iso_Chisq``` and ```iso_Chisq!``` solve the corresponding problem. Their input and output are as follows.

Input 
+ ``x = (x_1,\dots,x_n)``
+ ``w = (d_1,\dots,d_2)``
Output
+ ``\hat{\eta} = (\hat{\eta}_1,\dots,\hat{\eta}_n) = (2\hat{s}_1,\dots,2\hat{s}_n)``, which is monotone

```@docs
iso_Chisq
```

```@docs
iso_Chisq!
```

