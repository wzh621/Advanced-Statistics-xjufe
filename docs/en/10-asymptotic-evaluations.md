# Asymptotic Evaluations {#asymptotic-evaluations}

This chapter corresponds to *chap-10.tex* and its four LaTeX fragments. It covers consistency, asymptotic variance and efficiency, relative efficiency, bootstrap standard errors, Huber M-estimation, asymptotic LRT distributions, and approximate likelihood intervals. The source has no external figures or R files; two base-R examples illustrate bootstrap variation and robust loss.

## Consistent estimators {#consistent-estimators}

Asymptotic properties concern a sequence $W_n=W_n(X_1,\ldots,X_n)$ rather than one fixed-sample statistic.

::: {.definition}
The sequence $W_n$ is **consistent** for $\theta$ if, for every $\epsilon>0$ and $\theta\in\Theta$,
$$\lim_{n\to\infty}P_\theta(|W_n-\theta|<\epsilon)=1,$$
that is, $W_n\overset p\to\theta$.
:::

::: {.example .source-numbered}
**Example 10.1.2 (Consistency of the sample mean).** If $X_i\overset{\mathrm{iid}}\sim N(\theta,1)$, prove that $\bar X_n$ is consistent for $\theta$.
:::

<details class="course-details solution-details">
<summary><strong>View detailed solution</strong></summary>

Since $\bar X_n\sim N(\theta,1/n)$,
$$P_\theta(|\bar X_n-\theta|<\epsilon)
=P(-\epsilon\sqrt n<Z<\epsilon\sqrt n)\longrightarrow1.$$
Thus $\bar X_n$ is consistent.

</details>

## MLE regularity conditions and consistency {#mle-regularity-consistency}

The source lists a common sufficient condition set:

1. $X_1,\ldots,X_n$ are iid from $f(x\mid\theta)$;
2. the model is identifiable;
3. densities have common support and are differentiable in $\theta$;
4. the true $\theta_0$ is interior to an open part of the parameter space;
5. $f$ is three times differentiable in $\theta$, with differentiation under the integral allowed;
6. locally, $|\partial^3\log f(x\mid\theta)/\partial\theta^3|\le M(x)$ and $E_{\theta_0}M(X)<\infty$.

The first four support the source's consistency result; the last two are added for asymptotic normality and efficiency. Modern versions can use different weaker assumptions, but the source condition set is preserved.

::: {.theorem .source-numbered}
**Theorem 10.1.6 (Consistency of MLEs).** Under the consistency regularity conditions, if $\hat\theta$ is the MLE and $\tau$ is continuous, then
$$P_\theta\{|\tau(\hat\theta)-\tau(\theta)|\ge\epsilon\}\longrightarrow0.$$
Hence $\tau(\hat\theta)$ is consistent.
:::

<details class="course-details proof-details">
<summary><strong>View derivation</strong></summary>

The textbook gives only a proof outline. For fixed $\vartheta$, the strong law gives

$$\frac1n\log L(\vartheta\mid\mathbf X)
=\frac1n\sum_{i=1}^n\log f(X_i\mid\vartheta)
\longrightarrow E_\theta\{\log f(X\mid\vartheta)\}$$

almost surely. The Kullback--Leibler inequality and identifiability make the limit uniquely maximal at the true parameter $\theta$. Under the additional uniform-convergence and compactness conditions cited by the source, the maximizer $\hat\theta$ converges in probability to $\theta$. Continuity then gives
$\tau(\hat\theta)\overset p\to\tau(\theta)$.

The source refers the fully rigorous general proof to its regularity discussion and an external text, so this disclosure intentionally remains an outline.

</details>

## Limiting and asymptotic variance {#limiting-asymptotic-variance}

Consistency asks whether an estimator reaches its target; efficiency concerns its rate and normalized variation.

::: {.definition}
If $k_n\operatorname{Var}(T_n)\to\tau^2<\infty$, then $\tau^2$ is a **limit of variances**. If
$$a_n\{T_n-\tau(\theta)\}\overset d\longrightarrow N(0,\sigma^2),$$
then $\sigma^2$ is the **asymptotic variance** for that normalization.
:::

For example, $n\operatorname{Var}(\bar X_n)=\sigma^2$. Yet $1/\bar X_n$ can have infinite exact variance in a normal model, while the delta method gives
$$E(1/\bar X_n)\approx1/\mu,\qquad
\operatorname{Var}(1/\bar X_n)\approx\frac{\sigma^2}{n\mu^4}.$$
Thus a limit of finite-sample variances and the variance of a limiting distribution are distinct concepts.

## Asymptotic efficiency and MLEs {#asymptotic-efficiency}

::: {.definition}
If
$$\sqrt n\{W_n-\tau(\theta)\}\overset d\longrightarrow N\{0,v(\theta)\},$$
and
$$v(\theta)=\frac{[\tau'(\theta)]^2}{I_1(\theta)},\qquad
I_1(\theta)=E_\theta\left[\left\{\frac\partial{\partial\theta}\log f(X\mid\theta)\right\}^2\right],$$
then $W_n$ is asymptotically efficient.
:::

::: {.theorem .source-numbered}
**Theorem 10.1.12 (Asymptotic efficiency of MLEs).** Under the stated regularity conditions, if $\hat\theta$ is the MLE and $\tau$ has the required continuity and differentiability, then
$$\sqrt n\{\tau(\hat\theta)-\tau(\theta)\}\overset d\longrightarrow N\{0,v(\theta)\},$$
where $v(\theta)$ is the Cramér--Rao bound. Thus $\tau(\hat\theta)$ is consistent and asymptotically efficient.
:::

<details class="course-details proof-details">
<summary><strong>View derivation</strong></summary>

The textbook first treats $\hat\theta$. Write
$\ell_n(\theta)=\sum_{i=1}^n\log f(X_i\mid\theta)$. A Taylor expansion of the score around the true value $\theta_0$ gives

$$0=\ell_n'(\hat\theta)
=\ell_n'(\theta_0)+(\hat\theta-\theta_0)\ell_n''(\theta_0)
+o_p(\sqrt n\,|\hat\theta-\theta_0|).$$

Thus

$$\sqrt n(\hat\theta-\theta_0)
=\frac{n^{-1/2}\ell_n'(\theta_0)}{-n^{-1}\ell_n''(\theta_0)}+o_p(1).$$

The central limit theorem and weak law yield

$$n^{-1/2}\ell_n'(\theta_0)\Rightarrow N\{0,I_1(\theta_0)\},
\qquad -n^{-1}\ell_n''(\theta_0)\overset p\longrightarrow I_1(\theta_0).$$

Slutsky's theorem gives asymptotic variance $1/I_1(\theta_0)$. Applying the delta method to $\tau$ gives
$[\tau'(\theta_0)]^2/I_1(\theta_0)$, the Cramér--Rao bound. The source leaves this final extension to Exercise 10.7, so it is presented here as the indicated delta-method step.

</details>

## Asymptotic relative efficiency {#asymptotic-relative-efficiency}

::: {.definition}
If
$$\sqrt n\{W_n-\tau(\theta)\}\Rightarrow N(0,\sigma_W^2),\qquad
\sqrt n\{V_n-\tau(\theta)\}\Rightarrow N(0,\sigma_V^2),$$
then
$$\operatorname{ARE}(V_n,W_n)=\frac{\sigma_W^2}{\sigma_V^2}.$$
Under this convention, ARE above one means that $V_n$ has smaller asymptotic variance.
:::

## Bootstrap standard errors {#bootstrap-standard-errors}

Resample $n$ observations with replacement from the empirical distribution. Enumerating all $n^n$ resamples gives
$$\operatorname{Var}^*(\hat\theta)=\frac1{n^n-1}\sum_{i=1}^{n^n}(\hat\theta_i^*-\bar{\hat\theta}^*)^2.$$
With $B$ Monte Carlo resamples, use
$$\operatorname{Var}_B^*(\hat\theta)=\frac1{B-1}\sum_{i=1}^B(\hat\theta_i^*-\bar{\hat\theta}^*)^2.$$
One must distinguish Monte Carlo convergence as $B\to\infty$ for fixed data from bootstrap consistency as $n\to\infty$.


``` r
set.seed(1001)
x <- c(4.1, 5.3, 4.8, 6.2, 5.7, 4.9, 5.5, 6.0)
B <- 4000
boot_mean <- replicate(B, mean(sample(x, replace = TRUE)))
hist(boot_mean, breaks = 32, probability = TRUE,
     col = "#DCEAF4", border = "white",
     xlab = "Bootstrap sample mean", main = "")
abline(v = mean(x), col = "#C43C39", lwd = 2)
legend("topright", sprintf("Bootstrap SE = %.3f", sd(boot_mean)), bty = "n")
```

<div class="figure" style="text-align: center">
<img src="10-asymptotic-evaluations_files/figure-html/en-chap10-bootstrap-se-1.png" alt="Bootstrap distribution and standard error of a sample mean" width="90%" />
<p class="caption">(\#fig:en-chap10-bootstrap-se)Bootstrap distribution and standard error of a sample mean</p>
</div>

## Huber loss and M-estimation {#huber-m-estimators}

Huber loss is quadratic near zero and linear in the tails:
$$\rho_k(u)=\begin{cases}u^2/2,&|u|\le k,\\ k|u|-k^2/2,&|u|>k.\end{cases}$$
An M-estimator minimizes $\sum_i\rho_k(x_i-a)$. With $\psi=\rho'$, it solves
$\sum_i\psi(x_i-\hat\theta)=0$. If $E_{\theta_0}\psi(X-\theta_0)=0$, then
$$-\frac1{\sqrt n}\sum_{i=1}^n\psi(X_i-\theta_0)
\Rightarrow N\left(0,E_{\theta_0}[\psi(X-\theta_0)^2]\right),$$
and local linearization of the estimating equation yields asymptotic normality of $\hat\theta$.


``` r
u <- seq(-4, 4, length.out = 401)
k <- 1.5
huber <- ifelse(abs(u) <= k, u^2 / 2, k * abs(u) - k^2 / 2)
plot(u, u^2 / 2, type = "l", lwd = 2, col = "#1F77B4",
     xlab = "Residual u", ylab = "Loss")
lines(u, huber, lwd = 2, col = "#C43C39")
legend("topleft", c("Squared loss", "Huber loss"),
       col = c("#1F77B4", "#C43C39"), lty = 1, lwd = 2, bty = "n")
```

<div class="figure" style="text-align: center">
<img src="10-asymptotic-evaluations_files/figure-html/en-chap10-huber-loss-1.png" alt="Squared loss and Huber loss" width="90%" />
<p class="caption">(\#fig:en-chap10-huber-loss)Squared loss and Huber loss</p>
</div>

## Asymptotic distribution of LRTs {#asymptotic-lrt-distribution}

::: {.theorem .source-numbered}
**Theorem 10.3.1 (Asymptotic LRT distribution for a simple null).** For $H_0:\theta=\theta_0$ against $H_1:\theta\ne\theta_0$, if the iid model and MLE satisfy the regularity conditions, then under $H_0$,
$$-2\log\lambda(\mathbf X)\overset d\longrightarrow\chi_1^2.$$
:::

<details class="course-details proof-details">
<summary><strong>Expand proof</strong></summary>

Let $\ell_n(\theta)=\log L(\theta\mid\mathbf X)$. Expanding around $\hat\theta$ gives

$$\ell_n(\theta_0)=\ell_n(\hat\theta)
+\ell_n'(\hat\theta)(\theta_0-\hat\theta)
+\frac12\ell_n''(\tilde\theta)(\theta_0-\hat\theta)^2,$$

where $\tilde\theta$ lies between $\theta_0$ and $\hat\theta$. Since an interior MLE has $\ell_n'(\hat\theta)=0$,

$$-2\log\lambda(\mathbf X)
=\left\{\frac{-\ell_n''(\tilde\theta)}n\right\}
\{\sqrt n(\hat\theta-\theta_0)\}^2.$$

Under $H_0$, the first factor converges in probability to $I_1(\theta_0)$, while Theorem 10.1.12 gives
$\sqrt n(\hat\theta-\theta_0)\Rightarrow N\{0,1/I_1(\theta_0)\}$. Slutsky's theorem therefore makes the product converge to the square of a standard normal variable, namely $\chi_1^2$.

</details>

More generally, if the unrestricted model has dimension $p$ and the null has $q$ free dimensions, then under smooth interior regularity,
$$-2\log\lambda(\mathbf X)\Rightarrow\chi^2_{p-q}.$$
An approximate level $\alpha$ test rejects for $-2\log\lambda\ge\chi^2_{p-q,1-\alpha}$. Boundary parameters, nonidentifiability, and parameter-dependent support can invalidate this limit.

## Approximate maximum-likelihood intervals {#approximate-ml-intervals}

Observed information gives
$$\widehat{\operatorname{Var}}\{h(\hat\theta)\}
\approx\frac{[h'(\hat\theta)]^2}
{-\left.\dfrac{\partial^2}{\partial\theta^2}\log L(\theta\mid\mathbf x)\right|_{\theta=\hat\theta}}.$$
MLE asymptotic normality and Slutsky's theorem yield the Wald interval
$$h(\hat\theta)\pm z_{\alpha/2}
\sqrt{\widehat{\operatorname{Var}}\{h(\hat\theta)\}}.$$

## LRT intervals {#lrt-intervals}

Inverting the LRT gives the approximate $1-\alpha$ confidence set
$$\left\{\theta:-2\log\frac{L(\theta\mid\mathbf x)}{L(\hat\theta\mid\mathbf x)}
\le\chi^2_{1,1-\alpha}\right\}.$$

::: {.example .source-numbered}
**Example 10.4.3 (Binomial LRT interval).** If $Y=\sum_iX_i$ for independent Bernoulli$(p)$ observations, invert the asymptotic LRT to construct a confidence set for $p$.
:::

<details class="course-details solution-details">
<summary><strong>View detailed solution</strong></summary>

The unrestricted MLE is $\hat p=y/n$, so the confidence set is
$$\left\{p:-2\log\left[
\frac{p^y(1-p)^{n-y}}{\hat p^y(1-\hat p)^{n-y}}
\right]\le\chi^2_{1,1-\alpha}\right\}.$$
Its endpoints are normally obtained by one-dimensional root finding.

For $0<y<n$, strict concavity of the log likelihood usually makes this a single interval containing $\hat p$; solve the equality separately on either side. Boundary samples $y=0$ or $y=n$ require the corresponding one-sided boundary treatment.

</details>

## Chapter summary {#asymptotic-evaluations-summary}

Consistency, asymptotic normality, and efficiency describe large-sample accuracy and variation; bootstrap uses the empirical distribution to approximate sampling variation; M-estimation limits the influence of extreme residuals; and Wilks-type chi-square approximations connect LRTs, tests, and intervals, subject to regularity conditions.
