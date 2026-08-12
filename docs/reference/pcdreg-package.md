# pcdreg: Semiparametric Regression for Panel Count Data

Panel count data arise when a recurrent event process is observed only
at intermittent examination times, so that the number of events between
consecutive examinations is known but the event times themselves are
not.
[`panelrate()`](https://www.sundayu.me/pcdreg/reference/panelrate.md)
fits the semiparametric proportional rate model \$\$E\[dN(t) \mid X(t)\]
= \exp(\beta' X(t)) \\ d\Lambda(t)\$\$ to such data, allowing the
covariates \\X(t)\\ to vary over time.

Estimation treats the baseline cumulative rate \\\Lambda\\
nonparametrically and maximises the likelihood that a nonhomogeneous
Poisson process would imply, using an EM algorithm that augments the
observed counts with latent per-examination-time Poisson counts. The
Poisson assumption is a working device only: the estimator stays
consistent and asymptotically normal when it fails, and the default
covariance estimator is robust to that failure.

## References

Sun, D., Guo, Y., Li, Y., Tu, W. and Sun, J. (2024). A robust approach
for regression analysis of panel count data. *Bernoulli* **30**(4),
3251–3275. [doi:10.3150/23-BEJ1713](https://doi.org/10.3150/23-BEJ1713)

## See also

Useful links:

- <https://github.com/dayusun/pcdreg>

- <https://www.sundayu.me/pcdreg/>

- Report bugs at <https://github.com/dayusun/pcdreg/issues>

## Author

**Maintainer**: Dayu Sun <dayusun@iu.edu>
