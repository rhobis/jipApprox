jipApprox
======================================================

[![Travis-CI Build Status](https://travis-ci.org/rhobis/jipApprox.svg?branch=master)](https://travis-ci.org/rhobis/jipApprox)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/UPSvarApprox)](https://cran.r-project.org/package=jipApprox)
[![](https://cranlogs.r-pkg.org/badges/grand-total/jipApprox)](https://cran.r-project.org/package=jipApprox)

Description 
-----------------
This package to approximate joint-inclusion probabilities in Unequal Probability
Sampling, or to find Monte Carlo approximations of first and second-order inclusion
probabilities of a general sampling design.


The main functions are:

- `jip_approx()`: returns a matrix of approximated joint-inclusion probabilities for 
unequal probability sampling design with high entropy;
- `jip_MonteCarlo()`: produces a matrix of first and second order inclusion probabilities
for a given sampling design, approximated through Monte Carlo simulation. 
This method of approximation is more flexible but also computer-intensive.
- `HTvar()`: returns the Horvitz-Thompson or Sen-Yates-Grundy variance or their estimates,
computed using true inclusion probabilities or an approximation obtained by 
`jip_approx()` or `jip_MonteCarlo()`.


Installation
------------

Currently, the package can be installed only from GitHub:

``` r
# if not present, install 'devtools' package
install.packages("devtools")
devtools::install_github("rhobis/jipApprox")
```

Usage
-----

``` r
library(UPSvarApprox)

### Generate population data ---
N <- 20; n <- 5

set.seed(0)
x <- rgamma(500, scale=10, shape=5)
y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )

pik <- n * x/sum(x)

### Approximate joint-inclusion probabilities for high entropy designs ---
pikl <- jip_approx(pik, method='Hajek')
pikl <- jip_approx(pik, method='HartleyRao')
pikl <- jip_approx(pik, method='Tille')
pikl <- jip_approx(pik, method='Brewer1')
pikl <- jip_approx(pik, method='Brewer2')
pikl <- jip_approx(pik, method='Brewer3')
pikl <- jip_approx(pik, method='Brewer4')

### Approximate inclusion probabilities through Monte Carlo simulation ---
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "brewer")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "tille")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "poisson")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "maxEntropy")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "randomSystematic")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "systematic")
pikl <- jip_MonteCarlo(x=pik, n = n, replications = 100, design = "sampford")

```



More
----

- Please, report any bug or issue [here](https://github.com/rhobis/jipApprox/issues).
- For more information, please contact the manteiner at `roberto.sichera@unipa.it`. 

