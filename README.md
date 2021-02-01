# MR.RAPS: Robust statistical inference for Mendelian randomization with many weak instruments

## Setup
*mr.raps* is an R package for two-sample summary-data Mendelian randomization using the robust adjusted profile score (MR-RAPS). To install the most up-to-date version, run the following command in R

```
library(devtools)
install_github("qingyuanzhao/mr.raps")
```

The [CRAN version of this package](https://cran.r-project.org/web/packages/mr.raps/index.html) is currently not being maintained.

## Examples
The main function is *mr.raps*. You can find examples by running

```
library(mr.raps)
example(mr.raps) ## Recommended procedure
example(mr.raps.shrinkage) ## General function for empirical partially Bayes estimator
```

A more in-depth real-data example can be found [here](http://www.statslab.cam.ac.uk/~qz280/talks/mr_challenge_report_marked.pdf).

## Updates
In May 2018, a new general function *mr.raps.shrinkage* is added. You can choose whether weight shrinkage should be used by the option *shrinkage*. You can still use the original MR-RAPS procedure using the function *mr.raps.mle*, which usually gives similar results as setting *shrinkage = FALSE* in *mr.raps.shrinkage*.

We are updating the [multivariate branch](https://github.com/qingyuanzhao/mr.raps/tree/multivariate) to accomdate multivariable MR and sample overlap.

## Getting help
More information (including references and tutorials on Mendelian randomization) can be found on [this webpage](http://www.statslab.cam.ac.uk/~qz280/post/mr-software/). 

Please report issues and suggestions on the software using the [GitHub issues tracker](https://github.com/qingyuanzhao/mr.raps/issues).
