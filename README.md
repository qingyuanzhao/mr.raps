# MR-RAPS: Robust statistical inference for Mendelian randomization with many weak instruments

*mr.raps* is an R package for two-sample summary-data Mendelian randomization using the robust adjusted profile score (MR-RAPS). To install the package, run

```
library(devtools)
install_github("qingyuanzhao/mr.raps")
```

The main function is *mr.raps*. You can find examples by running

```
library(mr.raps)
example(mr.raps) ## Recommended procedure
example(mr.raps.shrinkage) ## General function for empirical partially Bayes estimator
```

In May 2018, a new general function *mr.raps.shrinkage* is added. You can choose whether weight shrinkage should be used by the option *shrinkage*. You can still use the original MR-RAPS procedure using the function *mr.raps.mle*, which usually gives similar results as setting *shrinkage = FALSE* in *mr.raps.shrinkage*.

References:
* Zhao, Qingyuan, Jingshu Wang, Jack Bowden, and Dylan S. Small. "Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score." arXiv:1801.09652 (2018).
* Zhao, Qingyuan, Yang Chen, Jingshu Wang, and Dylan S. Small. "Powerful genome-wide design and robust statistical inference in two-sample summary-data Mendelian randomization." arXiv:1804.07371 (2018).