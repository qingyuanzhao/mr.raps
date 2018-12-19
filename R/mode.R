## #' Mode estimator
## #'
## #' @keywords internal
## #'
## #' @details This function implements the a mode estimator. The idea belongs to Fernando Pires Hartwig and the author heard it from Jack Bowden. This estimator is used as a robust initialization of the RAPS estimator.
## #'
## #'
## mode.estimator <- function(b_exp, b_out, se_exp, se_out) {

##     fieller.ci <- function(a, b, v11, v22, v12, alpha = 0.05) {
##         g <- qnorm(alpha)^2 * v22 / b^2
##         cen <- a/b - g * v12 / v22
##         suppressWarnings(len <- abs(abs(qnorm(alpha)) / b * sqrt(v11 - 2 * a / b * v12 + a^2/b^2 * v22 - g * (v11 - v12^2/v22))))
##         c(cen - len, cen + len) / (1 - g)
##     }

##     ci <- do.call(rbind, lapply(1:length(b_exp), function(i) fieller.ci(b_out[i], b_exp[i], se_out[i]^2, se_exp[i]^2, 0)))
##     ci <- ci[complete.cases(ci), ]

##     candidate <- c(ci)
##     no.cover <- sapply(candidate, function(x) sum(x >= ci[, 1] & x <= ci[, 2]))
##     mean(candidate[which(no.cover == max(no.cover))])

## }
