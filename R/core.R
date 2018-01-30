#' Main function
#'
#' \code{mr.raps} is the main function.
#'
#' @param b_exp A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param b_out A vector of SNP effects on the outcome variable, usually obtained from a GWAS.
#' @param se_exp A vector of standard errors of \code{b_exp}.
#' @param se_out A vector of standard errors of \code{b_out}.
#' @param over.dispersion Should the model consider overdispersion (systematic pleiotropy)? Default is FALSE.
#' @param loss.function Either the squared error loss (\code{l2}) or robust loss functions/scores (\code{huber} or \code{tukey}).
#' @param diagnosis Should the function returns diagnostic plots and results? Default is FALSE
#' @param se.method How should the standard error be estimated? Either by sandwich variance formula (default and recommended) or the bootstrap.
#' @param k Threshold parameter in the Huber and Tukey loss functions.
#' @param B Number of bootstrap resamples
#' @param suppress.warning Should warning messages be suppressed?
#' @param initialization Method to initialize the robust estimator. "Mode" is not supported currently.
#' @param niter Maximum number of interations to solve the estimating equations.
#' @param tol Numerical precision.
#'
#' @return A list
#' \describe{
#' \item{beta.hat}{Estimated causal effect}
#' \item{beta.se}{Standard error of \code{beta.hat}}
#' \item{beta.p.value}{Two-sided p-value of \code{beta.hat}}
#' \item{tau2.hat}{Overdispersion parameter if \code{over.dispersion = TRUE}}
#' \item{tau2.se}{Standard error of \code{tau2.hat}}
#' \item{std.resid}{Standardized residuals of each SNP, returned if \code{diagnosis = TRUE}}
#' \item{beta.hat.loo}{Leave-one-out estimates of \code{beta.hat}, returned if \code{diagnosis = TRUE}}
#' \item{beta.hat.bootstrap}{Median of the bootstrap estimates, returned if \code{se.method = "bootstrap"}}
#' \item{beta.se.bootstrap}{Median absolute deviation of the bootstrap estimates, returned if \code{se.method = "bootstrap"}}
#' }
#'
#' @references Qingyuan Zhao, Jingshu Wang, Jack Bowden, Dylan S. Small. Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. \url{https://arxiv.org/abs/1801.09652}.
#'
#' @import stats
#' @export
#'
#' @examples
#'
#' data(bmi.sbp)
#' attach(bmi.sbp)
#'
#' ## All estimators
#' mr.raps.all(beta.exposure, beta.outcome, se.exposure, se.outcome)
#'
#' ## Diagnostic plots
#' res <- mr.raps(beta.exposure, beta.outcome, se.exposure, se.outcome,
#' diagnosis = TRUE)
#' res <- mr.raps(beta.exposure, beta.outcome, se.exposure, se.outcome,
#' TRUE, diagnosis = TRUE)
#' res <- mr.raps(beta.exposure, beta.outcome, se.exposure, se.outcome,
#' TRUE, "tukey", diagnosis = TRUE)
#'
#' detach(bmi.sbp)
#'
#' data(bmi.bmi)
#' attach(bmi.bmi)
#'
#' ## Because both the exposure and the outcome are BMI, the true "causal" effect should be 1.
#'
#' ## All estimators
#' mr.raps.all(beta.exposure, beta.outcome, se.exposure, se.outcome)
#'
#' detach(bmi.bmi)
#'
mr.raps <- function(b_exp, b_out, se_exp, se_out,
                    over.dispersion = FALSE,
                    loss.function = c("l2", "huber", "tukey"),
                    diagnosis = FALSE,
                    se.method = c("sandwich", "bootstrap"),
                    k = switch(loss.function[1], l2 = NULL, huber = 1.345, tukey = 4.685),
                    B = 1000,
                    suppress.warning = FALSE) {

    loss.function <- match.arg(loss.function, c("l2", "huber", "tukey"))
    se.method <- match.arg(se.method, c("sandwich", "bootstrap"))

    if (loss.function == "l2") {
        if (!over.dispersion) {
            fit <- mr.raps.simple(b_exp, b_out, se_exp, se_out, diagnosis = diagnosis)
        } else {
            fit <- mr.raps.overdispersed(b_exp, b_out, se_exp, se_out, diagnosis = diagnosis, suppress.warning = suppress.warning)
        }
    } else {
        if (!over.dispersion) {
            fit <- mr.raps.simple.robust(b_exp, b_out, se_exp, se_out, loss.function, k, diagnosis = diagnosis)
        } else {
            fit <- mr.raps.overdispersed.robust(b_exp, b_out, se_exp, se_out, loss.function, k, suppress.warning = suppress.warning, diagnosis = diagnosis)
        }
    }

    if (se.method == "bootstrap") {
        fit.bootstrap <- list()
        for (b in 1:B) {
            if (b %% round(B/10) == 0) {
                print(paste0("Bootstrap: ", round(b / B * 100), "% finished."))
            }
            s <- sample(1:length(b_exp), replace = TRUE)
            fit.bootstrap[[b]] <- tryCatch(unlist(mr.raps(b_exp[s], b_out[s], se_exp[s], se_out[s], over.dispersion, loss.function, se.method = "sandwich", k = k, suppress.warning = TRUE)), error = function(e) {print(e); NA})
        }
        fit.bootstrap <- data.frame(do.call(rbind, fit.bootstrap))
        fit <- c(fit,
                 list(beta.hat.bootstrap = median(fit.bootstrap$beta.hat),
                      beta.se.bootstrap = mad(fit.bootstrap$beta.hat)))
    }

    fit
}

#' \code{mr.raps.all}: Quick analysis with all six methods
#'
#' @describeIn mr.raps
#'
#' @export
#'
mr.raps.all <- function(b_exp, b_out, se_exp, se_out) {
    res <- data.frame()
    for (over.dispersion in c(FALSE, TRUE)) {
        for (loss.function in c("l2", "huber", "tukey")) {
            out <- mr.raps(b_exp, b_out, se_exp, se_out, over.dispersion, loss.function)
            out <- data.frame(out)
            out$over.dispersion <- over.dispersion
            out$loss.function <- loss.function
            res <- rbind(res, out[c("over.dispersion", "loss.function", "beta.hat", "beta.se")])
        }
    }
    res
}

#' \code{mr.raps.simple}: No overdispersion, l2 loss
#'
#' @describeIn mr.raps
#'
#' @export
#' @importFrom nortest ad.test
#' @import stats graphics
#'
mr.raps.simple <- function(b_exp, b_out, se_exp, se_out, diagnosis = FALSE) {

    profile.loglike <- function(beta) {
        - (1/2) * sum((b_out - b_exp * beta)^2 / (se_out^2 + se_exp^2 * beta^2)) # - (1/2) * sum(log(se_out^2 + beta^2 * se_exp^2))
    }

    bound <- quantile(abs(b_out / b_exp), 0.95) * 2

    beta.hat <- optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
    while (abs(beta.hat) > 0.95 * bound) {
        bound <- bound * 2
        beta.hat <- optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
    }
    score.var <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2 + se_exp^2 * se_out^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)
    I <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

    dif <- b_out - beta.hat * b_exp
    dif.var <- se_out^2 + beta.hat^2 * se_exp^2 # + (score.var / I^2 + 2 * beta.hat^2 * se_exp^2) * (b_exp^2 - se_exp^2)
    chi.sq.test <- sum((dif / sqrt(dif.var))^2)

    if (diagnosis) {
        std.resid <- (b_out - b_exp * beta.hat) / sqrt((se_out^2 + beta.hat^2 * se_exp^2))
        par(mfrow = c(1, 2))
        qqnorm(std.resid)
        abline(0, 1)
        ## Leave-one-out
        beta.hat.loo <- rep(NA, length(b_out))
        beta.hat.loo <- rep(NA, length(b_out))
        if (length(b_out) > 100) {
            a <- quantile(abs(b_exp / se_exp), 1 - 100/length(b_out))
        } else {
            a <- 0
        }
        for (i in 1:length(b_out)) {
            if (abs(b_exp[i] / se_exp[i]) > a) {
                beta.hat.loo[i] <- mr.raps.simple(b_exp[-i], b_out[-i], se_exp[-i], se_out[-i])$beta.hat
            }
        }
        plot(abs(b_exp / se_exp), beta.hat.loo)
        abline(h = beta.hat)
        ## print(ks.test(abs(rnorm(100000)), abs(b_out - b_exp * beta.hat) / sqrt((se_out^2 + beta.hat^2 * se_exp^2))))
        ## library(goftest)
        ## print(ad.test(abs(b_out - b_exp * beta.hat) / sqrt((se_out^2 + beta.hat^2 * se_exp^2)), null = function(x) 2 * pnorm(abs(x)) - 1))
        print(ad.test(std.resid))
        print(shapiro.test(std.resid))
    }

    out <- list(beta.hat = beta.hat,
                beta.se = sqrt(score.var/I^2),
                beta.p.value = min(1, 2 * (1 - pnorm(abs(beta.hat) / sqrt(score.var/I^2)))),
                naive.se = sqrt(1/I),
                chi.sq.test = chi.sq.test)
    if (diagnosis) {
        out$std.resid <- std.resid
        out$beta.hat.loo <- beta.hat.loo
    }
    out
}

#' \code{mr.raps.overdispersed}: Overdispersion, l2 loss
#'
#' @describeIn mr.raps
#'
#' @import stats graphics
#' @importFrom nortest ad.test
#' @export
#'
mr.raps.overdispersed <- function(b_exp, b_out, se_exp, se_out, initialization = c("simple", "mode"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5) {

    initialization <- match.arg(initialization, c("simple", "mode"))

    profile.loglike.fixbeta <- function(beta, tau2) {
        alpha.hat <- 0
        - (1/2) * sum(se_exp^2 * (log(tau2 + se_out^2 + se_exp^2 * beta^2))) - (1/2) * sum(se_exp^2 * (b_out - alpha.hat - b_exp * beta)^2 / (tau2 + se_out^2 + se_exp^2 * beta^2))
    }

    profile.loglike.fixtau <- function(beta, tau2) {
        alpha.hat <- 0
        - (1/2) * sum((b_out - alpha.hat - b_exp * beta)^2 / (tau2 + se_out^2 + se_exp^2 * beta^2))
    }

    bound.beta <- quantile(abs(b_out / b_exp), 0.95) * 10
    bound.tau2 <- quantile(se_out^2, 0.95) * 10

    ## Initialization
    if (initialization == "mode") {
        stop("Initialization by mode estimator is currently not supported.")
        ## beta.hat <- mode.estimator(b_exp, b_out, se_exp, se_out)
        ## tau2.hat <- 0
    } else {
        fit <- mr.raps.simple(b_exp, b_out, se_exp, se_out)
        beta.hat <- fit$beta.hat
        tau2.hat <- 0
    }

    for (iter in 1:niter) {
        ## print(c(beta.hat, tau2.hat))
        beta.hat.old <- beta.hat
        tau2.hat.old <- tau2.hat

        ## Estimate tau2
        tau2.hat <- optimize(function(tau2) profile.loglike.fixbeta(beta.hat, tau2), bound.tau2 * c(0, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
        ## while (abs(tau2.hat) > 0.95 * bound.tau2 && int.extend <= niter) {
        ##     int.extend <- int.extend + 1
        ##     bound.tau2 <- bound.tau2 * 2
        ##     tau2.hat <- optimize(function(tau2) profile.loglike.fixbeta(beta.hat, tau2), bound.tau2 * c(0, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
        ## }
        ## if (int.extend == niter) {
        ##     stop("Failed to find tau.")
        ## }
        if (tau2.hat > bound.tau2 * 0.95) {
            warning("Estimated overdispersion seems abnormaly large.")
        }

        ## Estimate beta
        beta.hat <- optimize(function(beta) profile.loglike.fixtau(beta, tau2.hat), bound.beta * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
        int.extend <- 0
        while (abs(beta.hat) > 0.95 * bound.beta && int.extend <= niter) {
            int.extend <- int.extend + 1
            bound.beta <- bound.beta * 2
            beta.hat <- optimize(function(beta) profile.loglike.fixtau(beta, tau2.hat), bound.beta * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
        }
        if (int.extend == niter) {
            stop("Failed to find beta.")
        }
        if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) + abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
            break
        }
    }

    if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) + abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) > tol && (!suppress.warning)) {
        warning("Did not converge when solving the estimating equations. Consider to increase niter or decrease tol.")
    }

    if ((tau2.hat <= min(se_out^2) / 5) && (!suppress.warning)) {
        warning("The estimated overdispersion parameter is quite small. Consider to use the simple model without overdispersion.")
    }

    score.var <- diag(
        c(sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) + (b_out^2 - tau2.hat - se_out^2) * se_exp^2 + se_exp^2 * (tau2.hat + se_out^2))/(tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
          sum(2 * se_exp^4 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2))
    )
    I <- matrix(c(
        - sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) + (b_out^2 - tau2.hat - se_out^2) * se_exp^2) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
        0,
        - sum(se_exp^2 * beta.hat / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
        - sum(se_exp^2/ (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)), 2, 2)

    asymp.var <- solve(I) %*% score.var %*% t(solve(I))

    if (diagnosis) {
        std.resid <- (b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2))
        par(mfrow = c(1, 2))
        qqnorm(std.resid)
        abline(0, 1)
        ## Leave-one-out
        beta.hat.loo <- rep(NA, length(b_out))
        if (length(b_out) > 100) {
            a <- quantile(abs(b_exp / se_exp), 1 - 100/length(b_out))
        } else {
            a <- 0
        }
        for (i in 1:length(b_out)) {
            if (abs(b_exp[i] / se_exp[i]) > a) {
                beta.hat.loo[i] <- mr.raps.overdispersed(b_exp[-i], b_out[-i], se_exp[-i], se_out[-i], suppress.warning = TRUE)$beta.hat
            }
        }
        plot(abs(b_exp / se_exp), beta.hat.loo)
        abline(h = beta.hat)
        ## print(ks.test(abs(rnorm(100000)), abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2))))
        ## library(goftest)
        ## print(ad.test(abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2)), null = function(x) 2 * pnorm(abs(x)) - 1))
        print(ad.test(std.resid))
        print(shapiro.test(std.resid))
    }

    out <- list(beta.hat = beta.hat,
                tau2.hat = tau2.hat,
                beta.se = sqrt(asymp.var[1, 1]),
                tau2.se = sqrt(asymp.var[2, 2]),
                beta.p.value = min(1, 2 * (1 - pnorm(abs(beta.hat) / sqrt(asymp.var[1, 1])))))

    if (diagnosis) {
        out$std.resid <- std.resid
        out$beta.hat.loo <- beta.hat.loo
    }

    out

}

#' Huber loss function and its derivatives
#'
#' @import stats
#' @keywords internal
#'
rho.huber <- function(r, k = 1.345, deriv = 0) {
    if (deriv == 0) {
        return(ifelse(abs(r) <= k,
                      r^2/2,
                      k * (abs(r) - k/2)))
    } else if (deriv == 1) {
        return(ifelse(abs(r) <= k,
                      r,
                      k * sign(r)))
    } else if (deriv == 2) {
        return(ifelse(abs(r) <= k,
                      1,
                      0))
    } else {
        stop("deriv must be 0, 1, or 2.")
    }
}

#' Tukey's beweight loss function and its derivatives
#'
#' @import stats
#' @keywords internal
#'
rho.tukey <- function(r, k = 4.685, deriv = 0) {
    if (deriv == 0) {
        pmin(1 - (1 - (r/k)^2)^3, 1)
    } else if (deriv == 1) {
        r * (1 - (r / k)^2)^2 * (abs(r) <= k)
    } else if (deriv == 2) {
        t <- (r/k)^2
        ifelse(t < 1, (1 - t) * (1 - 5 * t), 0)
    }
}

#' \code{mr.raps.simple.robust}: No overdispersion, robust loss
#'
#' @describeIn mr.raps
#'
#' @import stats graphics
#' @importFrom nortest ad.test
#' @export
#'
mr.raps.simple.robust <- function(b_exp, b_out, se_exp, se_out, loss.function = c("huber", "tukey"), k = switch(loss.function[1], huber = 1.345, tukey = 4.685), diagnosis = FALSE) {

    loss.function <- match.arg(loss.function, c("huber", "tukey"))
    rho <- switch(loss.function,
                  huber = function(r, ...) rho.huber(r, k, ...),
                  tukey = function(r, ...) rho.tukey(r, k, ...))

    delta <- integrate(function(x) x * rho(x, deriv = 1) * dnorm(x), -Inf, Inf)$value
    c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value
    c2 <- integrate(function(x) x^2 * rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c3 <- integrate(function(x) x^2 * rho(x, deriv = 2) * dnorm(x), -Inf, Inf)$value

    robust.loglike <- function(beta) {
        - sum(rho((b_out - b_exp * beta) / sqrt(se_out^2 + se_exp^2 * beta^2)))
    }

    bound <- quantile(abs(b_out / b_exp), 0.95) * 2

    beta.hat <- optimize(robust.loglike, bound * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
    while (abs(beta.hat) > 0.95 * bound) {
        bound <- bound * 2
        beta.hat <- optimize(robust.loglike, bound * c(-1, 1), maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
    }

    score.var <- c1 * sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2 + se_exp^2 * se_out^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)
    I <- delta * sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

    dif <- b_out - beta.hat * b_exp
    dif.var <- se_out^2 + beta.hat^2 * se_exp^2 # + (score.var / I^2 + 2 * beta.hat^2 * se_exp^2) * (b_exp^2 - se_exp^2)
    chi.sq.test <- sum((dif / sqrt(dif.var))^2)

    if (diagnosis) {
        std.resid <- (b_out - b_exp * beta.hat) / sqrt((se_out^2 + beta.hat^2 * se_exp^2))
        par(mfrow = c(1, 2))
        qqnorm(std.resid)
        abline(0, 1)
        ## Leave-one-out
        beta.hat.loo <- rep(NA, length(b_out))
        if (length(b_out) > 100) {
            a <- quantile(abs(b_exp / se_exp), 1 - 100/length(b_out))
        } else {
            a <- 0
        }
        for (i in 1:length(b_out)) {
            if (abs(b_exp[i] / se_exp[i]) > a) {
                beta.hat.loo[i] <- mr.raps.simple.robust(b_exp[-i], b_out[-i], se_exp[-i], se_out[-i], loss.function = loss.function, k = k)$beta.hat
            }
        }
        plot(abs(b_exp / se_exp), beta.hat.loo)
        abline(h = beta.hat)
        ## print(ks.test(abs(rnorm(100000)), abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2))))
        ## library(goftest)
        ## print(ad.test(abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2)), null = function(x) 2 * pnorm(abs(x)) - 1))
        print(ad.test(std.resid))
        print(shapiro.test(std.resid))
    }

    asymp.var <- solve(I) %*% score.var %*% t(solve(I))

    out <- list(beta.hat = beta.hat,
                beta.se = sqrt(score.var/I^2) ,
                naive.se = sqrt(1/I),
                chi.sq.test = chi.sq.test,
                beta.p.value = min(1, 2 * (1 - pnorm(abs(beta.hat) / sqrt(asymp.var[1, 1])))))

    if (diagnosis) {
        out$std.resid <- std.resid
        out$beta.hat.loo <- beta.hat.loo
    }

    out

}

#' \code{mr.raps.overdispersed.robust}: Overdispersed, robust loss
#'
#' @describeIn mr.raps
#'
#' @import stats graphics
#' @importFrom nortest ad.test
#' @export
#'
mr.raps.overdispersed.robust <- function(b_exp, b_out, se_exp, se_out, loss.function = c("huber", "tukey"), k = switch(loss.function[1], huber = 1.345, tukey = 4.685), initialization = c("l2", "mode"), suppress.warning = FALSE, diagnosis = FALSE, niter = 20, tol = .Machine$double.eps^0.5) {

    loss.function <- match.arg(loss.function, c("huber", "tukey"))
    initialization <- match.arg(initialization, c("l2", "mode"))
    rho <- switch(loss.function,
                  huber = function(r, ...) rho.huber(r, k, ...),
                  tukey = function(r, ...) rho.tukey(r, k, ...))

    delta <- integrate(function(x) x * rho(x, deriv = 1) * dnorm(x), -Inf, Inf)$value
    c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value
    c2 <- integrate(function(x) x^2 * rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c3 <- integrate(function(x) x^2 * rho(x, deriv = 2) * dnorm(x), -Inf, Inf)$value

    robust.loglike.fixtau <- function(beta, tau2) {
        alpha.hat <- 0
        - (1/2) * sum(rho((b_out - alpha.hat - b_exp * beta) / sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)))
    }

    robust.E <- function(beta, tau2) {
        t <- (b_out - beta * b_exp) / sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)
        se_exp^2 * (t * rho(t, deriv = 1) - delta) / (tau2 + se_out^2 + se_exp^2 * beta^2)
    }

    ## Initialize
    if (initialization == "mode") {
        stop("Initialization by mode estimator is currently not supported.")
        ## beta.hat <- mode.estimator(b_exp, b_out, se_exp, se_out)
        ## tau2.hat <- 0
    } else {
        fit <- mr.raps.overdispersed(b_exp, b_out, se_exp, se_out, suppress.warning = TRUE)
        beta.hat <- fit$beta.hat
        tau2.hat <- fit$tau2.hat
    }
    bound.beta <- quantile(abs(b_out / b_exp), 0.95) * 10
    bound.tau2 <- quantile(se_out^2, 0.95) * 10

    for (iter in 1:niter) {
        beta.hat.old <- beta.hat
        tau2.hat.old <- tau2.hat
        tau2.hat <- tryCatch(uniroot(function(tau2) sum(robust.E(beta.hat, tau2)), bound.tau2 * c(0, 1), extendInt = "yes", tol = bound.tau2 * .Machine$double.eps^0.25)$root, error = function(e) {warning("Did not find a solution for tau2."); 0})
        if (tau2.hat < 0) {
            tau2.hat <- 0
        }
        if (tau2.hat > bound.tau2 * 0.95) {
            warning("Estimated overdispersion seems abnormaly large.")
        }
        beta.hat <- optim(beta.hat, function(beta) robust.loglike.fixtau(beta, tau2.hat), method = "L-BFGS-B", lower = -bound.beta, upper = bound.beta, control = list(fnscale = -1))$par
        int.extend <- 0
        while (abs(beta.hat) > 0.95 * bound.beta && int.extend <= niter) {
            int.extend <- int.extend + 1
            bound.beta <- bound.beta * 2
            beta.hat <- optimize(function(beta) robust.loglike.fixtau(beta, tau2.hat), bound.beta * c(-1, 1), maximum = TRUE)$maximum
        }
        if (int.extend == niter) {
            stop("Failed to find beta.")
        }
        if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) + abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
            break
        }
    }

    if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) + abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) > tol) {
        warning("Did not converge when solving the estimating equations. Consider to increase niter or decrease tol.")
    }

    if ((tau2.hat <= min(se_out^2) / 5) && (!suppress.warning)) {
        warning("The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion.")
    }

    score.var <- diag(
        c(c1 * sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) + (b_out^2 - tau2.hat - se_out^2) * se_exp^2 + se_exp^2 * (tau2.hat + se_out^2))/(tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
          (c2 / 2) * sum(2 * se_exp^4 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2))
    )
    I <- matrix(c(
        - delta * sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) + (b_out^2 - tau2.hat - se_out^2) * se_exp^2) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
        0,
        - delta * sum(se_exp^2 * beta.hat / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
        - (delta + c3) / 2 * sum(se_exp^2/ (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)), 2, 2)

    asymp.var <- solve(I) %*% score.var %*% t(solve(I))

    if (diagnosis) {
        std.resid <- (b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2))
        par(mfrow = c(1, 2))
        qqnorm(std.resid)
        abline(0, 1)
        ## Leave-one-out
        beta.hat.loo <- rep(NA, length(b_out))
        if (length(b_out) > 100) {
            a <- quantile(abs(b_exp / se_exp), 1 - 100/length(b_out))
        } else {
            a <- 0
        }
        for (i in 1:length(b_out)) {
            if (abs(b_exp[i] / se_exp[i]) > a) {
                beta.hat.loo[i] <- mr.raps.overdispersed.robust(b_exp[-i], b_out[-i], se_exp[-i], se_out[-i], loss.function = loss.function, k = k, suppress.warning = TRUE)$beta.hat
            }
        }
        plot(abs(b_exp / se_exp), beta.hat.loo)
        abline(h = beta.hat)
        ## print(ks.test(abs(rnorm(100000)), abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2))))
        ## library(goftest)
        ## print(ad.test(abs(b_out - b_exp * beta.hat) / sqrt((tau2.hat + se_out^2 + beta.hat^2 * se_exp^2)), null = function(x) 2 * pnorm(abs(x)) - 1))
        print(ad.test(std.resid))
        print(shapiro.test(std.resid))
    }

    out <- list(beta.hat = beta.hat,
                tau2.hat = tau2.hat,
                beta.se = sqrt(asymp.var[1, 1]),# / sqrt(efficiency),
                tau2.se = sqrt(asymp.var[2, 2]),
                beta.p.value = min(1, 2 * (1 - pnorm(abs(beta.hat) / sqrt(asymp.var[1, 1])))))# / sqrt(efficiency),
    if (diagnosis) {
        out$std.resid <- std.resid
        out$beta.hat.loo <- beta.hat.loo
    }

    out

}

## #' Modified weights IVW
## #'
## #' This function implements the modified 2nd order weighted procedure.
## #'
## #' @inheritParams mr.raps
## #'
## #' @references Bowden, Jack, M. Fabiola Del Greco, Cosetta Minelli, Debbie Lawlor, Nuala Sheehan, John Thompson, and George Davey Smith. "Improving the accuracy of two-sample summary data Mendelian randomization: moving beyond the NOME assumption." bioRxiv (2017): 159442.
## #'
## #' @return A list
## #' \describe{
## #' \item{beta.hat}{Estimated causal effect}
## #' \item{beta.se}{Standard error of \code{beta.hat}}
## #' }
## #'
## #' @export
## #'
## ivw.modified <- function(b_exp, b_out, se_exp, se_out) {
##     BIV = b_out/b_exp
##     W1 = 1/(se_out^2/b_exp^2)
##     BIVw1 = BIV*sqrt(W1)
##     sW1 = sqrt(W1)

##     IVWfitR1 = summary(lm(BIVw1 ~ -1+sW1))
##     Bhat1 = IVWfitR1$coef[1]
##     DF = length(BIV)-1

##     W3 = 1/(se_out^2/b_exp^2 + (Bhat1^2)*se_exp^2/b_exp^2)
##     BIVw3 = BIV*sqrt(W3)
##     sW3 = sqrt(W3)
##     IVWfitR3 = summary(lm(BIVw3 ~ -1+sW3))
##     Bhat3 = IVWfitR3$coef[1]
##     phi_IVW3 = IVWfitR3$sigma^2
##     QIVW3 = DF*phi_IVW3 # Q statistic
##     Qp3 = 1-pchisq(QIVW3,DF) # p-value
##     Q3ind = W3*(BIV - Bhat3)^2 # individual Q contribution vector

##     return(list(beta.hat = Bhat3, beta.se = IVWfitR3$coef[2]))
## }
