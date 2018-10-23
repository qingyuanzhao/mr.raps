#' Fit a Gaussian mixture deconvolution model
#'
#' @param z A vector of z-scores.
#' @param n Number of mixture components.
#' @param ntry Number of random initializations.
#' @param force.mu.zero Should the means be forced to zero?
#' @param diagnostics Logical indicator for showing diagnostic plots.
#'
#' @return A list of \code{p} (mixture proportion), \code{mu} (mean), \code{sigma} (standard deviation).
#'
#' @details This function assumes that \eqn{z} is distributed as \eqn{N(\gamma, 1)} and \eqn{\gamma} follows a Gaussian mixture model. It fits this deconvolution model by maximum likelihood and outputs the estimated mixture distribution.
#'
#' @examples
#' z <- c(sqrt(2) * rnorm(900), sqrt(17) * rnorm(100))
#' ## So the correct sigma = (1, 4) and p = (0.9, 0.1)
#' fit.mixture.model(z)
#'
#' @export
#'
fit.mixture.model <- function(z, n = 2, ntry = 10, force.mu.zero = TRUE, diagnostics = FALSE) {

    loglike <- function(param, z) {
        n <- length(param) / 3
        p <- param[1:n]
        p <- p / sum(p)
        mu <- param[n + 1:n]
        sigma <- param[2*n + 1:n]

        l <- matrix(0, length(z), length(p))
        for (i in 1:length(p)) {
            l[, i] <- p[i] * dnorm(z, mu[i], sqrt(sigma[i]^2 + 1))
        }
        - sum(log(rowSums(l)))
    }

    get.random.init <- function(n = 2) {
        p <- rgamma(n, 1)
        p <- p / sum(p)
        mu <- rnorm(n)
        sigma <- c(0.1, rexp(n - 1) * 2)
        c(p, mu, sigma)
    }

    if (force.mu.zero) {
        mu.lower <- -1e-8
        mu.upper <- 1e-8
    } else {
        mu.lower <- -Inf
        mu.upper <- Inf
    }
    res <- list()
    for (try in 1:ntry) {
        try(res[[try]] <- optim(get.random.init(n),
                                function(param) loglike(param, z),
                                method = "L-BFGS-B",
                                lower = c(rep(0.01, n), rep(mu.lower, n), rep(0.01, n)),
                                upper = c(rep(0.99, n), rep(mu.upper, n), rep(Inf, n))))
    }
    i <- which.min(sapply(1:length(res), function(i) {tmp <- res[[i]]$value; if (is.null(tmp)) {Inf} else {tmp}}))
    res <- res[[i]]

    param <- res$par
    p <- param[1:n]
    p <- p / sum(p)
    mu <- param[n + 1:n]
    mu[abs(mu) < 1e-6] <- 0
    sigma <- param[2*n + 1:n]

    if (diagnostics) {

        print("Estimated mixture model:")
        print(paste0("p = ", signif(p, 3), ", ",
                     "mu = ", signif(mu, 3), ", ",
                     "sigma = ", signif(sigma, 3)))

        print("Generating diagnostic plots...")

        t <- seq(-10, 10, 0.1)
        fitted.cdf <- sapply(1:n, function(i) p[i] * pnorm(t, mu[i], sqrt(sigma[i]^2 + 1)))
        fitted.cdf <- rowSums(fitted.cdf)
        emp.cdf <- sapply(t, function(t) mean(z < t))

        par(mfrow = c(1, 2))
        plot(t, fitted.cdf, col = "red")
        lines(t, emp.cdf)
        plot(fitted.cdf, emp.cdf)
        abline(0,1)
    }

    list(p = p, mu = mu, sigma = sigma)
}

#' Compute the posterior mean under spike-and-slab Gaussian prior
#'
#' @param z a vector of z-scores
#' @param sigma a vector of standard deviations of \code{z} (if \code{sigma} is a single number, it is expanded to a vector)
#' @param p prior: mixture proportion
#' @param mu prior: mean
#' @param sigma.prior prior: standard deviation
#' @param deriv compute the posterior mean (\code{deriv = 0}) or its derivative (\code{deriv = 1})
#'
#' @details Similar to \code{fit.mixture.model}, this function assumes that \eqn{z} is distributed as \eqn{N(\gamma, 1)} and \eqn{\gamma} follows a Gaussian mixture model. The function computes the posterior mean \eqn{E[\gamma|z]}.
#'
#' @return a vector
#'
#' @examples
#'
#' require(mr.raps)
#' data(lipid.cad)
#' data <- subset(lipid.cad, lipid == "hdl" & restrict &
#' gwas.selection == "teslovich_2010" &
#' gwas.outcome == "cardiogramplusc4d")
#' z <- data$beta.exposure / data$se.exposure
#' prior.param <- fit.mixture.model(z)
#'
#' z.seq <- seq(-5, 5, 0.1)
#' gamma.hat <- posterior.mean(z.seq, 1, prior.param$p, prior.param$mu, prior.param$sigma)
#' gamma.hat.deriv <- posterior.mean(z.seq, 1, prior.param$p,
#' prior.param$mu, prior.param$sigma, deriv = 1)
#' par(mfrow = c(1, 2))
#' plot(z.seq, gamma.hat, type = "l")
#' plot(z.seq, gamma.hat.deriv, type = "l")
#'
#' @export
#'
posterior.mean <- function(z, sigma, p, mu, sigma.prior, deriv = 0) {
    if (length(z) > 1) {
        if (length(sigma) == 1) {
            sigma <- rep(sigma, length(z))
        }
        stopifnot(length(z) == length(sigma))
        return(sapply(1:length(z), function(k)
            posterior.mean(z[k], sigma[k], p, mu, sigma.prior, deriv)))
    }
    if (length(p) == 1){
        p <- c(p, 1 - p)
    }
    stopifnot(length(p) == length(mu) && length(p) == length(sigma.prior))
    stopifnot(deriv %in% c(0, 1))

    mu.tilde <- (z/sigma^2 + mu/sigma.prior^2) / (1/sigma^2 + 1/sigma.prior^2)
    ## sigma.tilde <- sqrt(1/(1/sigma^2 + 1/sigma.prior^2)) ## not needed
    p.tilde <- p * dnorm(z, mu, sqrt(sigma^2 + sigma.prior^2))

    if (deriv == 0) {
        return(sum(p.tilde * mu.tilde) / sum(p.tilde))
    } else {
        diff.mu.tilde <- (1/sigma^2) / (1/sigma^2 + 1/sigma.prior^2)
        diff.p.tilde <- - p * dnorm(z, mu, sqrt(sigma^2 + sigma.prior^2)) * (z - mu) / (sigma^2 + sigma.prior^2)
        return(sum(diff.p.tilde * mu.tilde + p.tilde * diff.mu.tilde) / sum(p.tilde) - sum(p.tilde * mu.tilde) * sum(diff.p.tilde) / sum(p.tilde)^2)
    }
}

#' Main function for RAPS (shrinkage weights)
#'
#' @inheritParams mr.raps.mle
#' @param shrinkage If shrinkage (empirical partially Bayes) should be used. Shrinkage does not affect the unbiasedness of the estimating equations and generally will increase the estimation accuracy. If TRUE, \code{prior.param} must be provided.
#' @param prior.param Parameters of the Gaussian spike-and-slab prior
#' @param num.init Number of initializations.
#' @param multiple.root.warning How to handle multiple roots of the estimating equations? When this happens, the results of \code{mr.raps.shrinkage} are less reliable. This parameter can take three values: 0---nothing will be done; 1---a warning is given; 2---an error is given. Default is 1.
#'
#' @details
#' \code{mr.raps.shrinkage} is the main function for RAPS in conjunction with empirical partially Bayes. It is more general than the first generation \code{mr.raps.mle} function and should be preferred in practice. With the option \code{shrinkage = TRUE}, it essentially reduces to \code{mr.raps.mle}. In that case, the main difference is that the standard errors in \code{mr.raps.shrinkage} are computed based on observed information (and also an empirical estimate of the variance of the score function). This is preferred over using the plugged-in Fisher information in \code{mr.raps.mle}. See Efron and Hinkley (1978) referenced below.
#'
#' Because the estimating equations are highly non-linear, it is possible that there are multiple roots. To overcome this issue, we use multiple initializations (controlled by \code{num.init}) around the \code{mr.raps.mle} point estimate. A warning is given if there seems to be another finite root, and no solution is returned if there are two roots close to the initialization.
#'
#' @references
#' Qingyuan Zhao, Q., Chen, Y., Wang, J., and Small, D. S. (2018). A genome-wide design and an empirical partially Bayes approach to increase the power of Mendelian randomization, with application to the effect of blood lipids on cardiovascular disease. <arXiv:1804.07371>.
#' Efron, B. and Hinkley, D. V. (1978). Assessing the accuracy of the maximum likelihood estimator: Observed versus expected Fisher information. Biometrika, 65(3), 457--483.
#'
#' @examples
#'
#' require(mr.raps)
#' data(lipid.cad)
#' data <- subset(lipid.cad, lipid == "hdl" & restrict &
#' gwas.selection == "teslovich_2010" & gwas.outcome == "cardiogramplusc4d")
#' z <- data$beta.exposure / data$se.exposure
#' prior.param <- fit.mixture.model(z)
#'
#' ## Results
#' mr.raps.shrinkage(data$beta.exposure, data$beta.outcome,
#' data$se.exposure, data$se.outcome, TRUE, "huber", shrinkage = FALSE)
#' \donttest{
#' mr.raps.shrinkage(data$beta.exposure, data$beta.outcome,
#' data$se.exposure, data$se.outcome, TRUE, "huber", shrinkage = TRUE,
#' prior.param = prior.param)
#' }
#'
#' @export
#'
#' @importFrom rootSolve multiroot
#'
mr.raps.shrinkage <- function(b_exp, b_out, se_exp, se_out, over.dispersion = FALSE, loss.function = c("l2", "huber", "tukey"), k = switch(loss.function[1], l2 = 2, huber = 1.345, tukey = 4.685), shrinkage = FALSE, prior.param = NULL, diagnostics = FALSE, se.method = c("sandwich", "bootstrap"), num.init = 10, multiple.root.warning = 1) {

    se.method <- match.arg(se.method)
    if (se.method == "bootstrap") {
        B <- 100
        beta.hat <- rep(0, B)
        tau2.hat <- rep(0, B)
        for (b in 1:B) {
            if (b %% 20 == 0) {
                print(paste0("Bootstrap iteration: ", b, "; Total: ", B, "."))
            }
            s <- sample.int(length(b_exp), replace = TRUE)
            res <- mr.raps.shrinkage(b_exp[s], b_out[s], se_exp[s], se_out[s], over.dispersion, loss.function, k, shrinkage, prior.param, diagnostics = FALSE, num.init = 10)
            beta.hat[b] <- res$beta.hat
            tau2.hat[b] <- res$tau2.hat
        }
        hist(beta.hat)
        return(list(beta.hat = mean(beta.hat, na.rm = TRUE),
                    beta.se = sd(beta.hat, na.rm = TRUE),
                    tau2.hat = mean(tau2.hat, na.rm = TRUE),
                    tau2.se = sd(tau2.hat, na.rm = TRUE)))
    }

    loss.function <- match.arg(loss.function, c("l2", "huber", "tukey"))
    rho <- switch(loss.function,
                  l2 = function(r, ...) rho.l2(r, k, ...),
                  huber = function(r, ...) rho.huber(r, k, ...),
                  tukey = function(r, ...) rho.tukey(r, k, ...))

    delta <- integrate(function(x) x * rho(x, deriv = 1) * dnorm(x), -Inf, Inf)$value
    c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value
    c2 <- integrate(function(x) x^2 * rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c3 <- integrate(function(x) x^2 * rho(x, deriv = 2) * dnorm(x), -Inf, Inf)$value

    get.t <- function(beta, tau2) {
        (b_out - b_exp * beta) / sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)
    }

    get.gamma.hat <- function(beta, tau2, deriv = c("0", "beta", "tau2")) {
        deriv <- match.arg(deriv, c("0", "beta", "tau2"))
        gamma.mle <- (b_exp/ se_exp^2 + beta * b_out/ (se_out^2 + tau2)) / (1 / se_exp^2 + beta^2 / (se_out^2 + tau2))
        if (deriv == "0") {
            if (shrinkage) {
                a <- posterior.mean(gamma.mle / se_exp,
                                    sqrt(1 / (1 + beta^2 * se_exp^2 / (se_out^2 + tau2))),
                                    prior.param$p, prior.param$mu, prior.param$sigma)
                gamma.hat <- a * se_exp
            } else {
                gamma.hat <- gamma.mle
            }
            return(gamma.hat)
        }

        if (deriv == "beta") {
            gamma.mle.deriv <- (b_out/ (se_out^2 + tau2)) / (1 / se_exp^2 + beta^2 / (se_out^2 + tau2))
        } else { ## deriv == "tau2"
            gamma.mle.deriv <- - beta * b_out / (se_out^2 + tau2)^2 / (1 / se_exp^2 + beta^2 / (se_out^2 + tau2))
        }
        if (shrinkage) {
            shrinkage.deriv <- posterior.mean(
                gamma.mle / se_exp,
                sqrt(1 / (1 + beta^2 * se_exp^2 / (se_out^2 + tau2))),
                prior.param$p, prior.param$mu, prior.param$sigma,
                deriv = 1)
            return(gamma.mle.deriv * shrinkage.deriv)
        } else {
            return(gamma.mle.deriv)
        }
    }

    psi <- function(param) {
        if (!over.dispersion) {
            beta <- param[1]
            tau2 <- 0
        } else {
            beta <- param[1]
            tau2 <- param[2]
        }
        t <- get.t(beta, tau2)
        gamma.hat <- get.gamma.hat(beta, tau2)
        v <- beta^2 * se_exp^2 + se_out^2 + tau2
        psi1 <- sum(gamma.hat * rho(t, deriv = 1) / sqrt(v))
        if (over.dispersion) {
            psi2 <- sum((t * rho(t, deriv = 1) - delta) / v)
            return(c(psi1, psi2))
        } else {
            return(psi1)
        }
    }

    if (!over.dispersion) {
        res <- mr.raps.simple(b_exp, b_out, se_exp, se_out)
        init.param <- res$beta.hat
    } else {
        res <- mr.raps.mle(b_exp, b_out, se_exp, se_out, over.dispersion, loss.function, suppress.warning = TRUE)
        init.param <- c(res$beta.hat, res$tau2.hat)
    }

    beta.init <- init.param[1] + c(0, 2 * init.param[1] * rnorm(num.init - 1))
    res <- list()
    beta <- rep(NA, num.init)
    for (i in 1:num.init) {
        init.param[1] <- beta.init[i]
        suppressWarnings(res[[i]] <- multiroot(psi, init.param))
        if ((over.dispersion) && (!is.na(res[[i]]$root[2])) && (res[[i]]$root[2] > median(se_out) * 10)) {
            beta[i] <- NA
        } else {
            if (is.na(res[[i]]$estim.precis)) {
                beta[i] <- NA
            } else {
                beta[i] <- res[[i]]$root[1]
            }
        }
    }
    j <- which.min(abs(beta - beta.init[1]))
    if (length(j) == 0) {
        warning("Cannot find solution with finite over.dispersion. Using tau2 = 0.")
        return(mr.raps.shrinkage(b_exp, b_out, se_exp, se_out, FALSE, loss.function, k, shrinkage, prior.param, diagnostics))
        ## res <- multiroot(function(beta) psi(c(beta, 0))[1], init.param[1])
        ## estimated.param <- c(res$root, 0)
    }

    if (multiple.root.warning > 0) {
        for(i in 1:num.init) {
            if (!is.na(beta[i]) && abs(beta[i] - beta[j]) > 1e-4 && abs(beta[i] - beta.init[1]) < 100 * abs(beta[j] - beta.init[1])) {
                warning(paste("The estimating equations might have another finite root. The closest root is beta =", beta[j], "and the other root is beta =", beta[i], "and the initialization is beta =", beta.init[1]))
            }
            if (multiple.root.warning > 1) {
                if (!is.na(beta[i]) && abs(beta[i] - beta[j]) > 1e-4 && abs(beta[i] - beta.init[1]) < 5 * abs(beta[j] - beta.init[1])) {
                    stop(paste("Found two very close solutions: beta =", beta[j], "and", beta[i]))
                    ## return(list(beta.hat = NA, beta.se = NA, tau2.hat = NA, tau2.se = NA))
                }
            }
        }
    }
    res <- res[[j]]

    if (over.dispersion && res$root[2] < 0) {
        warning("Estimated overdispersion is negative. Using tau2 = 0.")
        return(mr.raps.shrinkage(b_exp, b_out, se_exp, se_out, FALSE, loss.function, k, shrinkage, prior.param, diagnostics))
        ## res <- multiroot(function(beta) psi(c(beta, 0))[1], init.param[1])
        ## estimated.param <- c(res$root, 0)
    } else {
        estimated.param <- res$root
    }

    if (!over.dispersion) {
        beta <- estimated.param[1]
        tau2 <- 0
        t <- get.t(beta, tau2)
        gamma.hat <- get.gamma.hat(beta, tau2)
        v <- beta^2 * se_exp^2 + se_out^2 + tau2
        V1 <- c1 * sum(gamma.hat^2 / v)
        V2 <- sum((get.gamma.hat(beta, tau2, deriv = "beta") * rho(t, deriv = 1) + gamma.hat * delta * (- b_exp) / sqrt(v)) / sqrt(v))
        beta.se <- sqrt(V1 / V2^2)
        tau2.se <- 0
    } else {
        beta <- estimated.param[1]
        tau2 <- estimated.param[2]
        t <- get.t(beta, tau2)
        gamma.hat <- get.gamma.hat(beta, tau2)
        v <- beta^2 * se_exp^2 + se_out^2 + tau2

        V1 <- diag(c(c1 * sum(gamma.hat^2 / v), c2 * sum(1 / v^2)))
        V2 <- matrix(c(sum((get.gamma.hat(beta, tau2, deriv = "beta") * rho(t, deriv = 1) + gamma.hat * delta * (- b_exp) / sqrt(v)) / sqrt(v)),
                       0,
                       sum(get.gamma.hat(beta, tau2, "tau2") * rho(t, deriv = 1) / sqrt(v)),
                       (delta + c3) / 2 * sum(1/v^2)), 2, 2)
        V <- solve(V2) %*% V1 %*% t(solve(V2))
        beta.se <- sqrt(V[1, 1])
        tau2.se <- sqrt(V[2, 2])
    }

    t <- get.t(beta, tau2)
    ## par(mfrow = c(1, 3))
    ## qqnorm(t)
    ## abline(0, 1)
    gamma.hat <- get.gamma.hat(beta, tau2)
    v <- 1 / (1 / se_exp^2 + beta^2 / (se_out^2 + tau2))
    gamma.hat.z <- gamma.hat / sqrt(v)
    t <- t * sign(gamma.hat.z)
    gamma.hat.z <- abs(gamma.hat.z)
    ## plot(gamma.hat.z, t, xlab = "Weight", ylab = "Standardized residual")
    ## plot(rank(- gamma.hat.z), t, xlab = "Rank of weight", ylab = "Standardized residual")
    ## print("Test of independence:")
    ## print(summary(lm(t ~ poly(gamma.hat.z, 2) - 1)))
    out <- list(beta.hat = beta, tau2.hat = tau2, beta.se = beta.se, tau2.se = tau2.se, t = t, gamma.hat.z = gamma.hat.z)
    class(out) <- "mr.raps"

    if (diagnostics) {
        plot(out)
    }
    out

}

#' @describeIn mr.raps.shrinkage Print
#' @param x a \code{mr.raps} object
#' @param ... further arguments (not supported)
#'
#' @export
print.mr.raps <- function(x, ...) {
    print(x[1:4])
}

#' @describeIn mr.raps.shrinkage Diagnostic plots
#' @inheritParams print.mr.raps
#'
#' @export
#'
#' @import ggplot2 gridExtra
#'
plot.mr.raps <- function(x, ...) {

    qhnorm <- function(p) {
        - qnorm(p / 2)
    }

    df <- data.frame(t = x$t, w = x$gamma.hat.z)
    grid.arrange(
        ggplot(df) + aes(x = w, y = t, shape) + geom_point() + geom_smooth(method = "loess") + xlab("Absolute weight") + ylab("Standardized residual"),
        ggplot(df) + aes(x = qhnorm(ppoints(length(t)))[order(order(- abs(t)))], y = abs(t)) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("Theoretical") + ylab("Sample"),
        ncol = 2)
}

#' Recommended \code{mr.raps} procedure
#'
#' @param data A data frame (see Details)
#' @param diagnostics Logical indicator for showing diagnostic plots.
#' @param ... Additional parameters to be passed to \code{mr.raps.shrinkage} (default is \code{shrinkage=FALSE}).
#'
#' @details
#' This function calls \code{mr.raps.shrinkage} with \code{overdispersion = TRUE}, \code{loss.function = "huber"}, \code{shrinkage = TRUE}. The input data frame should contain the following variables:
#' \enumerate{
#' \item beta.exposure
#' \item beta.outcome
#' \item se.exposure
#' \item se.outcome
#' }
#'
#' @import splines
#' @export
#'
#' @examples
#'
#' mr.raps(bmi.sbp)
#'
#' \donttest{
#' require(mr.raps)
#' data(lipid.cad)
#' data <- subset(lipid.cad, lipid == "hdl" & restrict &
#' gwas.selection == "teslovich_2010" & gwas.outcome == "cardiogramplusc4d")
#' mr.raps(data)
#' }
#'
mr.raps <- function(data, diagnostics = TRUE, ...) {
    prior.param <- fit.mixture.model(data$beta.exposure / data$se.exposure)
    out <- mr.raps.shrinkage(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, TRUE, "huber", prior.param = prior.param, ...)

    if (diagnostics) {
        cat(paste0("Estimated causal effect: ", signif(out$beta.hat, 3), ", standard error: ", signif(out$beta.se, 3), ", p-value: ", signif(pnorm(-abs(out$beta.hat / out$beta.se)) * 2, 3), ".\n"))
        cat(paste0("Estimated pleiotropy variance: ", signif(out$tau2.hat, 3), ", standard error: ", signif(out$tau2.se, 3), ", p-value: ", signif(pnorm(-abs(out$tau2.hat / out$tau2.se)) * 2, 3), ".\n"))

        cat(paste0("ANOVA test: are the weights and residuals independent? \n"))
        weights <- out$gamma.hat.z
        std.resids <- out$t
        df <- max(round(length(weights) / 50), 3)
        lm.test <- lm(std.resids ~ bs(weights, df) - 1)
        print(anova(lm.test))

        cat("Showing diagnostic plot ...\n")
        plot(out)
    }

    out
}

#' Generate diagnostic plots for publishing
#'
#' @keywords internal
#'
#' @import ggplot2 splines gridExtra
#'
mr.raps.publish <- function(data) {

    ## require(splines)
    ## require(ggplot2)

    prior.param <- fit.mixture.model(data$beta.exposure / data$se.exposure)
    out1 <- mr.raps.shrinkage(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, TRUE, "huber", shrinkage = FALSE)
    out2 <- mr.raps.shrinkage(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, TRUE, "huber", shrinkage = TRUE, prior.param = prior.param)

    weights <- out1$gamma.hat.z
    std.resids <- out1$t
    df <- max(round(length(weights) / 50), 3)
    lm.test <- lm(std.resids ~ bs(weights, df) - 1)
    p1 <- anova(lm.test)[[5]][1]

    weights <- out2$gamma.hat.z
    std.resids <- out2$t
    df <- max(round(length(weights) / 50), 3)
    lm.test <- lm(std.resids ~ bs(weights, df) - 1)
    p2 <- anova(lm.test)[[5]][1]

    df <- data.frame(SNP = rep(data$SNP, 2),
                     t = c(out1$t, out2$t),
                     w = c(out1$gamma.hat.z, out2$gamma.hat.z),
                     pval.selection = rep(data$pval.selection, 2),
                     weight.method = rep(c("MLE", "Shrinkage"), each = nrow(data)))
    df.label <- data.frame(p = c(p1, p2),
                           beta.hat = c(out1$beta.hat, out2$beta.hat),
                           weight.method = rep(c("MLE", "Shrinkage")))
    out.plot <- ggplot(df) + aes(x = w, y = t) + geom_point(aes(shape = (pval.selection < 5e-8), color = (pval.selection < 5e-8), size = (pval.selection < 5e-8)), alpha = 0.7)
    out.plot <- out.plot + geom_text(x = max(df$w) * 0.5, y = max(df$t) * 1, aes(label = paste("Estimated effect:", as.character(signif(beta.hat, 2)))), data = df.label, size = 2)
    out.plot <- out.plot + geom_text(x = max(df$w) * 0.5, y = max(df$t) * 0.8, aes(label = paste("Heterogeneity p-value:", as.character(signif(p, 2)))), data = df.label, size = 2)
    out.plot <- out.plot + coord_cartesian(ylim = range(df$t) * 1.1) + facet_grid(weight.method ~ .) + geom_smooth(method = "loess", span = 1/3) + xlab("Absolute weight") + ylab("Standardized residual") + scale_shape_discrete(guide = FALSE) + scale_color_discrete(guide = FALSE) + scale_size_discrete(guide = FALSE, range = c(1.5, 2.5)) + theme_bw(base_size = 18)

    qhnorm <- function(p) {
        - qnorm(p / 2)
    }

    qq.plot <- ggplot(df) + aes(x = qhnorm(ppoints(length(t)))[order(order(- abs(t)))], y = abs(t)) + geom_point(aes(shape = (pval.selection < 5e-8), color = (pval.selection < 5e-8), size = (pval.selection < 5e-8))) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + facet_grid(weight.method ~ .) + xlab("Theoretical") + ylab("Sample") + scale_shape_discrete(guide = FALSE) + scale_color_discrete(guide = FALSE) + scale_size_discrete(guide = FALSE, range = c(1.5, 3)) + theme_bw(base_size = 18)

    ## require(gridExtra)
    grid.arrange(out.plot, qq.plot, ncol = 2, widths = c(2, 1))

}
