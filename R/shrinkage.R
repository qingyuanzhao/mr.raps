#' Fit a Gaussian mixture deconvolution model
#'
#' @param z A vector of z-scores.
#' @param n Number of mixture components.
#' @param ntry Number of random initializations.
#' @param force.mu.zero Should the means be forced to zero?
#' @param diagnosis Logical indicator for showing diagnostic plots.
#'
#' @return A list of \code{p} (mixture proportion), \code{mu} (mean), \code{sigma} (standard deviation).
#'
#' @details This function assumes that \eqn{z} is distributed as \eqn{N(\gamma, 1)} and \eqn{\gamma} follows a Gaussian mixture model. It fits this deconvolution model by maximum likelihood and outputs the estimated mixture distribution.
#'
#'
fit.mixture.model <- function(z, n = 2, ntry = 10, force.mu.zero = TRUE, diagnosis = FALSE) {

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

    if (diagnosis) {

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
