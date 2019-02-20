#' Robust multivariate MR estimation
#'
#' @param b_exp A matrix of size \code{p * k} for the effect sizes of \code{p} number of independent SNPs on \code{k} risk factors
#' @param b_out A vector of length \code{p} for the effect sizes of the \code{p} SNPs on the outcome
#' @param se_exp A matrix of size \code{p * k} for the standard deviations of the effect sizes in \code{b_exp}
#' @param se_out A vector of length \code{p} for the standard deviations of the effect sizes in \code{b_out}
#' @param cor.mat Either NULL or a \code{k + 1} by \code{k + 1} symmetric matrix. The correlation matrix of estimated effect sizes on the \code{k} risk factors and the outcome. Default is NULL, for the identity matrix
#' @param loss.function Loss function used, one of "l2", "huber" and "tukey"
#' @param k Tuning parameters of the loss function, for loss "l2", it is NA, for loss "huber", default is 1.345 and for loss "tukey", default is 4.685
#' @param suppress.warning Whether suppress warning messages or not, default is FALSE
#' @param diagnosis Run diagnosis analysis based on the residuals or not, default is FALSE
#' @param niter Number of iterations for optimization. Default is 20
#' @param tol Tolerance for convergence, default is the square root of the smallest positive floating number depending on the machine R is running on
#' @import nortest
#' @export
grappleRobustEst <- function(b_exp, b_out,
                             se_exp, se_out,
                             cor.mat = NULL,
                             loss.function = c("l2", "huber", "tukey"),
                             k = switch(loss.function[1],
                                        l2 = NA, huber = 1.345,
                                        tukey = 4.685),
                             suppress.warning = FALSE,
                             diagnosis = FALSE,
                             niter = 20,
                             tol = .Machine$double.eps^0.5) {

    b_exp <- as.matrix(b_exp)
    se_exp <- as.matrix(se_exp)

    loss.function <- match.arg(loss.function, c("l2", "huber", "tukey"))
    if (is.null(cor.mat))
        cor.mat <- diag(rep(1, ncol(b_exp) + 1))
    rho <- switch(loss.function,
                  l2 = function(r, ...) rho.l2(r, ...),
                  huber = function(r, ...) rho.huber(r, k, ...),
                  tukey = function(r, ...) rho.tukey(r, k, ...))

                                        #    delta <- integrate(function(x) x * rho(x, deriv = 1) * dnorm(x), -Inf, Inf)$value
    delta <- integrate(function(x)  rho(x) * dnorm(x), -Inf, Inf)$value

    c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value
                                        #   c2 <- integrate(function(x) x^2 * rho(x, deriv = 1)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c2 <- integrate(function(x) rho(x)^2 * dnorm(x), -Inf, Inf)$value - delta^2
    c4 <- integrate(function(x) rho(x, deriv = 1) * x * dnorm(x), -Inf, Inf)$value
                                        #   c3 <- integrate(function(x) (rho(x, deriv = 2) * x^2 + rho(x, deriv = 1) * x)* dnorm(x), -Inf, Inf)$value
    c3 <- c4

    ## First, calculate t for each SNP
    t_fun <- function(beta, tau2) {
                                        #    print(tau2)
        upper <- b_out - b_exp %*% beta
        temp <- t(t(cbind(se_exp, se_out)) * c(-beta, 1))
        lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)

        return(upper/lower)
    }


    robust.optfun.fixtau <- function(beta, tau2) {
        -sum(rho(t_fun(beta, tau2)))
    }

    robust.E <- function(beta, tau2) {
        return(sum(rho(t_fun(beta, tau2))) - (length(b_out) - 1) * delta)
    }
                                        #   robust.E <- function(beta, tau2) {
                                        #     tt <- t_fun(beta, tau2)
                                        #     return(sum(rho(tt, deriv = 1) * tt) - length(b_out) * delta)
                                        #    }


    ## Initialize
    tau2.hat <- 0
    bound.beta <- apply(abs(b_out / b_exp), 2, function(v)quantile(v[is.finite(v)], probs = 0.95, na.rm = T)) * 2
    bound.tau2 <- quantile(se_out^2, 0.95, na.rm = T) * 2
    if (ncol(b_exp) == 1) {
        beta.seq <- seq(-bound.beta, bound.beta, length.out = 5000)
        beta.hat <- beta.seq[which.max(sapply(beta.seq, robust.optfun.fixtau, tau2 = 0))]
    } else {
        beta.hat <- as.vector(lm(b_out ~ b_exp + 0)$coef)
        temp.fun <- function(bb, idx) {
            beta <- rep(0, length(beta.hat))
            beta[idx] <- bb
            return(robust.optfun.fixtau(beta, tau2=0))
        }

        for (i in 1:length(beta.hat)) {
            beta.seq <- seq(-bound.beta[i], bound.beta[i], length.out = 5000)
            beta.hat[i] <- beta.seq[which.max(sapply(beta.seq, temp.fun, idx = i))]
                                        #  beta.hat[i] <- optimize(function(bb) temp.fun(bb, i),
                                        #                          bound.beta[i] * c(-1, 1), maximum = TRUE,
                                        #                          tol = .Machine$double.eps^0.5)$maximum
        }
    }



    for (iter in 1:niter) {
        print(beta.hat)
        print(tau2.hat)
        beta.hat.old <- beta.hat
        tau2.hat.old <- tau2.hat
        tau2.hat <- tryCatch(uniroot(function(tau2) sum(robust.E(beta.hat, tau2)),
                                     bound.tau2 * c(0, 1), extendInt = "downX",
                                     tol = bound.tau2 * .Machine$double.eps^0.25)$root,
                             error = function(e) {warning("Did not find a solution for tau2."); 0})

        if (tau2.hat < 0) {
            tau2.hat <- 0

        }
                                        #    print(tau2.hat)
        if (tau2.hat > bound.tau2 * 0.95) {
            warning("Estimated overdispersion seems abnormaly large.")
        }
        if (length(beta.hat) == 1) {
            beta.hat <- beta.seq[which.max(sapply(beta.seq, robust.optfun.fixtau, tau2 = tau2.hat))]
        } else
            beta.hat <- optim(beta.hat, function(beta) robust.optfun.fixtau(beta, tau2.hat),
                              method = "L-BFGS-B", lower = -bound.beta, upper = bound.beta,
                              control = list(fnscale = -1))$par
        int.extend <- 0
        while (abs(beta.hat) > 0.95 * bound.beta && int.extend <= niter) {
            int.extend <- int.extend + 1
            bound.beta <- bound.beta * 2
            beta.hat <- optimize(function(beta) robust.optfun.fixtau(beta, tau2.hat),
                                 bound.beta * c(-1, 1), maximum = TRUE)$maximum
        }
        if (int.extend == niter) {
            stop("Failed to find beta.")
        }
        if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) +
            abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
            break
        }
    }

    if (max(abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10)) +
        abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) > tol) {
        warning("Did not converge when solving the estimating equations.
              Consider to increase niter or decrease tol.")
    }


    ## final round of coordinate-wise grid-search optimization to make the non-convex optimization
    ## more likely reach a global optimal point
    for(s in 1:5) {
        tau2.hat <- tryCatch(uniroot(function(tau2) sum(robust.E(beta.hat, tau2)),
                                     bound.tau2 * c(0, 1), extendInt = "downX",
                                     tol = bound.tau2 * .Machine$double.eps^0.25)$root,
                             error = function(e) {warning("Did not find a solution for tau2."); 0})
        if (tau2.hat < 0) {
            tau2.hat <- 0

        }

        temp.fun <- function(bb, idx, beta.hat.temp = beta.hat) {
            beta.hat.temp[idx] <- bb
            return(robust.optfun.fixtau(beta.hat.temp, tau2=tau2.hat))
        }

        for (i in 1:length(beta.hat)) {
            bb <- beta.hat[i] * 50
            beta.seq <- seq(-abs(bb), abs(bb), length.out = 3000)
            beta.hat[i] <- beta.seq[which.max(sapply(beta.seq, temp.fun, idx = i, beta.hat.temp = beta.hat))]
        }


    }

    if ((tau2.hat <= min(se_out^2) / 5) && (!suppress.warning)) {
        warning("The estimated overdispersion parameter is very small.
                Consider using the simple model without overdispersion.")
    }

### calculating asymptotic variance
    ## calculate the derivatives of t_fun at beta.hat, tau2.hat
    r <- ncol(b_exp)
    p <- nrow(b_exp)
    res <- as.vector(b_out - b_exp %*% beta.hat) ## size p
    temp <- t(t(cbind(se_exp, se_out)) * c(-beta.hat, 1))  ## size p * (r + 1)
    var.res <- rowSums((temp %*% cor.mat) * temp) + tau2.hat ## size p
    t.values <- res/sqrt(var.res) ## size p
    ## calculate each row as Sigma_X_j beta - Sigma_{X_jY_j}
    coefs <- se_exp * t(cor.mat[1:r, 1:r] %*% (beta.hat * t(se_exp))) -
        t(cor.mat[1:r, r + 1] * t(se_exp * se_out)) ## size p * r
    ## matrix of u_j
    u.values <- - (var.res * b_exp + res * coefs) / var.res^{3/2}
    S <- t(u.values) %*% u.values

    ## calculate B
    B <- matrix(0, r + 1, r + 1)
    B[1:r, 1:r] <- c1 * S
    B[r + 1, r + 1] <- p * c2

    ## calculate A
    A <- matrix(0, r + 1, r + 1)
    A[r + 1, r+1] <- - c3/2 * sum(1/var.res)
    A[1:r, r+1] <- c4 * colSums(coefs / var.res^2)


    rho.dev.var <- rho(t.values, deriv = 1) / var.res^{3/2} ## size p
    temp <- array(0, c(p, r, r))
    rho.se <- rho.dev.var * se_exp
    for (j in 1:p)
        temp[j, , ] <- res[j] * t(rho.se[j, ] * cor.mat[1:r, 1:r]) * se_exp[j, ]
    temp <- colSums(temp, dim = 1)
    temp1 <- t(rho.dev.var * b_exp) %*% coefs
    A[1:r, 1:r] <- c4 * S  - (2 * temp1 - t(temp1) + temp)

    A.inv <- solve(A)
    asymp.var <- A.inv %*% B %*% t(A.inv)

    if (diagnosis) {
        std.resid <- res / sqrt(var.res)
        qqnorm(std.resid)
        abline(0, 1)
        if (length(std.resid) > 7) {
            require(nortest)
            print(ad.test(std.resid))
            print(shapiro.test(std.resid[sample(length(std.resid),
                                                min(5000, length(std.resid)))]))
        }
    }


    out <- list(beta.hat = beta.hat,
                tau2.hat = tau2.hat,
                beta.var = asymp.var[1:ncol(b_exp), 1:ncol(b_exp)],
                tau2.se = sqrt(asymp.var[nrow(asymp.var), nrow(asymp.var)]), # / sqrt(efficiency),
                beta.p.value = pmin(1, 2 * pnorm(abs(beta.hat) / sqrt(diag(asymp.var)[1:r]),
                                                 lower.tail = F)))# / sqrt(efficiency),
    if (diagnosis) {
        out$std.resid <- std.resid
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

#' L2 loss function and its derivatives
#'
#' @import stats
#' @keywords internal
#'
rho.l2 <- function(r, deriv = 0) {
    if (deriv == 0) {
        return(r^2/2)
    } else if (deriv == 1) {
        return(r)
    } else if (deriv == 2) {
        return(1)
    } else {
        stop("deriv must be 0, 1, or 2.")
    }
}
