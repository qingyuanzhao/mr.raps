#' Modal plot to detect heterogeneity
#'
#' @param data Alternatively, dataset can be passed by the argument \code{data}, which must be a data frame with columns \code{beta.exposure}, \code{beta.outcome}, \code{se.exposure}, \code{se.outcome}.
#' @param k Locality of the robust likelihood (smaller \code{k} has more sensitivity for mode detection)
#' @param beta.range range of beta in the plot
#'
#' @export
#'
#' @examples
#' data(lipid.cad)
#' data <- subset(lipid.cad, lipid == "hdl" & restrict &
#' gwas.selection == "teslovich_2010" &
#' gwas.outcome == "cardiogramplusc4d_1000genome" &
#' pval.selection < 1e-5)
#'
#' modal.plot(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, k = 1)
#'
#' data <- subset(lipid.cad, lipid == "ldl" & restrict &
#' gwas.selection == "teslovich_2010" &
#' gwas.outcome == "cardiogramplusc4d_1000genome" &
#' pval.selection < 1e-5)
#'
#' modal.plot(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
#'
modal.plot <- function(b_exp = NULL, b_out = NULL, se_exp = NULL, se_out = NULL, data = NULL, k = 1.5, weight.option = c("MLE", "shrinkage"), beta.range = NULL) {

    weight.option <- match.arg(weight.option)

    rho <- function(r, ...) rho.tukey(r, k, ...)

    ## robust.loglike <- function(beta) {
    ##     - sum(rho((data$beta.outcome - data$beta.exposure * beta) / sqrt(data$se.outcome^2 + data$se.exposure^2 * beta^2)))
    ## }

    if (!is.null(data)) {
        b_exp <- data$beta.exposure
        b_out <- data$beta.outcome
        se_exp <- data$se.exposure
        se_out <- data$se.outcome
    } else {
        if (is.null(b_exp) || is.null(b_out) || is.null(se_exp) || is.null(se_out)) {
            stop("No data input!")
        }
    }

    get.t <- function(beta) {
        (b_out - b_exp * beta) / sqrt(se_out^2 + se_exp^2 * beta^2)
    }

    get.gamma.hat <- function(beta, deriv = c("0", "beta")) {
        deriv <- match.arg(deriv, c("0", "beta"))
        gamma.mle <- (b_exp/ se_exp^2 + beta * b_out/ (se_out^2)) / (1 / se_exp^2 + beta^2 / (se_out^2))
        gamma.mle[se_exp == 0] <- b_exp
        if (deriv == "0") {
            if (shrinkage) {
                a <- posterior.mean(gamma.mle / se_exp,
                                    sqrt(1 / (1 + beta^2 * se_exp^2 / (se_out^2))),
                                    prior.param$p, prior.param$mu, prior.param$sigma)
                gamma.hat <- a * se_exp
            } else {
                gamma.hat <- gamma.mle
            }
            return(gamma.hat)
        }

        if (deriv == "beta") {
            gamma.mle.deriv <- (b_out/ (se_out^2)) / (1 / se_exp^2 + beta^2 / (se_out^2))
        }
        if (shrinkage) {
            shrinkage.deriv <- posterior.mean(
                gamma.mle / se_exp,
                sqrt(1 / (1 + beta^2 * se_exp^2 / (se_out^2))),
                prior.param$p, prior.param$mu, prior.param$sigma,
                deriv = 1)
            return(gamma.mle.deriv * shrinkage.deriv)
        } else {
            return(gamma.mle.deriv)
        }
    }

    psi <- function(beta) {
        if (length(beta) > 1){
            return(sapply(beta, psi))
        }
        t <- get.t(beta)
        gamma.hat <- get.gamma.hat(beta)
        v <- beta^2 * se_exp^2 + se_out^2
        psi1 <- sum(gamma.hat * rho(t, deriv = 1) / sqrt(v))
        return(psi1)
    }

    if (weight.option == "shrinkage") {
        shrinkage <- TRUE
        prior.param <- fit.mixture.model(data$beta.exposure / data$se.exposure)
    } else {
        shrinkage <- FALSE
    }

    fit <- mr.raps.mle(b_exp, b_out, se_exp, se_out)
    if (is.null(beta.range)) {
        beta.min <- fit$beta.hat - 25 * fit$beta.se
        beta.max <- fit$beta.hat + 25 * fit$beta.se
    } else {
        beta.min <- beta.range[1]
        beta.max <- beta.range[2]
    }
    beta.seq <- c(-Inf, seq(beta.min, beta.max, length = 50))

    rll <- rep(0, length(beta.seq) - 1)
    for (i in 2:length(beta.seq)) {
        rll[i-1] <- integrate(psi, beta.seq[i - 1], beta.seq[i])$value
    }

    plot(beta.seq[-1],
         cumsum(rll),
         type = "l",
         xlab = expression(beta),
         ylab = "Robust log-likelihood",
         main = paste("k =", k))

}
