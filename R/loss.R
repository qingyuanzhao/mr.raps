#' Squared error loss and its derivatives
#'
#' @import stats
#' @keywords internal
#'
rho.l2 <- function(r, k = 2, deriv = 0) {
    if (deriv == 0) {
        r^2 / k
    } else if (deriv == 1) {
        r
    } else if (deriv == 2) {
        1
    }
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
