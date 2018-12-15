#' Scatter plot with annotation
#'
#' @param data A data frame (see \code{\link{mr.raps}}).
#' @param annotate Annotating the points? (rsid, chromosome, position)
#' @param annotate.genes Further annotation of closest genes? See example.
#' @param rank.method How to select strongest SNPs for plot?
#' @param num.snps How many SNPs are shown?
#' @param fit A \code{mr.raps} fit.
#'
#' @import rsnps
#' @import ggplot2 ggrepel
#'
#' @examples
#' data(bmi.sbp)
#' mr.raps.scatterplot(bmi.sbp)
#'
#' \donttest{
#' require(bumphunter)
#' require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' mr.raps.scatterplot(bmi.sbp, annotate.genes = genes)
#' }
#'
#' @export
#'
mr.raps.scatterplot <- function(data, annotate = TRUE, annotate.genes = NULL, rank.method = c("pval.both", "pval.selection", "pval.exposure"), num.snps = 10, fit = mr.raps(data, FALSE)) {

    rank.method <- match.arg(rank.method)

    if (!"pval.exposure" %in% names(data)) {
        data$pval.exposure <- 2 * pnorm(-abs(data$beta.exposure / data$se.exposure))
    }

    if (!"pval.exposure" %in% names(data)) {
        message("Use rank.method = pval.exposure because pval.selection is not available")
        rank.method <- "pval.exposure"
    }

    r <- switch(rank.method,
                pval.both = rank(data$pval.selection * data$pval.exposure),
                pval.selection = rank(data$pval.selection),
                pval.exposure = rank(data$pval.exposure))
    data <- data[r <= num.snps, ]

    if (annotate) {

        info <- ncbi_snp_query(data$SNP)

        snps <- data.frame(chr = paste0("chr", info$Chromosome),
                           start = info$BP,
                           end = info$BP)
        if (!is.null(annotate.genes)) {
            tab <- bumphunter::matchGenes(snps, genes)
            info <- cbind(info, tab)
        }

        data <- merge(data, info, by.x = "SNP", by.y = "Marker")
    }

    data$beta.outcome <- data$beta.outcome * sign(data$beta.exposure)
    data$beta.exposure <- abs(data$beta.exposure)

    p <- ggplot(data) + aes(x = beta.exposure, y = beta.outcome,
                            xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure, ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome) + geom_point(color = "#F8766D", size = 2) + geom_errorbar(alpha = 0.3, width = 0, color = "#F8766D") + geom_errorbarh(alpha = 0.3, height = 0, color = "#F8766D") + expand_limits(x = 0, y = 0) + xlab(paste("SNP effect on exposure")) + ylab(paste("SNP effect on outcome")) + theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4)

    p <- p + geom_abline(intercept = 0, slope = fit$beta.hat, alpha = 0.3, color = "#F8766D", size = 1, linetype = "solid")

    if (annotate) {
        if (!is.null(annotate.genes)) {
            p <- p + aes(label = paste(SNP, name, sep = "\n")) + geom_text_repel(alpha = 0.6, size = 3.5, force = 5)
        } else {
            p <- p + aes(label = paste(SNP, paste(Chromosome, BP, sep = ":"), sep = "\n")) + geom_text_repel(alpha = 0.6, size = 3.5, force = 5)
        }
    }

    p

}
