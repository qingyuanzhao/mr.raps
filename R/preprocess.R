#' SNP selection and clumping based on GWAS summary statistics
#'
#'
#' @param sel.files A vector of length \code{k} of the GWAS summary statistics file names of the \code{k} risk factors SNP selection. Each GWAS file is an ".rda" or ".rData" file containing an R object of name "dat". The "dat" object is a data frame that at least has a column "SNP" for the SNP ids, "pval" for the p-values.
#' @param exp.files A vector of length \code{k} of the GWAS summary statistics file names of the \code{k} risk factors for getting the effect sizes and standard deviations. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele of the SNP (capitalized letters) and "other_allele" for the other allele of the SNP (capitalized letters).
#' @param out.file The GWAS summary statistics file name for the outcome data. Each GWAS file should have a column "SNP" for the SNP ids, "beta" for the effect sizes, "se" for the standard deviation, "effect_allele" for the effect allele of the SNP (capitalized letters) and "other_allele" for the other allele of the SNP (capitalized letters).
#' @param p.thres The upper threshold of the selection p-values for a SNP to be selected before clumping. It only requires that at least one of the p-values of the risk factors of the SNPs to be below the threshold.
#' @param clump.directly Whether the clump is done directly on the union of selection files SNPs, or the intersection of the SNPs in any of the selection files, the SNPs in any of the exposure files and the SNPs in the outcome GWAS data file. Default is FALSE
#' @param clump_r2 The clumping r2 threshold in PLINK. Default is set to 0.001 for selection of independent SNPs, and can be set to a higher value when select candidate SNPs for the marker SNPs.
#' @param keep.pval0.snp Whether keep all the SNPs whose minimum selection p-values is 0. Useful for selection of marker SNPs when set to TRUE. Default is FALSE
#' @param plink_exe The excutable file of PLINK
#' @param plink_refdat The reference files (.bed, .bim, .fam) for PLINK
#'
#' @return A list of selected summary statistics
#' @export
getInput <- function(sel.files,
                     exp.files,
                     out.file,
                     plink_exe, plink_refdat,
                     p.thres = 0.01, clump.directly = F,
                     clump_r2 = 0.001, keep.pval0.snps = F) {
    sel.SNPs <- NULL
    pvals <- NULL
    for (file in sel.files) {
        print(paste("loading data for selection:", file, "..."))
        load(file)
        sel.snps <- dat$SNP[dat$pval < p.thres]


        if (is.null(sel.SNPs)) {
            sel.SNPs <- sel.snps
            pvals <- dat$pval[dat$pval < p.thres]
            names(pvals) <- sel.SNPs
        } else {
            ori.snps <- sel.SNPs
            sel.SNPs <- union(sel.SNPs, sel.snps)
            temp <- rep(1, length(sel.SNPs))
            names(temp) <- sel.SNPs
            temp[ori.snps] <- pvals
            temp1 <- rep(1, length(sel.SNPs))
            names(temp1) <- sel.SNPs
            temp1[sel.snps] <- dat$pval[dat$pval < p.thres]
            pvals <- pmin(temp, temp1)
                                        #  pvals <- temp
            names(pvals) <- sel.SNPs
            rm(temp)
        }
        rm(dat)
    }

    pvals <- pmin(pvals * length(sel.files), 1) # PLINK reports error when p-value > 1

    if (clump.directly) {
        data.sel <- data.frame(SNP = sel.SNPs, pval = pvals)
        paste("Clumping using PLINK")
        data.sel <- plink_clump(data.sel, plink_exe, plink_refdat, clump_r2 = clump_r2)
        sel.SNPs <- data.sel$SNP
        pvals <- data.sel$pval
        names(pvals) <- sel.SNPs
        rm(data.sel)
        gc()
    }


    beta_exp <- NULL
    se_exp <- NULL
    ref.data.exp <- NULL
    for (exp.file in exp.files) {
        print(paste("loading data from exposure:", exp.file, "..."))

        load(exp.file)
        sel.SNPs <- intersect(sel.SNPs, dat$SNP)
        temp <- dat[dat$SNP %in% sel.SNPs, ]
        rm(dat)


        ## harmonize
        if (!is.null(ref.data.exp)) {
            temp <- formatData(temp, "outcome")
            temp.result <- suppressMessages(harmonise_data(ref.data.exp, temp))
            temp.result <- temp.result[temp.result$mr_keep, ]
            temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
            colnames(temp) <- gsub(".outcome", "", colnames(temp))
            ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
                                        #    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
            rm(temp.result)
            gc()
            beta_exp <- beta_exp[as.character(ref.data.exp$SNP), , drop = F]
            se_exp <- se_exp[as.character(ref.data.exp$SNP), , drop = F]
            flip <- rep(1, nrow(beta_exp))
            flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
                                        #    print(sum(flip == -1))
            beta_exp <- flip * beta_exp
            sel.SNPs <- rownames(beta_exp)
        } else {
            ref.data.exp <- formatData(temp, "exposure")
        }

        if (is.null(beta_exp)) {
            beta_exp <- data.frame(temp$beta)
            rownames(beta_exp) <- temp$SNP
            se_exp <- data.frame(temp$se)
            rownames(se_exp) <- temp$SNP
        } else {
            beta_exp <- beta_exp[as.character(temp$SNP), , drop = F]
            se_exp <- se_exp[as.character(temp$SNP), , drop = F]
            beta_exp <- cbind(beta_exp, temp$beta)
            se_exp <- cbind(se_exp, temp$se)
        }
        rm(temp)
    }

    ## load outcome file
    print(paste("loading data for outcome:", out.file, "..."))
    load(out.file)
    sel.SNPs <- intersect(sel.SNPs, dat$SNP)
    temp <- dat[dat$SNP %in% sel.SNPs, ]
    rm(dat)

    ## harmonize
    temp <- formatData(temp, "outcome")
    temp.result <- suppressMessages(harmonise_data(ref.data.exp, temp))
    temp.result <- temp.result[temp.result$mr_keep, ]
    temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
    colnames(temp) <- gsub(".outcome", "", colnames(temp))
    ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
                                        #    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
    rm(temp.result)
    gc()
    beta_exp <- beta_exp[as.character(ref.data.exp$SNP), , drop = F]
    se_exp <- se_exp[as.character(ref.data.exp$SNP), , drop = F]

    flip <- rep(1, nrow(beta_exp))
    flip[sign(beta_exp[, 1]) != sign(ref.data.exp$beta.exposure)] <- -1
                                        #  print(sum(flip == -1))
    beta_exp <- flip * beta_exp
    sel.SNPs <- rownames(beta_exp)

    data_out <- temp
    rownames(data_out) <- data_out$SNP
    pvals <- pvals[sel.SNPs]

    if (!clump.directly) {
                                        #  print(quantile(pvals))

        if (keep.pval0.snps) {
            idx <- which(pvals == 0)
            keep.snps <- sel.SNPs[idx]
        } else {
            keep.snps <- c()
        }
        data.sel <- data.frame(SNP = sel.SNPs, pval = pvals)


        data.sel <- plink_clump(data.sel, plink_exe, plink_refdat, clump_r2 = clump_r2)
        sel.SNPs <- union(as.character(keep.snps), as.character(data.sel$SNP))
        beta_exp <- beta_exp[as.character(sel.SNPs), , drop = F]
        se_exp <- se_exp[as.character(sel.SNPs), , drop = F]
        data_out <- data_out[as.character(sel.SNPs), , drop = F]

        pvals <- pvals[sel.SNPs]
        names(pvals) <- sel.SNPs
        rm(data.sel)
        gc()
    }

    return(list(beta_exp = beta_exp, se_exp = se_exp,
                data_out = data_out, sel.pvals = pvals))
}




#' Calculate the \code{(k + 1)} by \code{(k + 1)} correlation matrix between the \code{k} risk factors and the outcome
#'
#' We use the SNPs that has no/very weak effect on the risk factors to estimate the shared correlation matrix acrooss SNPs
#'
#' @inheritParams getInput
#' @param p.thres The lower threshold of the p-values for a SNP to be used in calculating the correlation matrix
#' @param seed.vec A length \code{m} vector of seeds used. In order to get a stable estimate of the correlation matrix, we take the average of \code{m} numbers of estimated correlation matrix weather each estimate use a randomly sampled clumped set of SNPs.
#'
#' @export
calCor <- function(sel.files,
                   exp.files,
                   out.file,
                   plink_exe, refdat,
                   p.thres = 0.5, clump.directly = F,
                   seed.vec = 1:20) {
    sel.SNPs <- NULL
    for (file in sel.files) {
        print(paste("loading data from selection:", file, "..."))
        load(file)
        sel.snps <- dat$SNP[dat$pval > p.thres]
        if (is.null(sel.SNPs))
            sel.SNPs <- sel.snps
        else
            sel.SNPs <- intersect(sel.SNPs, sel.snps)
        rm(dat)
    }

    if (clump.directly) {
        data.sel <- data.frame(SNP = sel.SNPs, pval = runif(length(sel.SNPs)))
        paste("Clumping using PLINK")
        data.sel <- plink_clump(data.sel, plink_exe, refdat)
        sel.SNPs <- data.sel$SNP
        rm(data.sel)
        gc()
    }


    data.exp <- NULL
    ref.data.exp <- NULL
    for (exp.file in exp.files) {
        print(paste("loading data from exposure:", exp.file, "..."))

        load(exp.file)
        sel.SNPs <- intersect(sel.SNPs, dat$SNP)
        temp <- dat[dat$SNP %in% sel.SNPs, ]
        rm(dat)

        ## harmonize
        if (!is.null(ref.data.exp)) {
            temp <- formatData(temp, "outcome")
            temp.result <- suppressMessages(harmonise_data(ref.data.exp, temp))
            temp.result <- temp.result[temp.result$mr_keep, ]
            temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
            colnames(temp) <- gsub(".outcome", "", colnames(temp))
            ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
                                        #    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
            rm(temp.result)
            gc()
            data.exp <- data.exp[as.character(ref.data.exp$SNP), , drop = F]
            flip <- rep(1, nrow(data.exp))
            flip[sign(data.exp[, 1]) != sign(ref.data.exp$beta.exposure/ref.data.exp$se.exposure)] <- -1
                                        #   print(sum(flip == -1))
            data.exp <- flip * data.exp
            sel.SNPs <- rownames(data.exp)
        } else {
            ref.data.exp <- formatData(temp, "exposure")
        }

        temp$z <- temp$beta/temp$se
        if (is.null(data.exp)) {
            data.exp <- data.frame(temp$z)
            rownames(data.exp) <- temp$SNP
        } else {
            data.exp <- data.exp[as.character(temp$SNP), , drop = F]
            data.exp <- cbind(data.exp, temp$z)
        }
        rm(temp)
    }
    colnames(data.exp) <- paste("exposure", 1:length(exp.files))

    ## load outcome file
    print(paste("loading data for outcome:", out.file, "..."))
    load(out.file)
    sel.SNPs <- intersect(sel.SNPs, dat$SNP)
    temp <- dat[dat$SNP %in% sel.SNPs, ]
    rm(dat)
    ## harmonize
    temp <- formatData(temp, "outcome")
    temp.result <- suppressMessages(harmonise_data(ref.data.exp, temp))
    temp.result <- temp.result[temp.result$mr_keep, ]
    temp <- temp.result[, c(1, grep(".outcome", colnames(temp.result)))]
    colnames(temp) <- gsub(".outcome", "", colnames(temp))
    ref.data.exp <- temp.result[, c(1, grep(".exposure", colnames(temp.result)))]
                                        #    colnames(data.sel) <- gsub(".exposure", "", colnames(data.sel))
    rm(temp.result)
    gc()
    data.exp <- data.exp[as.character(ref.data.exp$SNP), , drop = F]
    flip <- rep(1, nrow(data.exp))
    flip[sign(data.exp[, 1]) != sign(ref.data.exp$beta.exposure/ref.data.exp$se.exposure)] <- -1
                                        #  print(sum(flip == -1))
    data.exp <- flip * data.exp
    sel.SNPs <- rownames(data.exp)

    temp$z <- temp$beta/temp$se
    data.exp <- data.exp[as.character(temp$SNP), , drop = F]
    data.exp <- cbind(data.exp, temp$z)
    colnames(data.exp) <- c(colnames(data.exp)[1:(ncol(data.exp) - 1)], "Outcome")
    rm(temp)
    gc()


    if (!clump.directly) {
        corr.list <- lapply(seed.vec, function(seed) {
            set.seed(seed)
            print(paste("Seed is:", seed))
            data.sel <- data.frame(SNP = sel.SNPs, pval = runif(length(sel.SNPs)))
            data.sel <- suppressMessages(plink_clump(data.sel, plink_exe, refdat))
            data.exp1 <- data.exp[as.character(data.sel$SNP), , drop = F]
            rm(data.sel)
            corr <- cor(data.exp1)
            rm(data.exp1)
            gc()
            return(corr)
        })
        corr <- Reduce("+", corr.list)/length(corr.list)
    } else
        corr <- cor(data.exp)


    return(corr)
}
