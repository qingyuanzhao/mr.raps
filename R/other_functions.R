#' @keywords internal
#' 
plink_clump <- function(dat,
                        plink_exe, 
                        refdat,
                        clump_kb = 10000,
                        clump_r2 = 0.001,
                        clump_p1 = 1,
                        clump_p2 = 1,
                        tempdir = "temp")
{
                                        # Make textfile
    snps <- dat$SNP
    if (!("pval" %in% colnames(dat)))
        dat$pval <- dat$pval.exposure
    pvals <- dat$pval
    dir.create(file.path(tempdir), showWarnings = FALSE)

    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    fn <- tempfile(tmpdir = tempdir)
    require(data.table)
    fwrite(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F, sep = " ")

    fun2 <- paste0(
      # shQuote(plink_exe),
        plink_exe,
        # shQuote("plink"),
        #    " --bfile ", shQuote(refdat, type=shell),
    #    " --clump ", shQuote(fn, type=shell),
        " --bfile ", refdat,
        " --clump ", fn,
        " --clump-p1 ", clump_p1,
        " --clump-p2 ", clump_p2,
        " --clump-r2 ", clump_r2,
        " --clump-kb ", clump_kb,
        " --out ", shQuote(fn, type=shell)
    )
    system(fun2)
    a <- fread(paste(fn, ".clumped", sep=""), he=T)
    ## a <- fread(paste(fn, sep=""), he=T)
    unlink(paste(fn, "*", sep=""))
    a <- a[, c(3, 5)]
    a$temp.p <- round(log10(a$P))
    dat$temp.p <- round(log10(dat$pval))
    a <- merge(a, dat, by = c("SNP", "temp.p"))

    a <- a[, -(2:3)]

    return(a)
}


#' @keywords internal
#'
formatData <- function(dat, cls = c("exposure", "outcome")) {
    cls <- match.arg(cls, c("exposure", "outcome"))
    for (nm in c("eaf", "id")) {
        if (!(nm %in% colnames(dat)))
            dat[, nm] <- rep("NR", nrow(dat))
    }
    idx <- which(colnames(dat) == "SNP")
    colnames(dat)[-idx] <- paste(colnames(dat)[-idx], ".", cls, sep = "")
    return(dat)
}


