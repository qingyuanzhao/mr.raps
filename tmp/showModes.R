#' Use the multiple modes of the robust profile likelihood function to find out multiple pathways in MR and their marker SNPs.
#'
#' @param b_exp A vector of length \code{p} for the effect sizes of \code{p} number of independent SNPs on the risk factor
#' @param b_out A vector of length \code{p} for the effect sizes of the \code{p} SNPs on the outcome
#' @param se_exp A vector of length \code{p} for the standard deviations of the effect sizes in \code{b_exp} 
#' @param se_out A vector of length \code{p} for the standard deviations of the effect sizes in \code{b_out} 
#' @param b_exp_st A vector of length \code{q} for the effect sizes of \code{q} candidates SNPs on the risk factor where the marker SNPs are selected from. Default is the same as \code{b_exp}
#' @param b_out_st A vector of length \code{q} for the effect sizes of the \code{q} same list of SNPs as \code{b_exp_st} on the outcome, default is the same as \code{b_out}
#' @param se_exp_st A vector of length \code{q} for the standard deviations of the effect sizes in \code{b_exp_st}, default is the same as \code{se_exp} 
#' @param se_out_st A vector of length \code{q} for the standard deviations of the effect sizes in \code{b_out_st}, default is the same as \code{se_out}
#' @param mode_lmts The range of \code{beta} that the modes are searched from. Default is \code{c(-2, 2)}
#' @param cor.mat Either NULL or a 2 by 2 symmetric matrix. The correlation matrix of estimated effect sizes on the isk factor and the outcome. Default is NULL, for the identity matrix
#' @param loss.function Loss function used, one of "l2", "huber" and "tukey". Default is "tukey". One should use "tukey" or "huber" in order to find multiple modes
#' @param k Tuning parameters of the loss function, for loss "l2", it is NA, for loss "huber", default is 1.345 and for loss "tukey", default is 3. 
#' @param beta.mode Allow providing values of the modes to find out marker SNPs for any given modes
#' @param include.thres Absolute value upper threshold of the standardized test statistics of one SNP on one mode for the SNP to be included as a marker for that mode, default is 1.4
#' @param exclude.thres Absolute value lower threshold of the standardized test statistics of one SNP on other modes for the SNP to be included as a marker for that mode, default is \code{qnorm(0.975)}
#' @param map.marker Whether map each marker to the earist gene or not. Default is TRUE if multiple markers are found. It is always FALSE if there is just one mode.
#'
#' @return A list of the modes, RAP plot and the markers
#'
#' @import ggplot2 
#' @import biomaRt 
#' @import rsnps 
#' @import GenomicRanges 
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene 
#' @import GenomicFeatures
#' @import IRanges
#' @export
findModes <- function(b_exp, b_out, 
                      se_exp, se_out, 
                      b_exp_st = b_exp, b_out_st = b_out, 
                      se_exp_st = se_exp, se_out_st = se_out, 
                      mode.lmts = c(-2, 2),
                      cor.mat = NULL, 
                      loss.function = c("tukey", "l2", "huber"), 
                      k = switch(loss.function[1], 
                                 l2 = NA, huber = 1.345, 
                                 tukey = 3), 
                      beta.mode = NULL,
                      include.thres = 1.4, exclude.thres = qnorm(0.975),
                      map.marker = T) {

  tau2 <- 0
  b_exp <- as.matrix(b_exp)
  se_exp <- as.matrix(se_exp)
  b_exp_st <- as.matrix(b_exp_st)
  se_exp_st <- as.matrix(se_exp_st)

  loss.function <- match.arg(loss.function, c("tukey", "l2", "huber"))
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

  #  delta <- integrate(function(x) rho(x) * dnorm(x), -Inf, Inf)$value
  #  delta <- integrate(function(x) rho(x) * dnorm(x), -Inf, Inf)$value


  ## First, calculate t for each SNP
  t_fun <- function(beta, tau2) {
    upper <- b_out - b_exp %*% beta
    temp <- t(t(cbind(se_exp, se_out)) * c(-beta, 1))
    lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)
    return(upper/lower)
  }

  t_fun_st <- function(beta, tau2) {
    upper <- b_out_st - b_exp_st %*% beta
    temp <- t(t(cbind(se_exp_st, se_out_st)) * c(-beta, 1))
    lower <- sqrt(rowSums((temp %*% cor.mat) * temp) + tau2)
    return(upper/lower)
  }



  robust.optfun.fixtau <- function(beta, tau2) {
    -sum(rho(t_fun(beta, tau2)))
  }

  ## Take 5000 equally spaced points to do grid search to check for modes
  beta.seq <- seq(mode.lmts[1], mode.lmts[2], length.out = 5000)
  val <- sapply(beta.seq, function(beta) robust.optfun.fixtau(beta, tau2))
  mode.pos <- findLocalModes(val)
  temp.data <- data.frame(beta = beta.seq, likelihood = val)


  p <- ggplot2::ggplot(aes(x = beta, y = likelihood), data = temp.data) + geom_line() + 
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = beta.seq[mode.pos], color = "red", linetype = "dashed") + 
  labs(y = "RAP lieklihood") +
  annotate("text", x = mode.lmts[1] * 0.75 + mode.lmts[2] * 0.25,
           y = sum(range(temp.data$likelihood) * c(0.25, 0.75)), label = paste(length(b_out), "SNPs")) + 
      theme(axis.line = element_line(linetype = "solid"),          
      #   axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      #    axis.title.x= element_blank(),
      #   axis.title.y="sss",
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

  if (is.null(beta.mode))
    beta.mode <- beta.seq[mode.pos]

  print(paste("The modes of beta are:", paste(beta.mode, collapse = ",")))
  res.mat <- sapply(beta.mode, function(beta) t_fun_st(beta,0))
  if (length(beta.mode) == 1)
    res.mat <- as.matrix(res.mat)
  colnames(res.mat) <- beta.mode

  rownames(res.mat) <- rownames(b_exp_st)

  ss <- which(rowSums(abs(res.mat) < exclude.thres) == 1 & rowSums(abs(res.mat) < include.thres) > 0)

  markers <- as.data.frame(abs(res.mat[ss, , drop = F]) < include.thres)

  if (ncol(markers) == 1)
    map.marker <- F
  if (map.marker) {

 #   require(biomaRt)

    #    snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
    print("loading ensembl data ...")
 #   data(sysdata)
    #  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    print("ensembl data loaded!")


    snp_ids = rownames(markers)

 #   require(rsnps)
 #   require(GenomicRanges)
    ## find the location of SNPs
    snp_locs0 <- rsnps::ncbi_snp_query(snp_ids)
    snp_locs <- try(GenomicRanges::GRanges(seqnames = paste("chr", snp_locs0$Chromosome, sep = ""), 
                            ranges = IRanges::IRanges(start = snp_locs0$BP, width = 1)))

    if (class(snp_locs) != "try-error") {

  #    require(GenomicFeatures)
  #    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
      ## find the nearest gene
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      refseq.genes<- GenomicFeatures::genes(txdb)
      rr <- data.frame(GenomicRanges::nearest(snp_locs, refseq.genes, ignore.strand = T, select = "all"), 
                       stringsAsFactors = F)
      genes <- refseq.genes[rr[, 2]]
      genes <- data.frame(genes, stringsAsFactors=F)$gene_id
      #  print(genes)

      ## map the gene IDs
      ans <- unique(biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene"),    
                          filters = "entrezgene",
                          values = genes,
                          mart = ensembl) )
      #print(ans)
      rownames(ans) <- ans$entrezgene
      ans <- ans[genes, ]
      rr <- data.frame(rr, ans$hgnc_symbol, genes, stringsAsFactors = F)
      hgnc.names <- sapply(unique(rr[, 1]), 
                           function(k) {
                             temp <- rr[rr[, 1] == k, 3]
                             paste(temp[!is.na(temp)], collapse = "/")})
      entrez.names <-  sapply(unique(rr[, 1]), 
                              function(k) {
                                temp <- rr[rr[, 1] == k, 4]
                                paste(temp[!is.na(temp)], collapse = "/")})
    } else {
      hgnc.names <- rep("NA", nrow(markers))
      entrez.names <- hgnc.names
    }

    markers$nearest_gene <- hgnc.names
    markers$Chromosome <- as.numeric(snp_locs0$Chromosome)
    markers$BP <- as.numeric(snp_locs0$BP)
    markers$inside_gene <- snp_locs0$Gene 
    #   markers <- cbind(markers, as.data.frame(snp_locs0[, -1]))
    markers.full <- markers
    markers.full$entrez_gene <- entrez.names
  } else
    markers.full <- markers

  #   print(markers)

  markers <- markers[order(markers[, 1]), , drop = F]



  return(list(fun = robust.optfun.fixtau,
              modes = beta.seq[mode.pos],
              p = p,
              markers = markers,
            #  markers.full = markers.full,
              res.mat = res.mat))

}


#' findLocal maximum modes of a series of points
#'
#' @keywords internal
#' 
findLocalModes <- function(v, npts = 21) {
  v.temp <- c(rep(0, npts), v, rep(0, npts))

  idx <- c(1:npts, 1:npts + npts + 1)
  #  print(idx)
  judge <- sapply(idx, function(i) v - v.temp[i:(length(v) + i -1)] > 0)
  modes <- rowSums(judge) == ncol(judge)
  return(which(modes))
}




