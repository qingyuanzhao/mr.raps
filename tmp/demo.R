source("showModes.R")

directory.jw <- "~/Dropbox/two_sample/profile_likelihood_data/"
directory.qz <- "~/Dropbox (Penn)/profile_likelihood_data/"
directory <- directory.qz

setwd(directory)

sel.files <- paste0("ldl/mc_hdl.rda")
exp.files <- paste0("ldl/hdl_2010_gwas.rda")
out.file <- paste0("cad/ukbb_exome_cad.rda")
plink_refdat <- paste0("ld_files/data_maf0.01_rs")

plink_exe.jw <- paste0("ld_files/plink1.90")
plink_exe.qz <- paste0("plink")
plink_exe <- plink_exe.qz

corr <- calCor(sel.files, exp.files, out.file,
               plink_exe, plink_refdat)

dat.list <- getInput(sel.files, exp.files, out.file,
                     plink_exe, plink_refdat, p.thres = 1e-5)
strong.dat.list <- getInput(exp.files, exp.files, out.file, plink_exe, plink_refdat,
                            p.thres = 1e-5,
                            clump_r2 = 0.05, keep.pval0.snps = T)

result <- grappleRobustEst(dat.list$beta_exp,
                           dat.list$data_out$beta,
                           dat.list$se_exp,
                           dat.list$data_out$se, cor.mat = corr,
                           loss.function = "tukey")

modes <- findModes(dat.list$beta_exp,
                   dat.list$data_out$beta,
                   dat.list$se_exp,
                   dat.list$data_out$se,
                   strong.dat.list$beta_exp,
                   strong.dat.list$data_out$beta,
                   strong.dat.list$se_exp,
                   strong.dat.list$data_out$se,
                   cor.mat = corr, mode.lmts = c(-2, 2))
