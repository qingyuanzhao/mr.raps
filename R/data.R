#' Effect of Body Mass Index (BMI) on Systolic Blood Pressure (SBP)
#'
#' Summary data obtained by combining three genome-wide association studies:
#' \enumerate{
#' \item{BMI-FEM:}{ BMI in females by the Genetic Investigation of ANthropometric Traits (GIANT) consortium (sample size: 171977).}
#' \item{BMI-MAL:}{ BMI in males in the same study by the GIANT consortium (sam- ple size: 152893)}
#' \item{SBP-UKBB:}{ SBP using the United Kingdom BioBank (UKBB) data (sample size: 317754)}
#' }
#'
#' The BMI-FEM dataset is used for SNP selection (column \code{pval.selection}). The BMI-MAL dataset estimates the SNPs' effect on BMI and the SBP-UKBB dataset estimates the SNPs' on SBP.
#'
#' @docType data
#'
#' @usage data(bmi.sbp)
#'
#' @format A \code{data.frame} with 160 rows and 29 variables.
#'
#' @keywords datasets
#'
#'
"bmi.sbp"

#' "Effect" of Body Mass Index (BMI) on Body Mass Index (BMI)
#'
#' Summary data obtained by combining three genome-wide association studies:
#' \enumerate{
#' \item{BMI-GIANT:}{ BMI in the Genetic Investigation of ANthropometric Traits (GIANT) consortium (sample size: 339224).}
#' \item{BMI-UKBB-1:}{ BMI in a half of the United Kingdom BioBank (UKBB) data (sample size: 234070)}
#' \item{SBP-UKBB-2:}{ BMI in the other half of the UKBB data (sample size: 234070)}
#' }
#'
#' The BMI-GIANT dataset is used for SNP selection (column \code{pval.selection}). The BMI-UKBB-1 dataset estimates the SNPs' effects on BMI (columns \code{beta.exposure} and \code{se.exposure}) and the BMI-UKBB-2 dataset provides independent estimates of the same effects (columns \code{beta.outcome} and \code{se.outcome}).
#'
#' @docType data
#'
#' @usage data(bmi.bmi)
#'
#' @format A \code{data.frame} with 812 rows and 28 variables.
#'
#' @keywords datasets
#'
#'
"bmi.bmi"

#' Effect of blood lipids (LDL cholesterol, HDL cholesterol, Triglycerides) on Cardiovascular disease risk
#'
#' This dataset is created from four genome-wide association studies:
#' \enumerate{
#' \item A 2010 GWAS of blood lipids (Teslovich et al.\, 2010), named "teslovich_2010" in the dataset.
#' \item The MetaboChip (MC) data in a 2013 GWAS of blood lipids (Willer et al.\, 2013), named "mc_2013" in the dataset.
#' \item The CARDIoGRAMplusC4D meta-analysis of coronary artery disease (CARDIoGRAMplusC4D Consortium, 20135, named "cardiogramplusc4d_1000genome" in the dataset.
#' \item The UK BioBank GWAS of self reported heart attach (interim release by the Neale lab), named "ukbb_6150_round2" in the dataset.
#' }
#'
#' \code{lipid.cad} contains in total 24 sub-datasets, each is suitable for a Mendelian randomization study. To obtain a sub-dataset, you must decide on
#' \describe{
#' \item{lipid}{Which lipid trait to consider? Either \code{ldl}, \code{hdl}, or \code{tg}.}
#' \item{gwas.selection}{Which GWAS is used for selection? Either \code{teslovich_2010} or \code{mc_2013}.}
#' \item{gwas.exposure}{Which GWAS is used for exposure? Either \code{teslovich_2010} or \code{mc_2013} and must be different from \code{gwas.selection}.}
#' \item{gwas.outcome}{Which GWAS is used for outcome? Either \code{cardiogramplusc4d} or \code{ukbb_self_report_heart}.}
#' \item{restrict}{Should we use SNPs that are not associated with the other lipids? For example, if we are studying the effect of HDL cholesterol (so \code{lipid} is "hdl") and \code{restrict} is TRUE, then the SNPs are not associated with LDL cholesterol and triglycerides (p-value > 0.01 in the \code{gwas.selection} data).}
#' }
#'
#' @docType data
#'
#' @usage data(lipid.cad)
#'
#' @format A \code{data.frame} with 12026 rows and 24 variables.
#'
#' @keywords datasets
#'
#'
"lipid.cad"
