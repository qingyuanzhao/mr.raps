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
#' @format A data.frame.
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
#' @format A data.frame.
#'
#' @keywords datasets
#'
#'
"bmi.bmi"
