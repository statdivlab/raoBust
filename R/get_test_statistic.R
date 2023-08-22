#' Just in case anyone wants to invert a p-value to recover the chi-squared distributed test statistic
#'
#' @param pvalue the p-value to invert
#' @param df the degrees of freedom for the test
#'
#' @importFrom stats qchisq
#'
#' @export
#'
get_test_statistic <- function(pvalue, df) {

  qchisq(1 - pvalue, df)


}
