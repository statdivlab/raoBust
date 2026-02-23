#' Presence and absence of COG 4251 with function Bacteriophytochrome (light-regulated signal transduction histidine kinase) within a set of *Prochlorococcus* genomes
#'
#' @format A data frame object with 31 observations of 5 covariates.
#' \describe{
#' \item{sample}{Sample id}
#' \item{cog}{COG, all observations are from COG 4251}
#' \item{presence}{1 if COG 4251 is present in the sample, 0 otherwise}
#' \item{clade}{light level, combined with category I or II}
#' \item{light}{HL for high light condition, LL for low light condition}
#' }
#' @references https://statdivlab.github.io/blog/articles/ebame10_regression_pt2.html
"light_presence"