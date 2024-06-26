#' Jackknife standard errors 
#'
#' @param object The fitted object under the alternative.
#' @param dat The data used to fit the model
#' @param id Observations with the same id are in the same cluster. If not included, independence between 
#' observations is assumed.
#'
#' @export
jackknife_se <- function(object, dat, id = NULL){
  if (is.null(id)) {
    id <- 1:nrow(dat)
  }
  if (!is.numeric(id)) {
    id <- as.factor(id)
  }
  unique_id <- unique(id)
  parm <- sapply(unique_id, update_model, object, dat, id)
  parm <- t(parm)
  parm.mean <- apply(parm, 2, mean)
  
  parm.cent <- sapply(1:nrow(parm),
                      function(i) {
                        parm[i, ] - parm.mean
                      })
  parm.cent <- t(parm.cent) 
  
  jack.var <- ((length(unique_id) - 1) / length(unique_id)) * t(parm.cent) %*% parm.cent
  jack.se <- sqrt(diag(jack.var))
  jack.se
}

update_model <- function(id_i, object, dat, id) {
  remove_ind <- which(id == id_i)
  dat.i <- dat[-remove_ind, ]
  object$call <- call_modify(object$call, id = id[-remove_ind])
  coef(stats::update(object, data = dat.i))
}
