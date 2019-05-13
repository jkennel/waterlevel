#' predict_terms
#'
#' predict by term for 'mlm' or 'lm' object
#'
#' @param object model fit from lm object
#' @param ... currently not used
#'
#' @return return matrix of predictions
#' 
#' @importFrom stats terms model.matrix
#' 
#' @export
#'
predict_terms <-
  function(object, ...)
  {
    tt <- terms(object)
    if(!inherits(object, "lm"))
      warning("calling predict.lm(<fake-lm-object>) ...")
    
    mm <- X <- model.matrix(object)
    offset <- object$offset
    
    n <- length(object$residuals) # NROW(qr(object)$qr)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- if(p) qr(object)$pivot[p1]
    beta <- object$coefficients
    
    aa <- attr(mm, "assign")
    ll <- attr(tt, "term.labels")
    hasintercept <- attr(tt, "intercept") > 0L
    if (hasintercept) ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
    asgn <- split(order(aa), aaa)
    if (hasintercept) {
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(mm)
      termsconst <- colSums(avx[piv] * beta[piv, , drop = FALSE])
    }
    nterms <- length(asgn)
    if(nterms > 0) {
      predictor <- matrix(ncol = nterms * ncol(beta), nrow = NROW(X))
      if (ncol(beta) > 1) {
        cn <- apply(expand.grid(colnames(beta), names(asgn)), 1, paste, collapse = '_') #jrk
      } else {
        cn <- names(asgn)
      }
      dimnames(predictor) <- list(rownames(X), cn)
      if(hasintercept)
        X <- sweep(X, 2L, avx, check.margin=FALSE)
      unpiv <- rep.int(0L, NCOL(X))
      unpiv[piv] <- p1
      ## Predicted values will be set to 0 for any term that
      ## corresponds to columns of the X-matrix that are
      ## completely aliased with earlier columns.
      for (i in seq.int(1L, nterms, length.out = nterms)) {
        iipiv <- asgn[[i]]      # Columns of X, ith term
        ii <- unpiv[iipiv]      # Corresponding rows of Rinv
        iipiv[ii == 0L] <- 0L
        inds <- 1:ncol(beta) + (i-1) * ncol(beta)             #jrk
        predictor[, inds] <-                                  #jrk
          if(any(iipiv > 0L)) X[, iipiv, drop = FALSE] %*% beta[iipiv, , drop = FALSE]
        else 0
      }
      attr(predictor, 'constant') <- if (hasintercept) termsconst else 0
    }
    
    return(predictor)
  }


