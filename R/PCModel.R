## Function to fit a partial credit modell via nlm.  See e.g. Masters & Wright (1984)
PCModel.fit <- function (y, weights = NULL, nullcats = c("keep", "downcode", "ignore"),
                         start = NULL, reltol = 1e-10, deriv = c("sum", "diff"),
                         hessian = TRUE, maxit = 100L, full = TRUE, ...)
{
  ## argument matching
  deriv <- match.arg(deriv)

  ## process input and data
  nullcats <- match.arg(nullcats)
  diff <- deriv == "diff"
  y <- as.matrix(y)
  n_org <- nrow(y)
  m <- ncol(y)
  ident_items <- rep.int(TRUE, m)

  ## process weights
  if (is.null(weights)) weights <- rep.int(1L, n_org) else stopifnot(length(weights) == n_org)
  y <- y[weights > 0, , drop = FALSE]
  weights_org <- weights
  weights <- weights[weights > 0]
  n <- nrow(y)

  ## move categories to zero if necessary, get number of categories per item
  mincat <- apply(y, 2, min, na.rm = TRUE)
  if(any(mincat > 0)) {
    warning("Minimum score is not zero for all items (", paste("I", which(mincat > 0), sep = "", collapse = ", "), "). These items are scaled to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }
  oj <- apply(y, 2, max, na.rm = TRUE)
  oj_vec <- lapply(oj, seq, from = 0)

  ## check for null categories and warn...
  nl_cats <- lapply(1:m, function (j) !(oj_vec[[j]] %in% y[, j]))
  nl_cats_items <- which(sapply(nl_cats, any))
  any_nl_cat_item <- length(nl_cats_items) > 0
  if (any_nl_cat_item) {
    warning("There are items with null categories (", paste("I", nl_cats_items, sep = "", collapse = ", "), ").")
    ## .. then treat according to nullcats (ignore = do nothing)
    if (nullcats == "downcode") {
      oj <- oj - sapply(nl_cats, sum)
      oj_vec <- lapply(1:m, function (j) oj_vec[[j]][1:(oj[j]+1)])
      for (j in nl_cats_items) {
        missing_scores <- sort(which(nl_cats[[j]]), decreasing = TRUE) - 1
        for (score in missing_scores) y[!is.na(y[, j]) & y[, j] > score, j] <-  y[!is.na(y[, j]) & y[, j] > score, j] - 1
      }
    }
    if (nullcats == "keep") {
      all_par <- rep.int(NA, sum(oj))
      est_par <- !unlist(lapply(nl_cats, "[", -1))
      if (est_par[1]) { ## if first parameter not null, restrict it
        all_par[1] <- 0
        est_par[1] <- FALSE        
      } else { ## otherwise restrict next parameter without null categories
        rp <- min(which(est_par))
        all_par[rp] <- 0
        est_par[rp] <- FALSE
      }
      oj <- oj - sapply(nl_cats, sum)
      oj_vec <- mapply(function(x, y) x[!y], oj_vec, nl_cats)
    }
  }
  oj_max <- sapply(oj_vec, max)         # maximum category number (different from number of categories if nullcats == "keep")

  ## remove unidentified items (only one category used)
  unid <- oj == 0
  nunid <- sum(unid)
  if (nunid) {
    y <- y[, !unid, drop = FALSE]
    oj <- oj[!unid]
    oj_vec <- oj_vec[!unid]
    oj_max <- oj_max[!unid]
    m <- m - nunid
    ident_items[unid] <- FALSE
    warning("There were unidentified items (only one category used, ", paste("I", which(unid), sep = "", collapse = ", "), ").")
  }
  mv <- 1:m

  ## check for missings
  y_na <- is.na(y)
  any_na <- any(y_na)

  ## calculate category totals and set starting values from these (calculation accounts for NAs)
  ## (generalisation of first proposal of Fischer & Molenaar, p. 50. If oj = 1: equal to starting values in RaschModel.fit())
  ctot <- vector("list", length = m)
  for (j in seq_len(m)) ctot[[j]] <- as.vector(tapply(weights, factor(y[, j], levels = oj_vec[[j]]), sum))
  if (is.null(start)) {
    start <- lapply(ctot, function (x) - cumsum(diff.default(log(x)))) # delta_jk = log(x_j(k-1) - log(x_jk), beta_jk = cumsum_k=1^k(delta_jk)
    start <- unlist(start)
    start <- start[-1] - start[1]
    start[is.na(start)] <- 0
  }

  ## calculate person totals, catch missing person scores
  rs <- rowSums(y, na.rm = TRUE)
  ptot <- as.vector(tapply(weights, factor(rs, levels = 0:sum(oj_max)), sum))
  ptot[is.na(ptot)] <- 0

  ## unlist ctot, remove first category total (since beta_j0 = 0 for all 0), catch null categories
  ctot <- unlist(lapply(ctot, "[", -1))
  ctot[is.na(ctot)] <- 0
    
  ## build log-likelihood
  if(!any_na) {
    
    ## objective function: conditional log-likelihood
    cloglik <- function (par) {

      ## include epsilons for esf when there are nullcats & nullcats == "keep" (see Wilson, 1993)
      if (any_nl_cat_item && nullcats == "keep") {
        all_par[est_par] <- par
        esf_par <- all_par
      } else {
        esf_par <- c(0, par)
      }
      esf_par <- split(esf_par, rep.int(mv, oj_max))

      ## finally: the cll
      cll <- - sum(ctot * c(0, par)) - sum(ptot * log(elementary_symmetric_functions(par = esf_par, order = 0, diff = diff)[[1]]))

      ## catch degenerated cases (typically caused by non-finite gamma)
      if (is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax
      
      return(-cll)
    }

    ## static helper variables for analytical gradient
    parindex <- unlist(lapply(oj_vec, "[", -1))
    itemindex <- rep.int(mv, oj)
    xmat <- matrix(FALSE, nrow = n, ncol = sum(oj))
    for (i in 1:n) xmat[i, ] <- y[i, itemindex] == parindex 

    ## analytical gradient
    agrad <- function (par) {

        ## elementary symmetric functions
        if (any_nl_cat_item && nullcats == "keep") {
            all_par[est_par] <- par
            esf_par <- all_par
        } else {
            esf_par <- c(0, par)
        }
        esf <- elementary_symmetric_functions(par = split(esf_par, rep.int(mv, oj_max)), order = 1, diff = diff)
        gamma0 <- esf[[1]][rs + 1]
        gamma1 <- apply(esf[[2]], 2, "[", rs + 1)
        if (any_nl_cat_item && nullcats == "keep") gamma1 <- gamma1[, est_par, drop = FALSE]

        ## calculate and return aggregated gradient
        return(- colSums((weights * (- xmat + (gamma1 / gamma0)))[, -1, drop = FALSE]))
    }

  } else {
    
    ## fetch different NA-patterns like Achim does, setup position vector of parameters
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    lev_na_patterns <- levels(na_patterns)
    par_pos <- rep.int(mv, oj_max)

    ## setup na_pattern lists for various elements of loglik
    m_i <- mv_i <- rs_i <- ptot_i <- ctot_i <- par_i <- oj_max_i <- vector("list", length(lev_na_patterns))
    na_obs_i <- weights_i <- est_par_i <- xmat_i <- vector("list", length(lev_na_patterns))

    ## loop over observed NA patterns, calculate constant things once
    for(i in seq_along(lev_na_patterns)) {

      ## from pattern i: get NA item(s), observations of and weights with this pattern
      na_items_i <- as.integer(strsplit(lev_na_patterns[i], "\r")[[1]])
      n_na_items_i <- length(na_items_i)
      na_obs_i[[i]] <- which(na_patterns == lev_na_patterns[i])
      weights_i[[i]] <- weights[na_obs_i[[i]]]
        
      ## select subset
      if(n_na_items_i < 1) {            # no missings
        y_i <- y[na_obs_i[[i]], , drop = FALSE]
        par_i[[i]] <- rep.int(TRUE, sum(oj_max))
        m_i[[i]] <- m
        mv_i[[i]] <- mv
        oj_vec_i <- oj_vec
        oj_max_i[[i]] <- oj_max
      } else {                          # specific NA pattern
        y_i <- y[na_obs_i[[i]], -na_items_i, drop = FALSE]
        par_i[[i]] <- !(par_pos %in% na_items_i)
        m_i[[i]] <- m - n_na_items_i
        mv_i[[i]] <- mv[-na_items_i]
        oj_vec_i <- oj_vec[-na_items_i]
        oj_max_i[[i]] <- oj_max[-na_items_i]
      }

      ## calculate category totals and person totals for NA-group i
      ctot_i[[i]] <- vector("list", length = m_i[[i]])
      for (j in seq_len(m_i[[i]])) ctot_i[[i]][[j]] <- as.vector(tapply(weights_i[[i]], factor(y_i[, j], levels = oj_vec_i[[j]]), sum))
      ctot_i[[i]] <- unlist(lapply(ctot_i[[i]], "[", -1))
      ctot_i[[i]][is.na(ctot_i[[i]])] <- 0
      rs_i[[i]] <- rowSums(y_i)
      ptot_i[[i]] <- as.vector(tapply(weights_i[[i]], factor(rs_i[[i]], levels = 0:sum(oj_max_i[[i]])), sum))
      ptot_i[[i]][is.na(ptot_i[[i]])] <- 0

      ## calculate helper variables for gradient
      parindex_i <- unlist(lapply(oj_vec_i, "[", -1))
      itemindex_i <- unlist(rep.int(mv_i[[i]], oj[mv_i[[i]]]))
      est_par_i[[i]] <- !unlist(lapply(nl_cats, "[", -1)[mv_i[[i]]])
      xmat_i[[i]] <- matrix(FALSE, nrow = length(na_obs_i[[i]]), ncol = sum(oj[mv_i[[i]]]))
      for (j in 1:length(na_obs_i[[i]])) xmat_i[[i]][j, ] <- (y[na_obs_i[[i]], , drop = FALSE])[j, itemindex_i] == parindex_i
    }

    ## fun for mapply (calculates cll contributions per na_group with variable par's)
    cll_i <- function (ctot_i, par_i, ptot_i, esf_par_i) {
      - sum(ctot_i * par_i) - sum(ptot_i * log(elementary_symmetric_functions(par = esf_par_i, order = 0, diff = diff)[[1]]))
    }
    
    ## objective function: conditional log-likelihood (build incrementally by looping over the NA-patterns
    cloglik <- function (par) {

      ## add first zero par, then select current parameters
      ## include epsilons for esf when there are nullcats & nullcats == "keep" (see Wilson, 1993)
      if (any_nl_cat_item && nullcats == "keep") {
        all_par[est_par] <- par
        parx <- all_par
      } else {
        parx <- c(0, par)
      }
      esf_par <- lapply(par_i, function (x) parx[x])
      par <- lapply(esf_par, function (x) x[!is.na(x)])
      esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))

      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, ctot_i, par, ptot_i, esf_par, SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## static helper variables for analytical gradient
    grad <- matrix(0, nrow = n, ncol = sum(oj))
    itemindex <- rep.int(mv, oj)

    ## analytical gradient
    agrad <- function (par) {
        
        ## elementary symmetric functions
        if (any_nl_cat_item && nullcats == "keep") {
            all_par[est_par] <- par
            esf_par <- all_par
        } else {
            esf_par <- c(0, par)
        }
        esf_par <- lapply(par_i, function (x) esf_par[x])
        esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))
        esf <- mapply(elementary_symmetric_functions, par = esf_par, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)

        ## loop over observed NA patterns and calculate gradients
        for(i in seq_along(lev_na_patterns)) {

            ## select gamma zero and first derivatives with rs_i
            gamma0_i <- esf[[i]][[1]][rs_i[[i]] + 1]
            gamma1_i <- apply(esf[[i]][[2]], 2, "[", rs_i[[i]] + 1)
            if (!is.matrix(gamma1_i)) gamma1_i <- matrix(gamma1_i, nrow = 1)

            ## finally: the gradient for NA group i
            if (any_nl_cat_item && nullcats == "keep") gamma1_i <- gamma1_i[, est_par_i[[i]], drop = FALSE]
            grad[na_obs_i[[i]], itemindex %in% mv_i[[i]]] <- weights_i[[i]] * (- xmat_i[[i]] + (gamma1_i/ gamma0_i))
        }

        return(- colSums(grad[, -1, drop = FALSE]))
    }
  }

  ## optimization
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = hessian, control = list(reltol = reltol, maxit = maxit, ...))

  ## final estimates ...
  est <- opt$par
  names(est) <- if (is.null(colnames(y))) {
    paste("I", rep.int(which(ident_items), oj), "-C", unlist(lapply(oj, seq_len)), sep = "")[-1] 
  } else paste(rep(colnames(y)[which(ident_items)], oj), "-C", unlist(lapply(oj, seq_len)), sep = "")[-1]
  
  ## ... and (if requested) esf of these ...
  if (full) {
    if (any_nl_cat_item && nullcats == "keep") {
      all_par[est_par] <- est
      parx <- all_par
    } else {
      parx <- c(0, est)
    }
    esf <- if (any_na) {
      esf_par <- lapply(par_i, function (x) parx[x])
      esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))
      mapply(elementary_symmetric_functions, par = esf_par, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)
    } else {
      elementary_symmetric_functions(par = split(parx, rep.int(mv, oj_max)), order = 1, diff = diff)
    }

    ## ... as well as variance-covariance matrix
    if (hessian) {
      vc <- opt$hessian
      vc <- solve(vc)
    } else {
      vc <- matrix(NA, nrow = length(est), ncol = length(est))
    }
    rownames(vc) <- colnames(vc) <- names(est)
  } else {
    esf <- NULL
    vc <- NULL
  }

  ## collect results, set class, and return
  res <- list(coefficients = est,
              vcov = vc,
              data = y,
              items = ident_items,
              categories = lapply(oj_vec, "[", -1),
              n = sum(weights_org > 0),
              n_org = n_org,
              weights = if (identical(as.vector(weights_org), rep.int(1L, n_org))) NULL else weights_org,
              na = any_na,
              nullcats = if (any_nl_cat_item && nullcats == "keep") lapply(nl_cats, "[", -1) else NULL,
              esf = esf,
              loglik = -opt$value,
              df = length(est),
              code = opt$convergence,
              iterations = tail(na.omit(opt$counts), 1L),
              reltol = reltol)
  class(res) <- "PCModel"
  return(res )
}

print.PCModel <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("PCM item category parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

coef.PCModel <- function (object, ...) object$coefficients

vcov.PCModel <- function (object, ...) object$vcov

logLik.PCModel <- function (object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.PCModel <- function (object, ...) if (is.null(object$weights)) rep.int(1L, object$n_org) else object$weights

summary.PCModel <- function (object, vcov. = NULL, ...) {
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if (is.null(vcov.)) 
    vc <- vcov(object)
  else {
    if (is.function(vcov.)) vc <- vcov.(object)
    else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.PCModel"
  return(object)
}

print.summary.PCModel <- function (x, digits = max(3, getOption("digits") - 3),
                                   signif.stars = getOption("show.signif.stars"), ...)
{
  if (is.null(x$call)) {
    cat("\nPartial credit model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if (any(!x$items)) cat("Excluded items:",
                         paste(names(x$items)[!x$items], collapse = ", "), "\n\n")

  cat("Item category parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
      "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

### plot.PCModel() draws a row of effect plots, which show for each item the most probable categories over the theta axis.
plot.PCModel <-  function(x, type = NULL, names = NULL, ref = NULL, main = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, ...)
{
  ## setup par, number of items, axis labels
  par(mar = c(3, 4, 1, 2.5))
  items <- x$items
  m <- length(items)
  if (is.null(names)) names  <- paste("Item", 1:m) else stopifnot(length(names) == m)

  ## fetch request threshold parameters
  tr <- threshold.PCModel(x, type = type, ref = ref, simplify = FALSE)
  
  ## setup axis range and positions
  if (is.null(ylim)) ylim <- extendrange(unlist(tr, use.names = FALSE), f = 0.25)
  xi <- 0:m + c(0:(m-1), m-1) * off
  xlim <- c(xi[1], xi[m+1])
  xi_id <- xi[items]

  ## setup graphical window
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", xaxs = "i", yaxs = "i", axes = FALSE, ylab = ylab)
  
  ## plot items 
  for (j in seq_along(tr)) {
    rect(xleft = xi_id[j], xright = xi_id[j] + 1, ybottom = c(ylim[1], tr[[j]]), ytop = c(tr[[j]], ylim[2]), col = gray.colors(length(tr[[j]])+1), ...)
  }
  
  ## indicate unordered parameters
  if ((is.null(type) || type == "mode")) {
    orgtr <- threshold.PCModel(x, type = "unmodified", ref = ref, simplify = FALSE)
    uo_items <- which(!sapply(mapply(all.equal, tr, orgtr, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
    for (j in uo_items) {
      uo_pars <- setdiff(orgtr[[j]], tr[[j]])
      lines(x = rep(c(xi_id[j], xi_id[j] + 1, NA), length(uo_pars)), y = rep(uo_pars, each = 3), col = "red", lty = 2)
    }
  }

  ## add axes
  axis(2, ...)
  axis(4, ...)
  axis(1, at=(xi[-(m+1)] + 0.5), labels = names, ...)
  box()
}

### plotCCC.PCModel() draws the category characteristic curves of a number of specified items for given model object
plotCCC.PCModel <- function (object, items = NULL, sum0 = TRUE, lyt = NULL, xlim = NULL, ylim = c(0, 1), col = NULL, n = 101,
                             main = NULL, xlab = "Latent Trait", ylab = "Probability", ...)
{
  ## setup number of items/categories, threshold parameters
  if (is.null(items)) items <- seq_len(sum(object$items)) else stopifnot(is.numeric(items) && all(items %in% which(object$items)))
  m <- length(items)
  mvec <- 1:m
  if (inherits(object, "RSModel")) {
    thresh <- coef(object, type = "PCM", sum0 = sum0, as.list = TRUE)[items]
    ojvec <- object$categories[items]  + 1
  } else {
    thresh <- coef(object, type = "threshold", sum0 = sum0, as.list = TRUE)[items]
    ojvec <- (sapply(object$categories, length) + 1)[items]
  }

  ## setup plotting parameters
  if (is.null(lyt)) lyt <- matrix(mvec, ncol = 1, nrow = m)
  if (is.null(xlim)) xlim <- extendrange(unlist(thresh), f = 2)
  if (is.null(col)) col <- rainbow(max(ojvec))
  theta <- seq(from = xlim[1], to = xlim[2], length.out = n)
  probs <- ppcm(theta = theta, delta = thresh)

  ## setup plotting area, setup par
  layout(lyt)
  par(mar = c(4, 4, 2.5, 0.5))

  ## loop through items and plot CCC
  for (j in mvec) {
    oj <-  1:ojvec[j]
    plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, main = if(is.null(main)) paste("Category Characteristic Curves Item", items[j]) else main, ylab = ylab, xaxs = "i")
    for (k in oj) lines(x = theta, y = probs[[j]][, k], col = col[k], ...)
    legend("topright", legend = paste("Cat.", oj), lty = 1, col = col[oj], horiz = TRUE, bty = "n")
  }
}

itempar.PCModel <- function (object, ref = NULL, vcov = TRUE, simplify = TRUE, ...) {
  
  ## extract cf and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  m <- length(cf)
  lbs <- names(cf)
  lbs[1] <- names(cf)[1] <- if (lbs[2] == "I1-C2") "I1-C1" else if (lbs[2] == "I2-C1") "I1-C1" else paste0(colnames(object$data)[1], "-C1")
  
  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with threshold parameter labels or a numeric vector with threshold parameter indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)

  ## impose restriction
  cf <- D %*% cf
  names(cf) <- lbs

  ## create adjusted vcov if requested
  if (vcov) {
    vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
    colnames(vc) <- rownames(vc) <- lbs
  }

  ## create list if requested
  if (!simplify) {
    cf <- relist(cf, object$categories)
    lbs <- relist(lbs, object$categories)
    cf <- mapply("names<-", cf, lbs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }

  ## return results
  rval <- structure(cf, class = "itempar", model = "PCM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

threshold.PCModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, simplify = TRUE, ...) {
  ## check type, extract item parameters, backup attributes
  type <- match.arg(type)
  ip <- itempar.PCModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  mdl <- attr(ip, "model")
  ref <- attr(ip, "ref")
  ip <- lapply(ip, function (j) diff.default(c(0, j))) # item category parameters -> item threshold parameters

  ## calculate threshold parameters (if necessary)
  if (type == "mode") {
    ## check if all threshold parameters are in order, if not, calculate sorted ones
    us <- sapply(ip, is.unsorted)
    if (any(us)) {
      usj <- which(us)
      for (j in usj) {
        ipj <- ip[[j]]
        nj <- length(ipj)
        
        ## check if there is a point with a biggest parameter, if yes, take mean
        for (i in 1:nj) {
          if (all(ipj[i] > ipj[(i+1):nj])) {
            ipj[i] <- mean(ipj[i:nj])
            ipj <- ipj[-(i+1:nj)]
            break
          }
        }
        
        ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
        while(is.unsorted(ipj)) {
          uo_pos <- which(diff(ipj) < 0)                             # locate unordered parameters, returns position of the first
          ipj[uo_pos] <- (ipj[uo_pos] + ipj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
          ipj <- ipj[-(uo_pos + 1)]                              # remove second
        }

        ip[[j]] <- ipj
      }
    }
  } else {
    ## backup labels
    lbs <- lapply(ip, names)
    oj <- sapply(ip, length)

    if (type == "median") {
      ## function to find locations on theta axis
      zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) {
        rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5
      }
      
      ## loop though items and find locations by means of zmedian() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(1:oj[j], function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = ip[[j]], geq = geq, ncat = oj[j] + 1)$root)
      }
    } else if (type == "mean") {
      ## function to find locations on theta axis
      xpct <- lapply(oj, function (oj) 1:oj - 0.5)
      zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
      
      ## loop though items and find locations by means of zexpct() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(xpct[[j]], function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = ip[[j]], expct = xp)$root)
      }
    }

    ## add labels again
    ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  }

  ## backup attributes, simplify if requested, then set attributes again
  if (simplify) ip <- unlist(ip)
  return(structure(ip, "class" = "threshold", "model" = mdl, "ref" = ref, "type" = type))
}

personpar.PCModel <- function (object, ref = NULL, vcov = TRUE, start = NULL, tol = 1e-6, ...) {

  ## extract item parameters and relevant informations #FIXME: constrained parameter included?
  ipar <- itempar.PCModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  m <- length(ipar)
  oj <- sapply(ipar, length)
  ojvl <- lapply(oj, seq)
  ojv <- unlist(ojvl)
  rng <- 1:(sum(oj) - 1)

  ## if vcov is requested, fetch relevant data for loglik
  if (vcov) {
    y <- object$data[weights(object) > 0, , drop = FALSE]
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = sum(oj) - 1)
    cs <- sapply(1:m, function (j) tabulate(y[, j], nbins = oj[j]))

    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    rf <- rf[rf != 0]
  }

  ## start values / range of person parameters
  if(is.null(start)) start <- if (vcov) numeric(length(rng)) else c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?

  if (vcov) {
    ## objective function
    cloglik <- function (ppar) {
      pparx <- lapply(ojvl, function (x) outer(x, ppar)) # l * theta_i
      - sum(rf * rng * ppar) + sum(unlist(mapply("*", cs, ipar))) + sum(rf * rowSums(log(1 + mapply(function (x, y) colSums(exp(x - y)), pparx, ipar))))
    }

    ## optimization & target variables
    opt <- nlm(cloglik, start, gradtol = tol, hessian = TRUE, check.analyticals = FALSE, ...)
    ppar <- opt$estimate
    vc <- solve(opt$hessian)
  } else {
    ## iterate over raw scores (from 1 until rmax-1)
    ppar <- sapply(rng, function(rawscore) {
      uniroot(function(ppar) rawscore - sum(ojv * exp(ojv * ppar - unlist(ipar)) / unlist(sapply(1:m, function (j) rep.int(1 + sum(exp(ojvl[[j]] * ppar - ipar[[j]])), oj[j])))), interval = start, tol = tol)$root
    })
  }

  ## return named vector of thetas
  return(structure(ppar, .Names = rng, class = "personpar", model = "PCM", vcov = if(vcov) vc else NULL))
}

predict.PCModel <- function (object, newdata = NULL, type = c("probability", "cumprobability", "mode", "median", "mean"),
                             ref = NULL, ...) {
  
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < ncol(object$data)]
    newdata <- personpar.PCModel(object, ref = ref, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }
  nms <- names(newdata)
  ## itempars, labels and raw probabilities
  ip <- itempar.PCModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  lbs <- lapply(ip, names)
  lbs <- unlist(lapply(lbs, function (j) c(gsub("(.*)-C1", "\\1-C0", j[1]), j)))
  probs <- ppcm(theta = newdata, delta = ip)

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = {
           ## create (named) matrix with probabilities
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "cumprobability" = {
           ## create (named) matrix with probabilities
           probs <- lapply(probs, function (probsj) t(apply(probsj, 1, function (x) rev(cumsum(rev(x))))))
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, gsub("(.*)-(.*)", "\\1>=\\2", lbs))
           probs
           },
         "mode" = {
           ## create (named) matrix with probabilities
           probs <- sapply(probs, apply, 1, which.max) - 1
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         },
         "median" = {
           ## create (named) matrix with probabilities
           probs <- sapply(probs, function (probsj) apply(probsj, 1, which.max) - 1)
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         },
         "mean" = {
           ## create (named) matrix with probabilities
           probs <- round(sapply(probs, function (probsj) apply(probsj, 1, function (j) sum(j * 0:(length(j) - 1)))))
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         })
}



### misc. internal functions

## ppcm: calculate response probabilities for given thetas and deltas under the PCM.
ppcm <- function(theta = NULL, delta = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(delta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, ppcm, delta = delta))
  if (is.list(delta)) return(lapply(delta, ppcm, theta = theta))

  ## calculate probabilities
  num <- cbind(0, outer(theta, delta, "-")) # all possible differences, 0 for category zero (\sum_0^0 \def 0)
  num <- t(exp(apply(num, 1, cumsum)))      # numerator: all possible cumulative sums
  denom <- rowSums(num)                     # denominator: sum over cumulative sums
  return(num/denom)
}

## rpcm: calculate response matrices for given thetas and deltas under the PCM.
rpcm <- function(theta = NULL, delta = NULL, nullcats = FALSE, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(delta))
  
  ## if list input, recurse... (for deltas: one list means several items, two means several groups of items)
  if (is.list(theta)) return(lapply(theta, rpcm, delta = delta, nullcats = nullcats, return_setting = return_setting))
  if (is.list(delta) && is.list(delta[[1]])) return(lapply(delta, rpcm, theta = theta, nullcats = nullcats, return_setting = return_setting))

  ## calculate response probabilities
  probs <- ppcm(theta = theta, delta = delta)
  if(!is.list(probs)) probs <- list(probs)

  ## calculate response matrices for given set of item
  rsp_item_i <- function (probmat_i) {  #inline function to calculate responses per item
    oj <- ncol(probmat_i)               #calculate rsp vector once
    rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
    if(!nullcats) {                     #recalculate if null categories not allowed
      while (length(unique.default(rsp)) != oj) {
        rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
      }
    }
    return(rsp-1)
  }
  res <- lapply(probs, rsp_item_i)
  res <- do.call(cbind, res)

  if (return_setting)
    return(list(delta = delta, theta = theta, data = res))
  else
    return(res)
}


### generic functions
plotCCC <- function (object, ...) UseMethod("plotCCC")
