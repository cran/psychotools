## estimate a gpcm in slope / intercept formulation using mirt (MML & EM)

gpcmodel <- function(y, weights = NULL, impact = NULL,
  type = c("GPCM", "PCM"), grouppars = FALSE, vcov = TRUE, nullcats = "downcode",
  start = NULL, method = "BFGS", maxit = 500L, reltol = 1e-5, ...)
{

  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for fitting a gpcmodel.", call. = FALSE)
  }

  ## handle arguments and defaults
  type <- match.arg(type, c("GPCM", "PCM"))

  if(nullcats != "downcode") stop("Currently only downcoding is supported for nullcats.")

  if(is.logical(vcov)) vcov <- c("none", "Oakes")[vcov + 1L]
  vcov <- match.arg(vcov, c("Oakes", "Richardson", "forward", "central", "crossprod",
    "Louis", "sandwich", "sandwich.Louis", "SEM", "Fisher", "complete", "none"))

  ## response matrix
  y <- as.matrix(y)
  M_org <- nrow(y)
  N <- ncol(y)
  if(is.null(colnames(y))) colnames(y) <- paste0("Item", gsub(" ", "0", format(1L:N)))
  any_na <- any(is.na(y))

  ## process weights
  if(is.null(weights)) weights <- rep.int(1L, M_org) else stopifnot(length(weights) == M_org)
  y <- y[weights > 0, , drop = FALSE]
  impact_org <- impact
  impact <- impact[weights > 0]
  weights_org <- weights
  weights <- weights[weights > 0]
  M <- nrow(y)

  ## impact (if any) needs to be a factor of suitable length
  ## with sensible reference category (enforced to be the largest group)
  if(!is.null(impact)) {
    if(!is.factor(impact)) impact <- as.factor(impact)
    stopifnot(M == length(impact))

    ## reorder the impact variable in decreasing order so the largest group
    ## is the reference group --> estimation stability, faster convergence
    ix <- which.max(table(impact))
    impact <- relevel(impact, ix)
    lvls <- levels(impact)
    nlvls <- length(lvls)

    ## multiple group models are estimated using nlminb
    method <- "nlminb"
  }

  ## move categories to zero if necessary, get number of categories per item
  mincat <- apply(y, 2L, min, na.rm = TRUE)
  if(any(mincat > 0)) {
    warning("Minimum score is not zero for all items (", paste(which(mincat > 0), collapse = ", "), "). These items are scaled to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }
  oj <- apply(y, 2L, max, na.rm = TRUE)
  oj_vec <- lapply(oj, seq, from = 0L)

  ## check for null categories and warn...
  nl_cats <- lapply(1L:N, function (j) !(oj_vec[[j]] %in% y[, j]))
  nl_cats_items <- which(sapply(nl_cats, any))
  any_nl_cat_item <- length(nl_cats_items) > 0
  if(any_nl_cat_item) {
    warning("There are items with null categories (", paste(nl_cats_items, sep = "", collapse = ", "), ").")
    ## ... then treat according to nullcats, currently mirt only allows for "downcode"
    if(nullcats == "downcode") {
      oj <- oj - sapply(nl_cats, sum)
      oj_vec <- lapply(1L:N, function (j) oj_vec[[j]][1L:(oj[j] + 1)])
      for(j in nl_cats_items) {
        missing_scores <- sort(which(nl_cats[[j]]), decreasing = TRUE) - 1L
        for(score in missing_scores) y[!is.na(y[, j]) & y[, j] > score, j] <-  y[!is.na(y[, j]) & y[, j] > score, j] - 1L
      }
    }
  }
  names(oj_vec) <- colnames(y)

  ## call mirt, mostly defaults
  ## if impact --> multiple group model, reference group mean 0 and var 1,
  ## other groups free, item parameters restricted across whole sample
  ## otherwise --> simple mirt model, group mean 0 and var 1
  model <-
  if(!is.null(impact)) {
    pars <- mirt::multipleGroup(data = y, model = 1L, group = impact,
      itemtype = "gpcm", invariance = c("free_means", "free_var", colnames(y)),
      pars = "values")

    constrain <- if(type == "PCM") list(which(pars$name == "a1")) else NULL

    if(!is.null(start)) {
      if(type == "PCM") {
        ## start should be a vector of item parameters ordered itemwise
        ## and means and vars appended groupwise
        ## a list can be provided containing a single slope,
        ## a list of intercept vectors itemwise, a vector of means and
        ## a vector of vars
        if(is.list(start)) {
          n <- length(start)
          start <- c(start[[1L]], unlist(start[[2L]]),
            as.vector(do.call(rbind, start[(n - 1L):n])))
        }
        if(length(start) != sum(pars$est[pars$group == lvls[1L]]) +
          (2L * (nlvls - 1L)) - (N - 1L)) {
          stop("Argument 'start' is misspecified (see ?gpcmodel for possible values).")
        }
        start_i <- start[1L:(length(start) - (2L * (nlvls - 1L)))]
        start_g <- start[-seq_along(start_i)]
        ## if the impact variable has been reordered, adjust the starting values
        ## of the means and variances accordingly
        if(ix != 1L) {
          mx <- seq(from = 1L, to = (2L * (nlvls - 1L)) - 1L, by = 2L)
          means <- c(0, start_g[mx])
          means <- (means - means[ix])[-nlvls]
          vx <- seq(from = 2L, to = (2L * (nlvls - 1L)), by = 2L)
          vars <- c(1, start_g[vx])
          vars <- (vars / vars[ix])[-nlvls]
          start_g <- as.vector(rbind(means, vars))
        }
        pars$value[pars$est & pars$name == "a1"] <- start_i[1L]
        pars$value[pars$est & pars$name != "a1" & pars$item != "GROUP"] <-
          rep.int(start_i[-1L], nlvls)
        pars$value[pars$est & pars$item == "GROUP"] <- start_g
      } else {
        ## start should be a vector of item parameters ordered itemwise
        ## and means and vars appended groupwise
        ## a list can be provided containing a vector of slopes,
        ## a list of intercept vectors itemwise, a vector of means and
        ## a vector of vars
        if(is.list(start)) {
          n <- length(start)
          start <- c(unlist(mapply(c, start[[1L]], start[[2L]], SIMPLIFY = FALSE)),
            as.vector(do.call(rbind, start[(n - 1L):n])))
        }
        if(length(start) != sum(pars$est[pars$group == lvls[1L]]) +
          (2L * (nlvls - 1L))) {
          stop("Argument 'start' is misspecified (see ?gpcmodel for possible values).")
        }
        start_i <- start[1L:(length(start) - (2L * (nlvls - 1L)))]
        start_g <- start[-seq_along(start_i)]
        ## if the impact variable has been reordered, adjust the starting values
        ## of the means and variances accordingly
        if(ix != 1L) {
          mx <- seq(from = 1L, to = (2L * (nlvls - 1L)) - 1L, by = 2L)
          means <- c(0, start_g[mx])
          means <- (means - means[ix])[-nlvls]
          vx <- seq(from = 2L, to = (2L * (nlvls - 1L)), by = 2L)
          vars <- c(1, start_g[vx])
          vars <- (vars / vars[ix])[-nlvls]
          start_g <- as.vector(rbind(means, vars))
        }
        pars$value[pars$est & pars$item != "GROUP"] <- rep.int(start_i, nlvls)
        pars$value[pars$est & pars$item == "GROUP"] <- start_g
      }
    }
    mirt::multipleGroup(data = y, model = 1L, group = impact,
      itemtype = "gpcm", invariance = c("free_means", "free_var", colnames(y)),
      SE = (vcov != "none"), SE.type = if(vcov == "none") "Oakes" else vcov,
      optimizer = method, calcNull = FALSE, TOL = reltol,
      technical = list(NCYCLES = maxit), verbose = FALSE, pars = pars,
      constrain = constrain, survey.weights = weights, ...)
  } else {
    grouppars <- FALSE
    pars <- mirt::mirt(data = y, model = 1L, itemtype = "gpcm", pars = "values")

    constrain <- if(type == "PCM") list(which(pars$name == "a1")) else NULL

    if(!is.null(start)) {
      if(type == "PCM") {
        ## start should be a vector of item parameters ordered itemwise
        ## a list can be provided containing a single slope and
        ## a list of intercept vectors itemwise
        if(is.list(start)) {
          start <- c(start[[1L]], unlist(start[[2L]]))
        }
        if(length(start) != sum(pars$est) - (N - 1L)) {
          stop("Argument 'start' is misspecified (see ?gpcmodel for possible values).")
        }
        pars$value[pars$est & pars$name == "a1"] <- start[1L]
        pars$value[pars$est & pars$name != "a1"] <- start[-1L]
      } else {
        ## start should be a vector of item parameters ordered itemwise
        ## a list can be provided containing a vector of slopes and
        ## a list of intercept vectors itemwise
        if(is.list(start)) {
          start <- unlist(mapply(c, start[[1L]], start[[2L]], SIMPLIFY = FALSE))
        }
        if(length(start) != sum(pars$est)) {
          stop("Argument 'start' is misspecified (see ?gpcmodel for possible values).")
        }
        pars$value[pars$est] <- start
      }
    }
    mirt::mirt(data = y, model = 1L, itemtype = "gpcm", SE = (vcov != "none"),
      SE.type = if(vcov == "none") "Oakes" else vcov, optimizer = method,
      calcNull = FALSE, TOL = reltol, technical = list(NCYCLES = maxit),
      verbose = FALSE, pars = pars, constrain = constrain, survey.weights = weights, ...)
  }

  ## check for negative slopes
  if(is.null(impact)) {
    a <- sapply(mirt::coef(model)[-(N + 1L)], function(x) {
      x[1L, ][grep("a1$", colnames(x))]
    })
  } else {
    a <- sapply(mirt::coef(model)[[1L]][-(N + 1L)], function(x) {
      x[1L, ][grep("a1$", colnames(x))]
    })
  }
  neg_a <- which(unlist(a) < 0)
  if(length(neg_a > 0)) {
    warning("Negative slopes were estimated (", paste(neg_a, sep = "", collapse = ", "), ").")
  }

  ## coefficients, names and useful info
  epars <- mirt::mod2values(model)
  if(!is.null(impact)) {
    epars_i <- epars[epars$group == levels(impact)[1L] & epars$est, ]
    epars_g <- epars[epars$group != levels(impact)[1L] & epars$item == "GROUP" &
      epars$est, ]
    est <- c(epars_i$value, epars_g$value)
    names(est) <- c(paste(epars_i$item, epars_i$name, sep = "-"),
      paste(epars_g$group, epars_g$name, sep = "-"))
  } else {
    est <- epars$value[epars$est]
    names(est) <- paste(epars$item[epars$est], epars$name[epars$est], sep = "-")
  }

  ## if grouppars = FALSE and multiple group model drop the group parameters
  if(!grouppars && !is.null(impact)) {
    estnms <- names(est)
    estnms <- estnms[-c(grep("MEAN", estnms), grep("COV", estnms))]
    est <- est[which(names(est) %in% estnms)]
  }

  ## handle type of model restriction
  if(type == "PCM") {
    est <- est[-grep("-a1$", names(est))[-1L]]
    names(est)[1L] <- "a1"
  }
  nest <- length(est)

  ## estimated variance-covariance matrix (if any)
  vc <- if(vcov != "none" && length(mirt::vcov(model)) != 1) {
    mirt::extract.mirt(model, "vcov")[1L:nest, 1L:nest]
  } else {
    matrix(NA_real_, nest, nest)
  }
  rownames(vc) <- colnames(vc) <- names(est)

  ## convergence code
  code <- if(mirt::extract.mirt(model, "converged")) {
    0L
  } else if(mirt::extract.mirt(model, "iterations") == maxit) {
    1L
  } else {
    1L  ## FIXME: it would be nice to have the optim convergence code here...
        ## but mirt currently does not provide this, hence resort to code 1
  }

  ## collect return values
  res <- list(
    coefficients = est,
    vcov = structure(vc, method = vcov),
    data = y,
    items = rep.int(TRUE, N),
    categories = lapply(oj_vec, "[", -1L),
    n = sum(weights_org > 0),
    n_org = M_org,
    weights = if(identical(as.vector(weights_org), rep.int(1L, M_org))) NULL else weights_org,
    na = any_na,
    nullcats = NULL,
    impact = impact,
    loglik = mirt::extract.mirt(model, "logLik"),
    df = mirt::extract.mirt(model, "nest"),
    code = code,
    iterations = structure(mirt::extract.mirt(model, "iterations"), names = "EM cycles"),
    reltol = reltol,
    method = method,
    grouppars = grouppars,
    type = type,
    mirt = model
  )
  class(res) <- "gpcmodel"
  return(res)
}



print.gpcmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if(x$grouppars) {
    cat(paste(x$type, "slopes, intercepts, and estimated group parameters:\n"))
  } else {
    cat(paste(x$type, "slopes and intercepts:\n"))
  }
  print(coef(x), digits = digits, ...)
  invisible(x)
}



coef.gpcmodel <- function(object, ...)
{
  object$coefficients
}



vcov.gpcmodel <- function(object, ...)
{
  object$vcov
}



logLik.gpcmodel <- function(object, ...)
{
  structure(object$loglik, df = object$df, class = "logLik")
}



weights.gpcmodel <- function(object, ...)
{
  if(is.null(object$weights)) rep.int(1L, object$n_org) else object$weights
}



summary.gpcmodel <- function(object, vcov. = NULL, ...)
{
  cf <- coef(object)
  vc <-
  if(is.null(vcov.)) {
    vcov(object)
  } else {
    if(is.function(vcov.)) {
      vcov.(object)
    } else {
      vcov.
    }
  }
  cf <- cbind(cf, sqrt(diag(vc)), cf / sqrt(diag(vc)),
    2 * pnorm(-abs(cf / sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf
  class(object) <- "summary.gpcmodel"
  return(object)
}



print.summary.gpcmodel <- function(x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat(sprintf("\n%s\n\n",
      if(x$type == "PCM") "Partial credit model" else "Generalized partial credit model"))
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")
  }
  if(x$grouppars) {
    cat("Slopes, intercepts and estimated group parameters:\n")
  } else {
    cat("Slopes and intercepts:\n")
  }
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
    na.print = "NA", ...)
  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n")
  cat(sprintf("Number of EM cycles: %s\n", x$iterations))
  cat(sprintf("M-step optimizer: %s\n", x$method))
  invisible(x)
}



plot.gpcmodel <- function(x, type = c("regions", "profile", "curves", "information", "piplot"), ...)
{
  type <- match.arg(type, c("regions", "profile", "curves", "information", "piplot"))
  switch(type, curves = curveplot(x, ...),
    regions = regionplot(x, ...), profile = profileplot(x, ...),
    information = infoplot(x, ...),
    piplot = piplot(x, ...))
}



predict.gpcmodel <- function(object, newdata = NULL,
  type = c("probability", "cumprobability", "mode", "median", "mean",
  "category-information", "item-information", "test-information"), ref = NULL, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for gpcmodels.", call. = FALSE)
  }

  if(!is.null(ref)) {
    stop("'ref' cannot be applied for models estimated via MML.")
  }
  type <- match.arg(type, c("probability", "cumprobability", "mode", "median",
    "mean", "category-information", "item-information", "test-information"))
  if(is.null(newdata)) {
    newdata <- personpar(object, personwise = TRUE)
    names(newdata) <- NULL
  }
  N <- sum(object$items)
  nms <- names(newdata)
  a <- as.list(discrpar(object, vcov = FALSE))
  b <- threshpar(object, type = "mode", vcov = FALSE)
  lbs <- unlist(mapply(paste, names(b), lapply(b, function(j) {
    c("C0", names(j))
  }), sep = "-", SIMPLIFY = FALSE), use.names = FALSE)
  probs <- pgpcm(theta = newdata, a = a, b = b)
  ## call mirt
  if(grepl("information", type)) {
    clnms <- names(b)
    if(type == "category-information") {
      info <- matrix(NA_real_, length(newdata), length(lbs))
      colnames(info) <- lbs
    } else {
      info <- matrix(NA_real_, length(newdata), N)
      colnames(info) <- clnms
    }
    for(j in 1L:N) {
      idx <- grepl(clnms[j], lbs)
      iteminfo <- mirt::iteminfo(mirt::extract.item(object$mirt, j, 1L),
        newdata)
      if(type == "category-information") {
        info[, idx] <- probs[[j]] * iteminfo
      } else {
        info[, j] <- iteminfo
      }
    }
  }
  switch(type,
    "probability" = {
      probs <- do.call("cbind", probs)
      dimnames(probs) <- list(nms, lbs)
      probs
    },
    "cumprobability" = {
      probs <- lapply(probs, function(probsj) {
        t(apply(probsj, 1L, function(x) {
          rev(cumsum(rev(x)))
        }))
      })
      probs <- do.call("cbind", probs)
      dimnames(probs) <- list(nms, gsub("(.*)-(.*)", "\\1>=\\2", lbs))
      probs
    },
    "mode" = {
      probs <- sapply(probs, apply, 1L, which.max) - 1L
      dimnames(probs) <- list(nms,
        unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
      probs
    },
    "median" = {
      probs <- sapply(probs, function(probsj) {
        apply(probsj, 1L, which.max) - 1L
      })
      dimnames(probs) <- list(nms,
        unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
      probs
    },
    "mean" = {
      probs <- round(sapply(probs, function(probsj) {
        apply(probsj, 1L, function(j) {
          sum(j * 0L:(length(j) - 1L))
        })
      }))
      dimnames(probs) <- list(nms,
        unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
      probs
    },
    "category-information" = info,
    "item-information" = info,
    "test-information" = matrix(rowSums(info), ncol = 1L))
}



itempar.gpcmodel <- function(object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  cf <- unlist(threshpar(object))
  N <- sum(object$items)
  B <- length(cf)
  oj <- sapply(object$categories, length)
  ojc <- cumsum(oj)
  ojc_l <- vector("list", N)
  ojc_t <- c(0L, ojc)
  for(l in 2L:length(ojc_t)) {
    ojc_l[[l - 1L]] <- (ojc_t[l - 1L] + 1L):ojc_t[l]
  }
  lbs <- colnames(object$data)
  if(is.null(ref)) {
    ref <- NULL
    ## scale defined via means and vars
  } else if(is.vector(ref) & is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if(is.vector(ref) & is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1L:N))
  } else if(is.matrix(ref) & is.numeric(ref)) {
    stopifnot(nrow(ref) == N & ncol(ref) == N)
  } else {
    stop("Argument 'ref' is misspecified (see ?itempar for possible values).")
  }
  ## D = ref, C = transformation
  if(is.matrix(ref)) {
    D <- ref
  } else {
    D <- diag(N)
    D[, ref] <- D[, ref] - 1 / length(ref)
  }
  C <- matrix(0, N, B)
  for(j in 1L:N) C[j, ojc_l[[j]]] <- 1 / oj[j]
  ## apply C
  cf <- as.vector(C %*% cf)
  ## delta method C
  if(vcov) {
    vc <- C %*% vcov(threshpar(object)) %*% t(C)
  }
  ## apply D
  cf <- as.vector(D %*% cf)
  ## delta method D
  vc <-
  if(vcov) {
    D %*% vc %*% t(D)
  } else {
    matrix(NA_real_, N, N)
  }
  names(cf) <- rownames(vc) <- colnames(vc) <- lbs
  if(!alias) {
    if(is.matrix(ref)) {
    ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
      stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
    } else if(!is.null(ref)) {
      cf <- cf[-ref[1L]]
      vc <- vc[-ref[1L], -ref[1L]]
      alias <- paste0("I", ref[1L])
      names(alias) <- lbs[ref[1L]]
    }
  }
  rv <- structure(cf, class = "itempar", model = "GPCM", ref = ref,
    alias = alias, vcov = vc)
  return(rv)
}



threshpar.gpcmodel <- function(object, type = c("mode", "median", "mean"),
  ref = NULL, alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = TRUE, ...)
{
  type <- match.arg(type, c("mode", "median", "mean"))
  N <- sum(object$items)
  oj <- sapply(object$categories, length)
  ojc <- cumsum(oj)
  K <- max(ojc)
  coefs <- coef(object)
  estnms <- names(coefs)
  ## go from slope / intercept to IRT
  a <- mapply(rep, as.list(discrpar(object, vcov = FALSE)), oj + 1L, SIMPLIFY = FALSE)
  d <- coefs[grep("-d[[:digit:]]$", estnms)]
  d <- lapply(split(d, rep(1L:N, oj)), function(x) c(0, x))
  d_t <- mapply(function(x, y) - (x / y), d, a, SIMPLIFY = FALSE)
  b <- lapply(d_t, diff)
  names(b) <- NULL
  lbs <- gsub("-d", "-C", names(unlist(b)))
  ilbs <- unique(gsub("(.*)\\-(.*)", "\\1", lbs))
  clbs <- lapply(oj, function(oj) paste0("C", 1L:oj))
  ## delta method for slope / intercept to IRT
  if(vcov) {
    vc <- vcov(object)
    sel <- which(estnms %in% estnms[c(grep("a1$", estnms),
      grep("-d[[:digit:]]$", estnms))])
    vc <- vc[sel, sel]
    estnms <- colnames(vc)
    if(object$type == "PCM") {
      ## generate a full duplicated vcov
      oj_c <- oj
      oj_c[1L] <- oj_c[1L] + 1L
      ojc_c <- cumsum(oj_c)
      sel_r <- mapply(":", ojc_c - (oj_c - 1L), ojc_c, SIMPLIFY = FALSE)
      ndup <- cumsum(oj + 1L)
      sel_c <- mapply(":", ndup - oj, ndup, SIMPLIFY = FALSE)
      tmp <- matrix(0, K + 1L, max(ndup))
      for(s in 1L:N) {
        ind_r <- sel_r[[s]]
        ind_c <- sel_c[[s]]
        if(s == 1L) {
          diag(tmp[ind_r, ind_c]) <- 1
        }
        else {
          if(is.matrix(tmp[ind_r, ind_c[-1L]])) {
            diag(tmp[ind_r, ind_c[-1L]]) <- 1L
          } else {
            tmp[ind_r, ind_c[-1L]] <- 1L
          }
          tmp[1L, ind_c[1L]] <- 1L
        }
      }
      vc <- t(tmp) %*% (vc %*% tmp)
      rownames(vc) <- colnames(vc) <- estnms[seq_along(estnms) %*% tmp]
    }
    tmp <- -1 / unlist(a)
    C <- diag(tmp, N + K, N + K) +
      diag(-tmp, N + K + 1L, N + K + 1L)[-(N + K + 1L), -1L]
    loc <- c(0L, cumsum(oj + 1L))
    for(j in 1L:N) {
      ind <- (loc[j] + 1L):loc[j + 1L]
      C[ind, ind][, 1L] <- c(1L, diff(d[[j]]) / (a[[j]][1L] ^ 2))
    }
    vc_t <- C %*% vc %*% t(C)
    rownames(vc_t) <- colnames(vc_t) <- colnames(vc)
    vc <- vc_t[grep("-d[[:digit:]]$", rownames(vc_t)), grep("-d[[:digit:]]$", colnames(vc_t))]
  }
  if(relative) {
    if(type == "mode") {
      if(!is.list(ref)) {
        if(is.null(ref)) {
          ref <- lapply(oj, seq)
        } else if(is.vector(ref) & is.character(ref)) {
          stopifnot(ref %in% unlist(clbs))
          ref <- lapply(clbs, function(j) which(j %in% ref))
        } else if(is.vector(ref) & is.numeric(ref)) {
          ref <- as.integer(ref)
          stopifnot(all(ref %in% 1L:min(oj)))
          ref <- split(rep(ref, N), rep(1L:N, length(ref)))
        } else if(is.matrix(ref) & is.numeric(ref)) {
          stopifnot(ncol(ref) == K & nrow(ref) == K)
          ref2 <- vector(mode = "list", length = N)
          for(j in 1L:N) ref2[[j]] <- ref
          ref <- ref2
        } else {
          stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
        }
      } else {
        if(length(ref) < N) {
          stop("Not enough restrictions provided in argument 'ref'.")
        } else {
          for(j in 1L:N) {
            if(is.null(ref[[j]])) {
              ref[[j]] <- 1L:oj[j]
            } else if(is.vector(ref[[j]]) & is.character(ref[[j]])) {
              stopifnot(ref %in% clbs[[j]])
              ref[[j]] <- which(clbs[[j]] %in% ref[[j]])
            } else if(is.vector(ref[[j]]) & is.numeric(ref[[j]])) {
              ref[[j]] <- as.integer(ref[[j]])
              stopifnot(ref[[j]] %in% 1L:oj[j])
            } else if(is.matrix(ref[[j]]) & is.numeric(ref[[j]])) {
              stopifnot(ncol(ref[[j]]) == oj[j] & nrow(ref[[j]]) == oj[j])
            } else {
              stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
            }
          }
        }
      }
      ## D = ref
      D <- diag(K)
      if(is.matrix(ref[[1L]])) {
        for(j in 1L:N) {
          D[c(0L, ojc)[j] + 1L:oj[j], c(0L, ojc)[j] + 1L:oj[j]] <- ref[[j]]
        }
      } else {
        for(j in 1L:N) {
          D[c(0L, ojc)[j] + 1L:oj[j], c(0L, ojc)[j] + ref[[j]]] <-
            D[c(0L, ojc)[j] + 1L:oj[j], c(0L, ojc)[j] + ref[[j]]] -
            1 / length(ref[[j]])
        }
      }
      ## apply D
      b <- as.vector(D %*% unlist(b))
      b <- split(b, rep(1L:N, oj))
      ## delta method D
      vc <-
      if(vcov) {
        D %*% vc %*% t(D)
      } else {
        matrix(NA_real_, K, K)
      }
      ## C = cumulative, apply and delta method
      if(cumulative) {
        b <- lapply(b, cumsum)
        C <- matrix(0, K, K)
        for(j in 1L:N) {
          for(k in 1L:oj[j]) C[c(0L, ojc)[j] + k, c(0L, ojc)[j] + 1L:k] <- 1
        }
        if(vcov) vc <- C %*% vc %*% t(C)
      }
    } else {
      stop("Relative threshold parameters not implemented for types other than mode.")
    }
  } else {
    if(type == "mode") {
      if(is.null(ref)) {
        ref <- NULL
        ## scale defined via means and vars
      } else if(is.vector(ref) & is.character(ref)) {
        stopifnot(all(ref %in% lbs))
        ref <- which(lbs %in% ref)
      } else if(is.vector(ref) & is.numeric(ref)) {
        ref <- as.integer(ref)
        stopifnot(all(ref %in% 1L:K))
      } else if(is.matrix(ref) & is.numeric(ref)) {
        stopifnot(nrow(ref) == K & ncol(ref) == K)
      } else {
        stop("Argument 'ref' can only be a character vector with threshold parameter labels or a numeric vector with threshold parameter indices.")
      }
      ## D = ref
      if(is.matrix(ref)) {
        D <- ref
      } else {
        D <- diag(K)
        D[, ref] <- D[, ref] - 1 / length(ref)
      }
      ## apply D
      b <- as.vector(D %*% unlist(b))
      b <- split(b, rep(1L:N, oj))
      ## delta method D
      vc <-
      if(vcov) {
        D %*% vc %*% t(D)
      } else {
        matrix(NA_real_, K, K)
      }
      ## C = cumulative, apply and delta method
      if(cumulative) {
        b <- lapply(b, cumsum)
        if(vcov) {
          C <- matrix(0, K, K)
          for(j in 1L:N) {
            for(k in 1L:oj[j]) C[c(0L, ojc)[j] + k, c(0L, ojc)[j] + 1L:k] <- 1
          }
          vc <- C %*% vc %*% t(C)
        }
      }
    } else {
      if(!is.null(ref)) {
        warning("Argument 'ref' is not processed for types other than mode.")
      }
      vc <- matrix(NA_real_, K, K)
      if(type == "median") {
        zmedian <- function(theta = NULL, a = NULL, b = NULL, geq = NULL,
          ncat = NULL) {
          p <- pgpcm(theta = theta, a = a, b = b)
          rowSums(p[, (geq + 1):ncat, drop = FALSE]) - 0.5
        }
        for(j in 1L:N) {
          b[[j]] <- sapply(1L:oj[j], function(geq) {
            uniroot(f = zmedian, interval = c(-10, 10), a = a[[j]][1],
              b = b[[j]], geq = geq, ncat = oj[j] + 1L)$root})
        }
      } else if(type == "mean") {
        xpct <- lapply(oj, function(oj) 1:oj - 0.5)
        zexpct <- function(theta = NULL, a = NULL, b = NULL,
          expct = NULL) {
          pgpcm(theta = theta, a = a, b = b) %*% 0:length(b) - expct
        }
        for(j in 1L:N) {
          b[[j]] <- sapply(xpct[[j]], function(xp) {
            uniroot(f = zexpct, interval = c(-10, 10), a = a[[j]][1L],
              b = b[[j]], expct = xp)$root})
        }
      }
      if(cumulative) b <- lapply(b, cumsum)
    }
  }
  names(b) <- ilbs
  rownames(vc) <- colnames(vc) <- lbs
  for(j in 1L:N) names(b[[j]]) <- paste0("C", 1L:oj[j])
  if(!alias & type == "mode") {
    if(is.matrix(ref)) {
    ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
      stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
    } else {
      if(relative) {
        i <- split(1L:K, rep(1L:N, oj))
        alias <- vector(mode = "list", length = N)
        for(j in 1L:N) {
          b[[j]] <- b[[j]][-ref[[j]][1L]]
          i[[j]] <- i[[j]][-ref[[j]][1L]]
          alias[[j]] <- ref[[j]][1L]
        }
        i <- unlist(i)
        vc <- vc[i, i]
        names(alias) <- ilbs
      } else {
        ref1 <- ref[1L]
        i <- split(1L:K, rep(1L:N, oj))
        item <- which(sapply(i, function(j) ref1 %in% j))
        b[[item]] <- b[[item]][-which(ref1 == i[[item]])]
        vc <- vc[-ref1, -ref1]
        alias <- paste0("I", item, "-b", which(ref1 == i[[item]]))
        names(alias) <- ilbs[item]
      }
    }
  }
  rv <- structure(b, class = "threshpar", model = "GPCM", type = type,
    ref = ref, relative = relative, cumulative = cumulative, alias = alias, vcov = vc)
  return(rv)
}



discrpar.gpcmodel <- function(object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  coefs <- coef(object)
  a <-
  if(object$type == "PCM") {
    rep(coefs[grep("a1$", names(coefs))], N)
  } else {
    coefs[grep("-a1$", names(coefs))]
  }
  lbs <- colnames(object$data)
  vc <-
  if(vcov) {
    tmp <- vcov(object)
    if(object$type == "PCM") {
      a_var <- tmp[1L, 1L]
      tmp <- matrix(1, N, N)
      diag(tmp) <- a_var
    } else {
      tmp[grep("-a1$", rownames(tmp)), grep("-a1$", colnames(tmp))]
    }
  } else {
    matrix(NA_real_, N, N)
  }
  ## handling of model type
  if(object$type == "PCM" & !is.null(ref)) {
    warning("Argument 'ref' is currently not processed.")
    ref <- NULL
  }
  ## handling of negative slopes
  if(any(a < 0)) {
    if(!is.null(ref)) {
      stop("'ref' cannot be applied due to negative slopes.")
    } else {
      names(a) <- rownames(vc) <- colnames(vc) <- lbs
    }
  } else {
    if(is.null(ref)) {
      ref <- NULL
      ## scale defined via means and vars
    } else if(is.vector(ref) & is.character(ref)) {
      stopifnot(all(ref %in% lbs))
      ref <- which(lbs %in% ref)
    } else if(is.vector(ref) & is.numeric(ref)) {
      ref <- as.integer(ref)
      stopifnot(all(ref %in% 1L:N))
    } else if(is.matrix(ref) & is.numeric(ref)) {
      stopifnot(nrow(ref) == N & ncol(ref) == N)
    } else {
      stop("Argument 'ref' is misspecified (see ?discrpar for possible values).")
    }
    ## D = ref is applied on the log scale --> ratio based ref
    if(is.matrix(ref)) {
      D <- ref
    } else {
      D <- diag(N)
      D[, ref] <- D[, ref] - 1 / length(ref)
    }
    ## delta method D
    if(vcov) {
      C <- diag(1 / a, N, N)
      tmp <- C %*% vc %*% t(C)
      tmp <- D %*% tmp %*% t(D)
      C <- diag(exp(as.vector(D %*% log(a))), N, N)
      vc <- C %*% tmp %*% t(C)
    }
    ## apply D
    a <- exp(as.vector(D %*% log(a)))
    names(a) <- rownames(vc) <- colnames(vc) <- lbs
    if(!alias) {
      if(object$type == "PCM") {
        alias <- rep(1, N)
        a <- numeric()
        vc <- matrix(0, 0L, 0L)
        names(alias) <- lbs
      } else if(is.matrix(ref)) {
        ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
        stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
      } else if(!is.null(ref)) {
        a <- a[-ref[1L]]
        vc <- vc[-ref[1L], -ref[1L]]
        alias <- paste0("I", ref[1L])
        names(alias) <- lbs[ref[1L]]
      }
    }
  }
  rv <- structure(a, class = "discrpar", model = "GPCM", ref = ref,
    alias = alias, vcov = vc)
  return(rv)
}



guesspar.gpcmodel <- function(object, alias = TRUE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  lbs <- colnames(object$data)
  if(alias) {
    g <- rep(0, N)
    vc <-
    if(vcov) {
      matrix(0, N, N)
    } else {
      matrix(NA_real_, N, N)
    }
    names(g) <- rownames(vc) <- colnames(vc) <- lbs
  } else {
    g <- numeric()
    vc <- matrix(0, 0L, 0L)
    alias <- rep(1, N)
    names(alias) <- lbs
  }
  rv <- structure(g, class = "guesspar", model = "GPCM", alias = alias, vcov = vc)
  return(rv)
}



upperpar.gpcmodel <- function(object, alias = TRUE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  lbs <- colnames(object$data)
  if(alias) {
    u <- rep(1, N)
    vc <-
    if(vcov) {
      matrix(0, N, N)
    } else {
      matrix(NA_real_, N, N)
    }
    names(u) <- rownames(vc) <- colnames(vc) <- lbs
  } else {
    u <- numeric()
    vc <- matrix(0, 0L, 0L)
    alias <- rep(1, N)
    names(alias) <- lbs
  }
  rv <- structure(u, class = "upperpar", model = "GPCM", alias = alias,
    vcov = vc)
  return(rv)
}



personpar.gpcmodel <- function(object, personwise = FALSE, vcov = TRUE,
  interval = NULL, tol = 1e-6, method = "EAP", ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for gpcmodels.", call. = FALSE)
  }

  ## dropped ref
  ## changed default method to EAP
  ## FIXME: full vcov not provided by mirt
  if(is.null(interval)) {
    interval <- c(-6, 6)
  }
  ## personwise matches each person with their parameter
  ## otherwise group means and vars
  if(personwise) {
    tmp <- mirt::fscores(object$mirt, method = method,
      full.scores = TRUE, theta_lim = interval, gradtol = tol, ...)[, 1L]
    nms <- seq_along(tmp)
    vc <- matrix(NA_real_, length(tmp), length(tmp))
    rownames(vc) <- colnames(vc) <- nms
    rv <- structure(tmp, .Names = nms, class = "personpar", model = "GPCM",
      vcov = vc, type = "personwise")
  } else {
    tmp <- mirt::mod2values(object$mirt)
    tmp <- tmp[tmp$class == "GroupPars", ]
    nms <- paste(tmp$group, tmp$name, sep = "-")
    if(vcov & (length(nms) > 2L) & attributes(vcov(object))$method != "none") {
      sel <- paste(tmp$name, tmp$parnum, sep = ".")[-(1L:2L)]
      vc <- mirt::vcov(object$mirt)
      vc <- vc[rownames(vc) %in% sel, colnames(vc) %in% sel]
      vc <- rbind(0, cbind(0, vc)) ## zero for reference variance
      vc <- rbind(0, cbind(0, vc)) ## zero for reference mean
    } else if(vcov & length(nms) == 2L) {
      vc <- matrix(0, length(tmp[, 6L]), length(tmp[, 6L]))
    } else {
      vc <- matrix(NA_real_, length(tmp[, 6L]), length(tmp[, 6L]))
    }
    rownames(vc) <- colnames(vc) <- nms
    rv <- structure(tmp[, 6L], .Names = nms, class = "personpar", model = "GPCM",
      vcov = vc, type = "normal")
  }
  return(rv)
}



nobs.gpcmodel <- function(object, ...)
{
  object$n
}



bread.gpcmodel <- function(x, ...)
{
  x$n * vcov(x, ...)
}



## unexported: estfun interface for mirt:::estfun.AllModelClass
estfun_gpcmodel <- function(x, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for gpcmodels.", call. = FALSE)
  }

  ## call mirt
  scores <- mirt::estfun.AllModelClass(x$mirt)
  if(!x$grouppars & !is.null(x$impact)) {
    scores <- scores[, seq_along(x$coefficients)]
  }
  colnames(scores) <- names(x$coefficients)
  return(scores)
}

estfun.gpcmodel <- function(x, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for gpcmodels.", call. = FALSE)
  }

  ## completely cleaned (lowest category zero, null cats treatment, weights) data
  dat <- x$data
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  ## groups
  impact <- x$impact
  g <-
  if(is.null(impact)) {
    as.factor(rep("all", x$n))
  } else {
    impact
  }
  G <- length(levels(g))
  g <- as.numeric(g)
  ## coefficients
  coefs <- coef(x)
  ## names of the estimated parameters
  estnms <- names(coefs)
  ## quadratures
  AX <- do.call(cbind, mirt::extract.mirt(x$mirt, "Prior"))
  q <- dim(AX)[1L]
  X <- seq(-6, 6, length.out = q)             ## no. quads, see mirt
  ## relevant infos
  M <- x$n                                    ## no. pers
  N <- sum(x$items)                           ## no. items
  p <- apply(dat, 2L, max, na.rm = TRUE) + 1L ## no. cats
  K <- max(p)                                 ## max no. cats
  ## indices:
  ## groups      a - G
  ## persons     i - M
  ## items       j - N
  ## categories  k - p or K, if two times needed also l - p or K
  ## quadratures f - q
  ## means and variances of the groups
  chars <- personpar(x, vcov = FALSE)
  means <- chars[seq(from = 1L, by = 2L, length.out = G)]
  vars <- chars[seq(from = 2L, by = 2L, length.out = G)]
  ## scoring values
  ak <- lapply(x$categories, function(ct) c(0, ct))
  ## estimated slopes
  aest <- discrpar(x, vcov = FALSE)
  ## estimated intercepts & d0
  dest <- coefs[grep("-d[[:digit:]]$", estnms)]
  dest <- lapply(split(dest, rep(1L:N, p - 1L)), function(d) c(0, d))
  ## pX posterior for a person showing the observed pattern given f
  ## LX placeholder for pattern probability
  ## scores_a empirical estimation function for a
  ## scores_d empirical estimation function for d
  pX <- LX <- matrix(0, M, q)
  scores_a <- scores_d <- vector("list", G)
  scores <- matrix(0, M, sum(p))
  for(a in 1L:G) {
    ## px_tmp_ probability for a person showing the observed pattern
    ## given f, computed iteratively
    px_tmp_ <- matrix(1, sum(g == a), q)
    scores_a[[a]] <- scores_d[[a]] <- vector("list", N)
    M_tmp <- sum(g == a)
    for(j in 1L:N) {
      ## relevant data and parameters
      pat <- dat[g == a, j] + 1
      ak_tmp <- ak[[j]]
      aest_tmp <- aest[j]
      dest_tmp <- dest[[j]]
      ## px_tmp probability for a person falling into any k of j given q
      exp_tmp <- exp(outer(ak_tmp, aest_tmp * X, "*") +
        matrix(dest_tmp, p[j], q))
      px_tmp <- exp_tmp / matrix(colSums(exp_tmp), p[j], q, TRUE)
      ones <- which(px_tmp == 1, arr.ind = TRUE)
      zeros <- which(px_tmp == 0, arr.ind = TRUE)
      if(length(ones) > 0) {
        px_tmp[ones] <- 1 - .Machine$double.neg.eps
      }
      if(length(zeros) > 0) {
        px_tmp[zeros] <- .Machine$double.neg.eps
      }
      ## select the matching k
      px_tmp_sel <- px_tmp[pat, ]
      ## handle NAs
      ex <- is.na(px_tmp_sel)
      if(any(ex)) {
        px_tmp_sel[which(ex, arr.ind = TRUE)] <- 1
      }
      ## update px_tmp_
      px_tmp_ <- px_tmp_ * px_tmp_sel
      ## evaluated 1st partial derivative of the model with respect to a
      px_da_tmp <- ((exp_tmp * ak_tmp * matrix(X, p[j], q, TRUE) *
        matrix(colSums(exp_tmp), p[j], q, TRUE)) - (exp_tmp *
        matrix(colSums(exp_tmp * outer(ak_tmp, X, "*")), p[j], q, TRUE))) /
        (matrix(colSums(exp_tmp), p[j], q, TRUE) *
        matrix(colSums(exp_tmp), p[j], q, TRUE))
      scores_a[[a]][[j]] <- px_da_tmp[pat, ] / px_tmp[pat, ]
      ## evaluated 1st partial derivative of the model with respect to d
      px_dd_tmp <- array(0, dim = c(p[j], p[j], q))
      for(k in 1L:p[j]) {
        px_dd_tmp[k, k, ] <- px_tmp[k, ] - (px_tmp[k, ] * px_tmp[k, ])
        px_dd_tmp[k, -k, ] <- - (exp_tmp[-k, ] *
          matrix(exp_tmp[k, ], p[j] - 1, q, TRUE)) /
          (matrix(colSums(exp_tmp), p[j] - 1, q, TRUE) *
          matrix(colSums(exp_tmp), p[j] - 1, q, TRUE))
      }
      scores_d[[a]][[j]] <- array(px_dd_tmp[pat, -1, ],
        dim = c(M_tmp, p[j] - 1, q)) / aperm(array(px_tmp[pat, ],
        dim = c(M_tmp, q, p[j] - 1)), c(1L, 3L, 2L))
    }
    ## calculate LX and pX
    LX[g == a, ] <- px_tmp_
    pX[g == a, ] <- px_tmp_ * matrix(AX[, a], M_tmp, q, TRUE)
    pX[g == a, ] <- pX[g == a, ] / rowSums(pX[g == a, , drop = FALSE])
    ## finally calculate the scores using pX
    scores_a[[a]] <-
    lapply(scores_a[[a]], function(sc_a) {
      rowSums(sc_a * pX[g == a, , drop = FALSE])
    })
    scores_d[[a]] <-
    lapply(scores_d[[a]], function(sc_d) {
      colSums(aperm(sc_d, c(3L, 1L, 2L)) *
        array(t(pX[g == a, , drop = FALSE]), dim = c(q, M_tmp, dim(sc_d)[2L])))
    })
    ## combine the scores into one matrix
    text <- paste0("cbind(", paste0(
      sapply(1L:N, function(j) {
        paste0(gsub("item", j, "scores_a[[a]][[item]]"),
          ",", gsub("item", j, "scores_d[[a]][[item]]"))
      }), collapse = ","), ")")
    scores[g == a, ] <- eval(parse(text = text))
  }
  ## restriction to pcmodel (scores are the rowsums over the restrictions)
  if(x$type == "PCM") {
    ind <- cumsum(p) - (p - 1)
    tmp <- scores
    scores <- tmp[, -ind[-1L]]
    scores[, 1L] <- rowSums(tmp[, ind])
  }
  ## grouppars
  if(x$grouppars) {
    ## PX posterior for a person showing the observed pattern over all fs
    PX <- numeric(M)
    for(a in 2L:G) {
      PX[g == a] <- rowSums(LX[g == a, ] *
        matrix(AX[, a], sum(g == a), q, TRUE))
    }
    for(a in 2L:G) {
      scores <- cbind(scores, matrix(0, M, 2L))
      m <- means[a]
      v <- vars[a]
      ## scores for the group mean, Baker & Kim, 2004, p. 274, Eq. (10.36)
      scores[g == a, dim(scores)[2L] - 1L] <- ((v ^ -1) * (PX[g == a] ^ -1) *
        colSums(matrix(X - m, q, sum(g == a)) * t(LX[g == a, , drop = FALSE]) *
          matrix(AX[, a], q, sum(g == a))))
      ## scores for the group variance, Baker & Kim, 2004, p. 274, Eq. (10.37)
      scores[g == a, dim(scores)[2L]] <- (- (1 / 2) * (PX[g == a] ^ -1) *
        colSums(matrix((v ^ -1) - ((X - m) ^ 2) * (v ^ -2), q, sum(g == a)) *
          t(LX[g == a, , drop = FALSE]) * matrix(AX[, a], q, sum(g == a))))
    }
  }
  ## handle weights
  scores <- scores * weights
  ## handle NAs
  scores[which(is.na(scores), arr.ind = TRUE)] <- 0
  ## rename the columns
  colnames(scores) <- estnms
  return(scores)
}



pgpcm <- function(theta = NULL, a = NULL, b = NULL)
{
  ## probabilities under the model (in IRT formulation) given theta
  stopifnot(!is.null(theta) & !is.null(a) & !is.null(b))
  stopifnot(mode(a) == mode(b))
  if(is.list(theta)) {
    return(lapply(theta, pgpcm, a = a, b = b))
  }
  if(is.list(a) & is.list(b)) {
    return(mapply(pgpcm, a, b, MoreArgs = list(theta = theta),
      SIMPLIFY = FALSE))
  }
  num <- cbind(0, a * outer(theta, b, "-"))
  num <- t(exp(apply(num, 1L, cumsum)))
  denom <- rowSums(num)
  return(num / denom)
}



rgpcm <- function(theta, a, b, nullcats = FALSE, return_setting = TRUE)
{
  ## sample data under the model (in IRT formulation) given theta
  stopifnot(mode(a) == mode(b))
  if(is.list(theta)) {
    return(lapply(theta, rgpcm, a = a, b = b, nullcats = nullcats,
      return_setting = return_setting))
  }
  probs <- pgpcm(theta = theta, a = a, b = b)
  if(!is.list(probs)) {
    probs <- list(probs)
  }
  rsp_item_i <- function(probmat_i) {
    oj <- ncol(probmat_i)
    rsp <- apply(probmat_i, 1L, function(p) sample.int(oj, 1L, prob = p))
    if(!nullcats) {
      while(length(unique.default(rsp)) != oj) {
        rsp <- apply(probmat_i, 1L, function(p) sample.int(oj, 1L, prob = p))
      }
    }
    return(rsp - 1)
  }
  res <- lapply(probs, rsp_item_i)
  res <- do.call(cbind, res)
  if(return_setting) {
    return(list(a = a, b = b, theta = theta, data = res))
  } else {
    return(res)
  }
}

