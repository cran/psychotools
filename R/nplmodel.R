# estimate an n-PL model in slope / intercept formulation using mirt (MML & EM)

plmodel <- ## for backward compatibility
nplmodel <- function(y, weights = NULL, impact = NULL,
  type = c("2PL", "3PL", "3PLu", "4PL", "1PL", "RM"), grouppars = FALSE, vcov = TRUE,
  start = NULL, method = "BFGS", maxit = 500L, reltol = 1e-5, ...)
{

  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for fitting a nplmodel.", call. = FALSE)
  }

  ## handle arguments and defaults
  type <- match.arg(type, c("2PL", "3PL", "3PLu", "4PL", "1PL", "RM"))
  if(type == "RM") type <- "1PL"

  itemtype <- if(type == "1PL") "2PL" else type

  if(is.logical(vcov)) vcov <- c("none", "Oakes")[vcov + 1L]
  vcov <- match.arg(vcov, c("Oakes", "Richardson", "forward", "central", "crossprod",
    "Louis", "sandwich", "sandwich.Louis", "SEM", "Fisher", "complete", "none"))

  ## response matrix
  y <- as.matrix(y)
  ## data should be dichotomous
  if(any(apply(y, 2, function(x) sum(unique(x) != "NA", na.rm = TRUE) > 2))) {
    stop("'y' should only contain dichotomous data.")
  }
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

  ## move categories to zero if necessary
  mincat <- apply(y, 2L, min, na.rm = TRUE)
  if(any(mincat > 0)) {
    warning("Minimum score is not zero for all items (", paste(which(mincat > 0), collapse = ", "), "). These items are scaled to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }

  ## call mirt, mostly defaults
  ## if impact --> multiple group model, reference group mean 0 and var 1,
  ## other groups free, item parameters restricted across whole sample
  ## otherwise --> simple mirt model, group mean 0 and var 1
  model <-
  if(!is.null(impact)) {
    pars <- mirt::multipleGroup(data = y, model = 1L, group = impact,
      itemtype = itemtype, invariance = c("free_means", "free_var", colnames(y)),
      pars = "values")

    constrain <- if(type == "1PL") list(which(pars$name == "a1")) else NULL

    if(!is.null(start)) {
      if(type == "1PL") {
        ## start should be a vector of item parameters ordered itemwise
        ## and means and vars appended groupwise
        ## a list can be provided containing a single slope,
        ## a vector of intercepts, a vector of means and a vector of vars
        if(is.list(start)) {
          n <- length(start)
          start <- c(start[[1L]], start[[2L]],
            as.vector(do.call(rbind, start[(n - 1L):n])))
        }
        if(length(start) != sum(pars$est[pars$group == lvls[1L]]) +
          (2L * (nlvls - 1L)) - (N - 1L)) {
          stop("Argument 'start' is misspecified (see ?nplmodel for possible values).")
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
        ## a vector of intercepts, a vector of guessing parameters,
        ## a vector of upper asymptotes (depending on the restriction),
        ## a vector of means and a vector of vars
        if(is.list(start)) {
          n <- length(start)
          start <- c(as.vector(do.call(rbind, start[-((n - 1L):n)])),
            as.vector(do.call(rbind, start[(n - 1L):n])))
        }
        if(length(start) != sum(pars$est[pars$group == lvls[1L]]) +
          (2L * (nlvls - 1L))) {
          stop("Argument 'start' is misspecified (see ?nplmodel for possible values).")
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
      itemtype = itemtype, invariance = c("free_means", "free_var", colnames(y)),
      SE = (vcov != "none"), SE.type = if(vcov == "none") "Oakes" else vcov,
      optimizer = method, calcNull = FALSE, TOL = reltol,
      technical = list(NCYCLES = maxit), verbose = FALSE, pars = pars,
      constrain = constrain, survey.weights = weights, ...)
  } else {
    grouppars <- FALSE
    pars <- mirt::mirt(data = y, model = 1L, itemtype = itemtype, pars = "values")

    constrain <- if(type == "1PL") list(which(pars$name == "a1")) else NULL

    if(!is.null(start)) {
      if(type == "1PL") {
        ## start should be a vector of item parameters ordered itemwise
        ## a list can be provided containing a single slope and
        ## a vector of intercepts
        if(is.list(start)) {
         start <- c(start[[1L]], start[[2L]])
        }
        if(length(start) != sum(pars$est) - (N - 1L)) {
          stop("Argument 'start' is misspecified (see ?nplmodel for possible values).")
        }
        pars$value[pars$est & pars$name == "a1"] <- start[1L]
        pars$value[pars$est & pars$name != "a1"] <- start[-1L]
      } else {
        ## start should be a vector of item parameters ordered itemwise
        ## a list can be provided containing a vector of slopes,
        ## a vector of intercepts, a vector of guessing parameters and
        ## a vector of upper asymptotes (depending on the restriction)
        if(is.list(start)) {
         start <- as.vector(do.call(rbind, start))
        }
        if(length(start) != sum(pars$est)) {
          stop("Argument 'start' is misspecified (see ?nplmodel for possible values).")
        }
        pars$value[pars$est] <- start
      }
    }
    mirt::mirt(data = y, model = 1L, itemtype = itemtype, SE = (vcov != "none"),
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
    warning(paste0("Negative slopes were estimated (", paste(neg_a, sep = "", collapse = ", "), ")."))
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

  ## g and u were estimated on the logit scale
  if(type %in% c("3PL", "4PL")) {
    gx <- grep("-g$", names(est))
    est[gx] <- qlogis(est[gx])
    names(est)[gx] <- gsub("-g$", "-logit\\(g\\)", names(est[gx]))
  }
  if(type %in% c("3PLu", "4PL")) {
    ux <- grep("-u$", names(est))
    est[ux] <- qlogis(est[ux])
    names(est)[ux] <- gsub("-u$", "-logit\\(u\\)", names(est[ux]))
  }

  ## if grouppars = FALSE and multiple group model drop the group parameters
  if(!grouppars && !is.null(impact)) {
    estnms <- names(est)
    estnms <- estnms[-c(grep("MEAN", estnms), grep("COV", estnms))]
    est <- est[which(names(est) %in% estnms)]
  }

  ## handle type of model restriction
  if(type == "1PL") {
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
    n = sum(weights_org > 0),
    n_org = M_org,
    weights = if(identical(as.vector(weights_org), rep.int(1L, M_org))) NULL else weights_org,
    na = any_na,
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
  class(res) <- "nplmodel"
  return(res)
}



print.nplmodel <- function(x, digits = max(3L, getOption("digits") - 3L), logit = FALSE, ...)
{
  msg <-
  if(x$type == "1PL" | x$type == "2PL") {
    list(mgm = "slopes, intercepts and estimated group parameters:\n",
      m = "slopes and intercepts:\n")
  } else if(x$type == "3PL") {
    list(mgm = "slopes, intercepts, guessing parameters and estimated group parameters:\n",
      m = "slopes, intercepts and guessing parameters:\n")
  } else if(x$type == "3PLu") {
    list(mgm = "slopes, intercepts, upper asymptotes and estimated group parameters:\n",
      m = "slopes, intercepts and upper asymptotes:\n")
  } else {
    list(mgm = "slopes, intercepts, guessing parameters, upper asymptotes and estimated group parameters:\n",
      m = "slopes, intercepts, guessing parameters and upper asymptotes:\n")
  }
  if(x$grouppars) {
    msg <- msg[[1L]]
  } else {
    msg <- msg[[2L]]
  }
  cat(paste(x$type, msg))
  print(coef(x, logit = logit), digits = digits, ...)
  invisible(x)
}



coef.nplmodel <- function(object, logit = FALSE, ...)
{
  cf <- object$coefficients
  if(!logit) {
    if(object$type %in% c("3PL", "4PL")) {
      gx <- grep("-logit\\(g\\)$", names(cf))
      cf[gx] <- plogis(cf[gx])
      names(cf)[gx] <- gsub("-logit\\(g\\)$", "-g", names(cf[gx]))
    }
    if(object$type %in% c("3PLu", "4PL")) {
      ux <- grep("-logit\\(u\\)$", names(cf))
      cf[ux] <- plogis(cf[ux])
      names(cf)[ux] <- gsub("-logit\\(u\\)$", "-u", names(cf[ux]))
    }
  }
  cf
}



vcov.nplmodel <- function(object, logit = TRUE, ...)
{
  vc <- object$vcov
  method <- attributes(object$vcov)$method
  ## delta method
  if(!logit) {
    cf <- object$coefficients
    estnms <- names(cf)
    C <- diag(1, length(cf), length(cf))
    if(object$type %in% c("3PL", "4PL")) {
      gx <- grep("-logit\\(g\\)$", estnms)
      diag(C)[gx] <- exp(cf[gx]) / ((1 + exp(cf[gx])) ^ 2)
      estnms[gx] <- gsub("-logit\\(g\\)$", "-g", estnms[gx])
    }
    if(object$type %in% c("3PLu", "4PL")) {
      ux <- grep("-logit\\(u\\)$", estnms)
      diag(C)[ux] <- exp(cf[ux]) / ((1 + exp(cf[ux])) ^ 2)
      estnms[ux] <- gsub("-logit\\(u\\)$", "-u", estnms[ux])
    }
    vc <- C %*% vc %*% t(C)
    rownames(vc) <- colnames(vc) <- estnms
  }
  structure(vc, method = method)
}



## based on stats::confint.default
confint.nplmodel <- function(object, parm, level = 0.95, logit = TRUE, logistic_bounds = FALSE, ...)
{
  cf <- coef(object, logit = logit)
  estnms <- names(cf)
  if(missing(parm)) {
    parm <- estnms
  } else if (is.numeric(parm)) {
    parm <- estnms[parm]
  }
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object, logit = logit)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  if(logit & logistic_bounds) {
    if(object$type %in% c("3PL", "4PL")) {
      gx <- grep("-logit\\(g\\)$", estnms)
      ci[gx, ] <- plogis(ci[gx, ])
      rownames(ci)[gx] <- gsub("-logit\\(g\\)$", "-g", estnms[gx])
    }
    if(object$type %in% c("3PLu", "4PL")) {
      ux <- grep("-logit\\(u\\)$", estnms)
      ci[ux, ] <- plogis(ci[ux, ])
      rownames(ci)[ux] <- gsub("-logit\\(u\\)$", "-u", estnms[ux])
    }
  }
  ci
}



logLik.nplmodel <- function(object, ...)
{
  structure(object$loglik, df = object$df, class = "logLik")
}



weights.nplmodel <- function(object, ...)
{
  if(is.null(object$weights)) rep.int(1L, object$n_org) else object$weights
}



summary.nplmodel <- function(object, vcov. = NULL, ...)
{
  pcf <- coef(object, logit = FALSE)
  vc <-
  if(is.null(vcov.)) {
    vcov(object, logit = TRUE)
  } else {
    if(is.function(vcov.)) {
      vcov.(object, logit = TRUE)
    } else {
      vcov.
    }
  }
  if(object$type %in% c("3PL", "3PLu", "4PL")) {
    lcf <- lcf_drop <- coef(object, logit = TRUE)
    lcf_drop[-grep("-logit\\(g\\)|-logit\\(u\\)", names(lcf_drop))] <- NA
    cf <- cbind(pcf, lcf_drop, sqrt(diag(vc)), lcf / sqrt(diag(vc)),
      2 * pnorm(-abs(lcf / sqrt(diag(vc)))))
    colnames(cf) <- c("Estimate", "Logit Estim.", "Std. Error", "z value", "Pr(>|z|)")
  } else {
    cf <- cbind(pcf, sqrt(diag(vc)), pcf / sqrt(diag(vc)),
      2 * pnorm(-abs(pcf / sqrt(diag(vc)))))
    colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
  object$coefficients <- cf
  class(object) <- "summary.nplmodel"
  return(object)
}



print.summary.nplmodel <- function(x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat(sprintf("\n%s\n\n",
      if(x$type == "1PL") {
        "One parameter logistic model"
      } else if(x$type == "2PL") {
        "Two parameter logistic model"
      } else if(x$type == "3PL") {
        "Three parameter logistic model"
      } else if(x$type == "3PLu") {
        "Three parameter logistic model (upper asymptotes)"
      } else if(x$type == "4PL") {
        "Four parameter logistic model"
      }))
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")
  }
  msg <-
  if(x$type == "1PL" | x$type == "2PL") {
    list(mgm = "Slopes, intercepts and estimated group parameters:\n",
      m = "Slopes and intercepts:\n")
  } else if(x$type == "3PL") {
    list(mgm = "Slopes, intercepts, guessing parameters and estimated group parameters:\n",
      m = "Slopes, intercepts and guessing parameters:\n")
  } else if(x$type == "3PLu") {
    list(mgm = "Slopes, intercepts, upper asymptotes and estimated group parameters:\n",
      m = "Slopes, intercepts and upper asymptotes:\n")
  } else {
    list(mgm = "Slopes, intercepts, guessing parameters, upper asymptotes and estimated group parameters:\n",
      m = "Slopes, intercepts, guessing parameters and upper asymptotes:\n")
  }
  msg <- if(x$grouppars) msg[[1L]] else msg[[2L]]
  cat(msg)
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
    na.print = "NA", ...)
  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n")
  cat(sprintf("Number of EM cycles: %s\n", x$iterations))
  cat(sprintf("M-step optimizer: %s\n", x$method))
  invisible(x)
}



plot.nplmodel <- function(x, type = c("regions", "profile", "curves", "information", "piplot"), ...)
{
  type <- match.arg(type, c("regions", "profile", "curves", "information", "piplot"))
  switch(type, curves = curveplot(x, ...),
    regions = regionplot(x, ...), profile = profileplot(x, ...),
    information = infoplot(x, ...),
    piplot = piplot(x, ...))
}



predict.nplmodel <- function(object, newdata = NULL,
  type = c("probability", "cumprobability", "mode", "median", "mean",
  "category-information", "item-information", "test-information"), ref = NULL, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for nplmodels.", call. = FALSE)
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
  a <- discrpar(object, vcov = FALSE)
  b <- itempar(object, vcov = FALSE)
  g <- guesspar(object, logit = FALSE, vcov = FALSE)
  u <- upperpar(object, logit = FALSE, vcov = FALSE)
  probs <- ppl(theta = newdata, a = a, b = b, g = g, u = u)
  colnames(probs) <- names(b)
  if(type %in% c("probability", "cumprobability", "category-information",
    "item-information", "test-information")) {
    clnms <- colnames(probs)
    rwnms <- rownames(probs)
    nc <- ncol(probs)
    probs0 <- matrix(0, nrow(probs), 2L * nc)
    probs0[, seq(from = 2L, by = 2L, length.out = nc)] <- probs
    probs0[, seq(from = 1L, by = 2L, length.out = nc)] <- 1 - probs
    if(type == "cumprobability") {
      probs0[, seq(from = 1L, by = 2L, length.out = nc)] <-
      probs0[, seq(from = 1L, by = 2L, length.out = nc)] + probs
    }
    probs <- probs0
    rownames(probs) <- rwnms
    colnames(probs) <- as.vector(t(outer(clnms, c("C0", "C1"), paste,
      sep = if(type == "probability") "-" else ">=")))
  }
  ## call mirt
  if(grepl("information", type)) {
    N <- sum(object$items)
    if(type == "category-information") {
      info <- matrix(NA_real_, nrow(probs), ncol(probs))
      colnames(info) <- as.vector(t(outer(clnms, c("C0", "C1"), paste,
        sep = "-")))
    } else {
      info <- matrix(NA_real_, nrow(probs), N)
      colnames(info) <- clnms
    }
    for(j in 1L:N) {
      idx <- grepl(clnms[j], colnames(probs))
      iteminfo <- mirt::iteminfo(mirt::extract.item(object$mirt, j, 1L),
        newdata)
      if(type == "category-information") {
        info[, idx] <- probs[, idx, drop = FALSE] * iteminfo
      } else {
        info[, j] <- iteminfo
      }
    }
  }
  switch(type,
    "probability" = probs,
    "cumprobability" = probs,
    "mode" = round(probs),
    "median" = round(probs),
    "mean" = round(probs),
    "category-information" = info,
    "item-information" = info,
    "test-information" = matrix(rowSums(info), ncol = 1L))
}



itempar.nplmodel <- function(object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  cf <- unlist(threshpar(object))
  N <- sum(object$items)
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
  ## apply D
  cf <- as.vector(D %*% cf)
  ## delta method D
  vc <-
  if(vcov) {
    D %*% vcov(threshpar(object)) %*% t(D)
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
  rv <- structure(cf, class = "itempar", model = "PL", ref = ref,
    alias = alias, vcov = vc)
  return(rv)
}



threshpar.nplmodel <- function(object, type = c("mode", "median", "mean"),
  ref = NULL, alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = TRUE, ...)
{
  type <- match.arg(type, c("mode", "median", "mean"))
  N <- sum(object$items)
  oj <- rep(1L, N)
  ojc <- cumsum(oj)
  K <- max(ojc)
  coefs <- coef(object, logit = FALSE)
  estnms <- names(coefs)
  ## go from slope / intercept to IRT
  a <- as.list(discrpar(object, vcov = FALSE))
  d <- as.list(coefs[grep("-d$", estnms)])
  b <- mapply(function(x, y) - (x / y), d, a, SIMPLIFY = FALSE)
  lbs <- gsub("-d", "-C1", names(unlist(b)))
  ilbs <- unique(gsub("(.*)\\-(.*)", "\\1", lbs))
  clbs <- lapply(oj, function(oj) paste0("C", 1L:oj))
  ## delta method for slope / intercept to IRT
  if(vcov) {
    vc <- vcov(object)
    sel <- which(estnms %in% estnms[c(grep("a1$", estnms),
      grep("-d$", estnms))])
    vc <- vc[sel, sel]
    estnms <- colnames(vc)
    if(object$type == "1PL") {
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
        } else {
          tmp[ind_r, ind_c[-1L]] <- 1
          tmp[1L, ind_c[1L]] <- 1
        }
      }
      vc <- t(tmp) %*% (vc %*% tmp)
      rownames(vc) <- colnames(vc) <- estnms[seq_along(estnms) %*% tmp]
    }
    C <- matrix(0, N + K, N + K)
    loc <- c(0L, cumsum(oj + 1L))
    for(j in 1L:N) {
      ind <- (loc[j] + 1L):loc[j + 1L]
      C[ind[1L], ind[1L]] <- 1
      C[ind[2L], ind[1L]] <- d[[j]] / ((a[[j]]) ^ 2)
      C[ind[2L], ind[2L]] <- -1 / a[[j]]
    }
    vc_t <- C %*% vc %*% t(C)
    rownames(vc_t) <- colnames(vc_t) <- colnames(vc)
    vc <- vc_t[grep("-d$", rownames(vc_t)), grep("-d$", colnames(vc_t))]
  }
  if(relative) {
    if(is.null(ref)) {
      ref <- 1
    } else if(is.vector(ref) & is.character(ref)) {
      stopifnot(all(ref %in% clbs))
    } else if(is.vector(ref) & is.numeric(ref)) {
      stopifnot(all(as.integer(ref) %in% 1L))
    } else if(is.matrix(ref) & is.numeric(ref)) {
      stopifnot(nrow(ref) == N & ncol(ref) == N)
    } else {
      stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
    }
    b <- rep(0, N)
    if(is.matrix(ref)) {
      b <- as.vector(ref %*% b)
    }
    b <- as.list(b)
    if(vcov) {
      vc <- matrix(0, N, N)
    } else {
      vc <- matrix(NA_real_, N, N)
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
    } else if(is.matrix(ref)) {
      stopifnot(nrow(ref) == N & ncol(ref) == N)
    } else {
      stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
    }
    ## D = ref
    if(is.matrix(ref)) {
      D <- ref
    } else {
      D <- diag(N)
      D[, ref] <- D[, ref] - 1 / length(ref)
    }
    ## apply D
    b <- as.list(D %*% unlist(b))
    ## delta method D
    vc <-
    if(vcov) {
      D %*% vc %*% t(D)
    } else {
      matrix(NA_real_, N, N)
    }
  }
  names(b) <- ilbs
  rownames(vc) <- colnames(vc) <- lbs
  for(j in 1L:N) names(b[[j]]) <- "C1"
  if(!alias) {
    if(is.matrix(ref)) {
    ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
      stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
    } else {
      if(relative) {
        b <- vector(mode = "list", length = N)
        for(j in 1L:N) b[[j]] <- numeric()
        vc <- matrix(0, 0L, 0L)
        alias <- as.list(rep(1, N))
        names(alias) <- ilbs
      } else {
        b <- b[-ref[1L]]
        vc <- vc[-ref[1L], -ref[1L]]
        alias <- paste0("I", ref[1L], "-C", ref[1L])
        names(alias) <- ilbs[ref[1L]]
      }
    }
  }
  rv <- structure(b, class = "threshpar", model = "PL", type = type,
    ref = ref, relative = relative, cumulative = cumulative, alias = alias, vcov = vc)
  return(rv)
}



discrpar.nplmodel <- function(object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  coefs <- coef(object, logit = FALSE)
  a <-
  if(object$type == "1PL") {
    rep(coefs[grep("a1$", names(coefs))], N)
  } else {
    coefs[grep("-a1$", names(coefs))]
  }
  lbs <- colnames(object$data)
  vc <-
  if(vcov) {
    tmp <- vcov(object)
    if(object$type == "1PL") {
      a_var <- tmp[1L, 1L]
      tmp <- matrix(a_var, N, N)
    } else {
      tmp[grep("-a1$", rownames(tmp)), grep("-a1$", colnames(tmp))]
    }
  } else {
    matrix(NA_real_, N, N)
  }
  ## handling of model type
  if(object$type == "1PL" & !is.null(ref)) {
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
      if(object$type == "1PL") {
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
  rv <- structure(a, class = "discrpar", model = "PL", ref = ref,
                  alias = alias, vcov = vc)
  return(rv)
}



guesspar.nplmodel <- function(object, alias = TRUE, logit = FALSE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  lbs <- colnames(object$data)
  if(object$type %in% c("1PL", "2PL", "3PLu")) {
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
  } else {
    coefs <- coef(object, logit = logit)
    gx <- grep("-g$|-logit\\(g\\)$", names(coefs))
    g <- coefs[gx]
    vc <-
    if(vcov) {
      vcov(object, logit = logit)[gx, gx]
    } else {
      matrix(NA_real_, N, N)
    }
    names(g) <- rownames(vc) <- colnames(vc) <- lbs
  }
  rv <- structure(g, class = "guesspar", model = "PL", alias = alias, logit = logit, vcov = vc)
  return(rv)
}



upperpar.nplmodel <- function(object, alias = TRUE, logit = FALSE, vcov = TRUE, ...)
{
  N <- sum(object$items)
  lbs <- colnames(object$data)
  if(object$type %in% c("1PL", "2PL", "3PL")) {
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
  } else {
    coefs <- coef(object, logit = logit)
    ux <- grep("-u$|-logit\\(u\\)$", names(coefs))
    u <- coefs[ux]
    vc <-
    if(vcov) {
      vcov(object, logit = logit)[ux, ux]
    } else {
      matrix(NA_real_, N, N)
    }
    names(u) <- rownames(vc) <- colnames(vc) <- lbs
  }
  rv <- structure(u, class = "upperpar", model = "PL", alias = alias, logit = logit, vcov = vc)
  return(rv)
}



personpar.nplmodel <- function(object, personwise = FALSE, vcov = TRUE,
  interval = NULL, tol = 1e-6, method = "EAP", ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for nplmodels.", call. = FALSE)
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
    rv <- structure(tmp, .Names = nms, class = "personpar", model = "PL",
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
    rv <- structure(tmp[, 6L], .Names = nms, class = "personpar", model = "PL",
      vcov = vc, type = "normal")
  }
  return(rv)
}



nobs.nplmodel <- function(object, ...)
{
  object$n
}



bread.nplmodel <- function(x, ...)
{
  x$n * vcov(x, ...)
}



## unexported: estfun interface for mirt:::estfun.AllModelClass
estfun_nplmodel <- function(x, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for nplmodels.", call. = FALSE)
  }

  ## call mirt
  scores <- mirt::estfun.AllModelClass(x$mirt)
  if(!x$grouppars & !is.null(x$impact)) {
    scores <- scores[, seq_along(x$coefficients)]
  }
  colnames(scores) <- names(x$coefficients)
  return(scores)
}

estfun.nplmodel <- function(x, ...)
{
  ## check for mirt
  if(!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' ist required for nplmodels.", call. = FALSE)
  }

  ## logit = TRUE (g and u parameters) as in mirt
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
  ## coefficients (not logit for g and u)
  coefs <- coef(x, logit = FALSE)
  ## names of the estimated parameters (logit for g and u)
  estnms <- names(coef(x, logit = TRUE))
  ## quadratures
  AX <- do.call(cbind, mirt::extract.mirt(x$mirt, "Prior"))
  q <- dim(AX)[1L]
  X <- seq(-6, 6, length.out = q)             ## no. quads, see mirt
  ## relevant infos
  M <- x$n                                    ## no. pers
  N <- sum(x$items)                           ## no. items
  p <- rep.int(2L, N)                         ## no. cats
  K <- 2L                                     ## max no. cats
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
  ## estimated slopes
  aest <- discrpar(x, vcov = FALSE)
  ## estimated intercepts
  dest <- coefs[grep("-d$", estnms)]
  ## (estimated) guessing parameters
  gest <- guesspar(x, logit = FALSE, vcov = FALSE)
  ## (estimated) upper asymptotes
  uest <- upperpar(x, logit = FALSE, vcov = FALSE)
  ## pX posterior for a person showing the observed pattern given f
  ## LX placeholder for pattern probability
  ## scores_a empirical estimation function for a
  ## scores_d empirical estimation function for d
  ## scores_g empirical estimation function for g
  ## scores_u empirical estimation function for u
  pX <- LX <- matrix(0, M, q)
  if(x$type == "3PL") {
    scores_a <- scores_d <- scores_g <- vector("list", G)
    nest <- 3L * N
  } else if(x$type == "3PLu") {
    scores_a <- scores_d <- scores_u <- vector("list", G)
    nest <- 3L * N
  } else if(x$type == "4PL") {
    scores_a <- scores_d <- scores_g <- scores_u <- vector("list", G)
    nest <- 4L * N
  } else {
    scores_a <- scores_d <- vector("list", G)
    nest <- 2L * N
  }
  scores <- matrix(0, M, nest)
  for(a in 1L:G) {
    ## px_tmp_ probability for a person showing the observed pattern
    ## given f, computed iteratively
    px_tmp_ <- matrix(1, sum(g == a), q)
    scores_a[[a]] <- scores_d[[a]] <- vector("list", N)
    if(x$type == "3PL") {
      scores_g[[a]] <- vector("list", N)
    } else if(x$type == "3PLu") {
      scores_u[[a]] <- vector("list", N)
    } else if(x$type == "4PL") {
      scores_g[[a]] <- scores_u[[a]] <- vector("list", N)
    }
    M_tmp <- sum(g == a)
    for(j in 1L:N) {
      ## relevant data and parameters
      pat <- dat[g == a, j] + 1
      aest_tmp <- aest[j]
      dest_tmp <- dest[j]
      gest_tmp <- gest[j]
      uest_tmp <- uest[j]
      ## px_tmp probability for a person falling into 0 or 1 of j given q
      exp_tmp <- exp(aest_tmp * X + dest_tmp)
      px_tmp <- gest_tmp + (uest_tmp - gest_tmp) * (exp_tmp / (1 + exp_tmp))
      if(any(px_tmp == 1)) {
        px_tmp[px_tmp == 1] <- 1 - .Machine$double.neg.eps
      }
      if(any(px_tmp == 0)) {
        px_tmp[px_tmp == 0] <- .Machine$double.neg.eps
      }
      px_tmp <- rbind(1 - px_tmp, px_tmp)
      ## select the matching k
      px_tmp_sel <- px_tmp[pat, ]
      ## handle NAs
      ex <- is.na(px_tmp_sel)
      if(any(ex)) {
        px_tmp_sel[which(ex, arr.ind = TRUE)] <- 1
      }
      ## update px_tmp_
      px_tmp_ <- px_tmp_ * px_tmp_sel
      ## calculate some intermediates that will be used repeatedly
      gdiffu <- gest_tmp - uest_tmp
      expone <- 1 + exp_tmp
      exponesq <- expone ^ 2
      ## evaluated 1st partial derivative of the model with respect to a
      px_da_tmp <- -(gdiffu * X * exp_tmp) / exponesq
      px_da_tmp <- rbind(-px_da_tmp, px_da_tmp)
      scores_a[[a]][[j]] <- px_da_tmp[pat, ] / px_tmp[pat, ]
      ## evaluated 1st partial derivative of the model with respect to d
      px_dd_tmp <- -(gdiffu * exp_tmp) / exponesq
      px_dd_tmp <- rbind(-px_dd_tmp, px_dd_tmp)
      scores_d[[a]][[j]] <- px_dd_tmp[pat, ] / px_tmp[pat, ]
      ## evaluated 1st partial derivative of the model with respect to g
      ## mirt internally logit transforms the guessing parameters
      ## for the scores this also has to be done to g as well as applying the
      ## logistic transformation prior to deriving the model after g
      if(x$type == "3PL" | x$type == "4PL") {
        gest_tmp_logit <- qlogis(gest_tmp)
        px_dg_tmp <- exp(gest_tmp_logit) /
          (((1 + exp(gest_tmp_logit)) ^ 2) * expone)
        px_dg_tmp <- rbind(-px_dg_tmp, px_dg_tmp)
        scores_g[[a]][[j]] <- px_dg_tmp[pat, ] / px_tmp[pat, ]
      }
      ## evaluated 1st partial derivative of the model with respect to u
      ## mirt internally logit transforms the upper asymptotes
      ## for the scores this also has to be done to u as well as applying the
      ## logistic transformation prior to deriving the model after u
      if(x$type == "3PLu" | x$type == "4PL") {
        uest_tmp_logit <- qlogis(uest_tmp)
        px_du_tmp <- exp(aest_tmp * X + dest_tmp + uest_tmp_logit) /
          (((1 + exp(uest_tmp_logit)) ^ 2) * expone)
        px_du_tmp <- rbind(-px_du_tmp, px_du_tmp)
        scores_u[[a]][[j]] <- px_du_tmp[pat, ] / px_tmp[pat, ]
      }
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
      rowSums(sc_d * pX[g == a, , drop = FALSE])
    })
    if(x$type == "3PL" | x$type == "4PL") {
      scores_g[[a]] <-
      lapply(scores_g[[a]], function(sc_g) {
        rowSums(sc_g * pX[g == a, , drop = FALSE])
      })
    }
    if(x$type == "3PLu" | x$type == "4PL") {
      scores_u[[a]] <-
      lapply(scores_u[[a]], function(sc_u) {
        rowSums(sc_u * pX[g == a, , drop = FALSE])
      })
    }
    ## combine the scores into one matrix
    text <-
    if(x$type == "4PL") {
      paste0("cbind(", paste0(
        sapply(1L:N, function(j) {
          paste0(gsub("item", j, "scores_a[[a]][[item]]"),
            ",", gsub("item", j, "scores_d[[a]][[item]]"),
            ",", gsub("item", j, "scores_g[[a]][[item]]"),
            ",", gsub("item", j, "scores_u[[a]][[item]]"))
        }), collapse = ","), ")"
      )
    } else if(x$type == "3PL") {
      paste0("cbind(", paste0(
        sapply(1L:N, function(j) {
          paste0(gsub("item", j, "scores_a[[a]][[item]]"),
            ",", gsub("item", j, "scores_d[[a]][[item]]"),
            ",", gsub("item", j, "scores_g[[a]][[item]]"))
        }), collapse = ","), ")"
      )
    } else if(x$type == "3PLu") {
       paste0("cbind(", paste0(
        sapply(1L:N, function(j) {
          paste0(gsub("item", j, "scores_a[[a]][[item]]"),
            ",", gsub("item", j, "scores_d[[a]][[item]]"),
            ",", gsub("item", j, "scores_u[[a]][[item]]"))
        }), collapse = ","), ")"
      )
    } else {
      paste0("cbind(", paste0(
        sapply(1L:N, function(j) {
          paste0(gsub("item", j, "scores_a[[a]][[item]]"),
            ",", gsub("item", j, "scores_d[[a]][[item]]"))
        }), collapse = ","), ")"
      )
    }
    scores[g == a, ] <- eval(parse(text = text))
  }
  ## restriction to raschmodel (scores are the rowsums over the restrictions)
  if(x$type == "1PL") {
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



ppl <- function(theta = NULL, a = NULL, b = NULL, g = NULL, u = NULL)
{
  ## probabilities under the model (in IRT formulation) given theta
  stopifnot(!is.null(theta) & !is.null(a) & !is.null(b) & !is.null(g) &
    !is.null(u))
  stopifnot(all.equal(mode(a), mode(b), mode(g), mode(u)))
  if(is.list(theta)) {
    return(lapply(theta, ppl, a = a, b = b, g = g, u = u))
  }
  if(is.list(a) | is.list(b | is.list(g) | is.list(u))) {
    return(mapply(ppl, a, b, MoreArgs = list(g = g, u = u, theta = theta),
      SIMPLIFY = FALSE))
  }
  a <- matrix(rep(a, length(theta)), length(theta), length(a), TRUE)
  g <- matrix(rep(g, length(theta)), length(theta), length(g), TRUE)
  u <- matrix(rep(u, length(theta)), length(theta), length(u), TRUE)
  return(g + (u - g) * plogis(a * outer(theta, b, "-")))
}



rpl <- function(theta, a = NULL, b, g = NULL, u = NULL, return_setting = TRUE)
{
  ## sample data under the model (in IRT formulation) given theta
  stopifnot(!is.null(theta) & !is.null(b))
  N <- length(b)
  if (is.null(a)) a <- rep(1, length = N)
  if (is.null(g)) g <- rep(0, length = N)
  if (is.null(u)) u <- rep(1, length = N)
  stopifnot(all.equal(mode(a), mode(b), mode(g), mode(u)))
  if(is.list(theta)) {
    return(lapply(theta, rpl, a = a, b = b, g = g, u = u,
      return_setting = return_setting))
  }
  if((is.list(a) | is.list(b) | is.list(g) | is.list(u))) {
    return(mapply(rpl, a, b, MoreArgs = list(g = g, u = u, theta = theta),
      return_setting = return_setting, SIMPLIFY = FALSE))
  }
  M <- length(theta)
  N <- length(b)
  probs <- ppl(theta = theta, a = a, b = b, g = g, u = u)
  res <- (matrix(runif(M * N), M, N) < probs) + 0
  if(return_setting) {
    return(list(a = a, b = b, g = g, u = u, theta = theta, data = res))
  } else {
    return(res)
  }
}

