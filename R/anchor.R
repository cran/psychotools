###################
## choose anchor ##
###################

## generic function anchor
anchor <- function(object, ...)
{  
  UseMethod("anchor")         
}

## formula interface for generic function anchor
anchor.formula <- function(formula, data = NULL, subset = NULL,
  na.action = NULL, weights = NULL, model = raschmodel, ...)
{
  ## build model frame
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- formula <- Formula::as.Formula(formula)
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract item responses, binary factor, and weights
  y <- model.matrix(~ 0 + ., Formula::model.part(formula, mf, lhs = 1L))
  x <- Formula::model.part(formula, mf, rhs = 1L)[[1L]]
  w <- model.weights(mf)
  
  ## fit separate Rasch models
  x0 <- sort(unique(x))
  if(length(x0) != 2L) stop("explanatory variable needs to be a binary factor")
  obj1 <- model(y[x == x0[1L], , drop = FALSE], weights = w[x == x0[1L]])
  obj2 <- model(y[x == x0[2L], , drop = FALSE], weights = w[x == x0[2L]])

  ## call default method
  anchor.default(obj1, obj2, ...)
}

## default method for generic function anchor
anchor.default <- function(object, object2, class = c("constant", "forward"),
  select = NULL, length = NULL, range = c(0.1, 0.8), ...)
{
  ## check if items are equal and finite in both groups
  cf1 <- coef(itempar(object))
  cf2 <- coef(itempar(object2))
  stopifnot(all.equal(names(cf1), names(cf2)))
  if (any(is.infinite(cf1)) | any(is.infinite(cf2))) {
    infitems <- unique(c(names(cf1)[is.infinite(cf1)], names(cf2)[is.infinite(cf2)]))
    stop("Infinite item parameter estimates are not allowed: (", paste(infitems, sep = "", collapse = ", "), ").")
  }

  ## anchor method
  class <- match.arg(class)
  if(is.null(select)) {
    select <- if(class == "constant" & (is.null(length) || length == 1L)) {
      "Gini"
    } else if(class == "constant") {
      "MPT"
    } else {
      "MTT"
    }
  }
  select <- match.arg(select, c("MTT", "MPT", "MT", "MP", "NST", "AO", "AOP", "Gini", "CLF", "GiniT", "CLFT"))

  ## default anchor length:
  ## - constant: 1 for alignment anchors, 4 for "classic" test-based anchors
  ## - forward: unrestricted (up to number of items)
  if(is.null(length)) length <- if(class == "constant") {
    if(select %in% c("Gini", "CLF", "GiniT", "CLFT")) 1 else 4
  } else {
    length(cf1)
  }

  ## warn about alignment selection strategy with length > 1 or forward class
  if(select %in% c("Gini", "CLF", "GiniT", "CLFT")) {
    if(class == "forward") stop("Anchor forward selection cannot be combined with Gini or CLF.")
    if(length > 1) warning("Selecting alignment anchors (Gini or CLF) with length > 1 might not yield optimal results.")
    ## FIXME: when alignanchor() becomes available:
    ## if(length > 1) warning("For selecting alignment anchors (Gini or CLF) with length > 1, please use the alignanchor() function.")
  }

  ## carry out auxiliary DIF tests
  aux_tests <- auxtests(obj1 = object, obj2 = object2, select)$aux_tests
  selection_results <- anchorselect(obj1 = object, obj2 = object2, aux_tests, select)

  anchor_items <- anchorclass(obj1 = object, obj2 = object2, ranking_order = selection_results$ranking_order, 
    class, length = length, range = range)$anchor_items
  results <- list(anchor_items = anchor_items, ranking_order = selection_results$ranking_order,
    criteria = selection_results$criteria, class = class, select = select, length = length, range = range, item_labels = names(cf1))
  class(results) <- "anchor"
  results
}


## methods for class 'anchor'
print.anchor <- function(x, ...)
{
  method <- switch(x$class,
    "constant" = sprintf("Anchor selection with %s criterion and constant length %s", x$select, x$length),
    "forward" = sprintf("Anchor forward selection with %s criterion and maximum length %s", x$select, x$length)
  )
  
  cat(method, "\n")
  cat(if(length(x$anchor_items) > 1L) "Anchor items:" else "Anchor item:",
    strwrap(paste(x$item_labels[x$anchor_items], collapse = ", ")),
    "\n")
  
  invisible(x)  
}

summary.anchor <- function(object, ...)
{
  class(object) <- c("summary.anchor", class(object))
  return(object)
}

print.summary.anchor <- function(x, ...)
{
  method <- switch(x$class,
    "constant" = sprintf("Anchor selection with %s criterion and constant length %s", x$select, x$length),
    "forward" = sprintf("Anchor forward selection with %s criterion and maximum length %s", x$select, x$length)
  )
  
  cat(method, "\n")

  cat(if(length(x$anchor_items) > 1L) "\nAnchor items:" else "\nAnchor item:", "\n")
  cat(strwrap(paste(x$item_labels[x$anchor_items], collapse = ", ")),
    "\n")

  cat("\nRanking order:\n",
    strwrap(paste(x$ranking_order, collapse = ", ")),
    "\n\n")
  cat("Criterion values (not sorted):\n")
  print(round(setNames(x$criteria, x$item_labels), digits = 2 + (x$select %in% c("Gini", "GiniT"))))

  invisible(x)  
}


#############################
## pairwise anchored tests ##
#############################

## generic function anchortest
anchortest <- function(object, ...)
{
  UseMethod("anchortest")
}

## formula interface for generic function anchortest
anchortest.formula <- function(formula, data = NULL, subset = NULL,
  na.action = NULL, weights = NULL, model = raschmodel, ...)
{
  ## build model frame
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- formula <- Formula::as.Formula(formula)
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract item responses, binary factor, and weights
  y <- model.matrix(~ 0 + ., Formula::model.part(formula, mf, lhs = 1L))
  x <- Formula::model.part(formula, mf, rhs = 1L)[[1L]]
  w <- model.weights(mf)
  
  ## fit separate Rasch models
  x0 <- sort(unique(x))
  if(length(x0) != 2L) stop("explanatory variable needs to be a binary factor")
  obj1 <- model(y[x == x0[1L], , drop = FALSE], weights = w[x == x0[1L]])
  obj2 <- model(y[x == x0[2L], , drop = FALSE], weights = w[x == x0[2L]])
  
  ## call default method
  anchortest.default(obj1, obj2, ...)
}

## default method for generic function anchor
anchortest.default <- function(object, object2,  class = c("constant", "forward", "all-other", "fixed"),
  select = NULL, test = TRUE, adjust = "none", length = NULL, range = c(0.1, 0.8), ...)
{
  ## check input:
  ## - default: constant-MPT  
  ## - if class="forward" -> select="MTT"
  ## - if select="numeric" -> class="fixed"
  if(is.null(select) || is.character(select)) {
    class <- match.arg(class)
    if(is.null(select)) {
      select <- if(class == "constant" & (is.null(length) || length == 1L)) {
        "Gini"
      } else if(class == "constant") {
        "MPT"
      } else {
        "MTT"
      }    
    }
    select <- match.arg(select, c("MTT", "MPT", "MT", "MP", "NST", "AO", "AOP", "Gini", "CLF", "GiniT", "CLFT"))
  } else {
    if(!missing(class) && match.arg(class) != "fixed") warning("if 'select' specifies the anchor items directly, 'class' needs to be 'fixed'")
    class <- "fixed"
  }
  adjust <- match.arg(adjust, c("single-step", "Shaffer", "Westfall", "free", p.adjust.methods))
  
  ## here: add all-other
  if (class == "all-other") {
    anchorres <- list(
      anchor_items = NA,
      ranking_order = NA,
      criteria = NA
    )
    final_tests <- allothertests(obj1 = object, obj2 = object2, adjust = adjust) 
  } else if (class == "fixed") {
    if(!is.numeric(select)) stop("if 'class' is 'fixed', 'select' needs to specify the anchor items directly")
    anchorres <- list(
      anchor_items = select,
      ranking_order = NA,
      criteria = NA
    )
    final_tests <- diftests(obj1 = object, obj2 = object2, anchor_items = select, adjust = adjust)
  } else {
    anchorres <- anchor(object = object, object2 = object2, class = class, select = select, length = length, range = range)
    final_tests <- diftests(obj1 = object, obj2 = object2, anchor_items = anchorres$anchor_items, adjust = adjust)
  }

  results <- list(anchor_items = anchorres$anchor_items, ranking_order = anchorres$ranking_order,
                  criteria = anchorres$criteria, anchored_item_parameters =  final_tests$itempars,
                  anchored_covariances = final_tests$vcovs, final_tests = if(test) final_tests$test else NULL,
                  item_labels = colnames(object$data))

  class(results) <- "anchortest"
  results
}

## methods for class anchortest
print.anchortest <- function(x, ...)
{
  cat("Anchor items:\n")
  writeLines(if(isTRUE(all.equal(x$anchor_items, NA))) "NA" else strwrap(paste(x$item_labels[x$anchor_items], collapse = ", ")))
  cat("\nFinal DIF tests:\n")
  print(x$final_tests, ...)
  invisible(x)  
}

summary.anchortest <- function(object, ...)
{
  class(object) <- c("summary.anchortest", class(object))
  return(object)
}

print.summary.anchortest <- function(x, ...)
{
  cat("Anchor items:\n")
  writeLines(if(isTRUE(all.equal(x$anchor_items, NA))) "NA" else strwrap(paste(x$item_labels[x$anchor_items], collapse = ", ")))
  cat("\nAnchored item parameters:\n")
  print(x$anchored_item_parameters, ...)
  cat("\nRanking order:\n")
  print(x$ranking_order)  
  cat("\nCriterion values (not sorted):\n")
  print(setNames(x$criteria, if(isTRUE(all.equal(x$criteria, NA))) NULL else x$item_labels))
  cat("\nFinal DIF tests:\n")
  print(x$final_tests, ...)
  invisible(x)  
}

plot.anchortest <- function(x, main = NULL, ...) {
  if(is.null(main)) main <- sprintf("Anchor item%s: %s",
    if(length(x$anchor_items) > 1L) "s" else "",
    paste(x$item_labels[x$anchor_items], collapse = ", "))
  plot(x$final_tests, main = main, ...)
}


###############
## utilities ##
###############

## function auxtests
auxtests <- function(obj1, obj2, select)
{
  k <- length(coef(itempar(obj1)))
  aux_tests <- matrix(NA, ncol = k, nrow = k)

  if (select %in% c("AO", "AOP")) {
    diag(aux_tests) <- allothertests(obj1 = obj1, obj2 = obj2, adjust = "none")$test$test$tstat
  } else if(select %in% c("Gini", "CLF")) {
    d <- coef(itempar(obj1)) - coef(itempar(obj2))
    aux_tests[] <- outer(d, d, "-")
  } else if (select %in% c("MTT", "MT", "GiniT", "CLFT")) {
    for(i in 1:k) {
      tempdt <- diftests(obj1, obj2, anchor_items = i, adjust = "none")$test     # contains test stats using single anchors
      aux_tests[i,-i] <- tempdt$test$tstat                                       # contains test statistics 
    }
  } else {
    for(i in 1:k){
      tempdt <- diftests(obj1, obj2, anchor_items = i, adjust = "none")$test # contains test stats using single anchors
      aux_tests[i,-i] <- tempdt$test$pvalue                                  # contains p-values
    }
  }

  return(list(aux_tests = aux_tests))

}

## function anchorselect
anchorselect <- function(obj1, obj2, aux_tests, select)
{  
  k <- length(coef(itempar(obj1)))
  mysort <- function(sortvar){ ## sort function, not first item fulfilling minimum, but randomly chosen item fulfilling minimum 
    sorttemp <- sort(sortvar, index.return = TRUE)
    sortlist <- unlist(tapply(sorttemp$ix,as.factor(sorttemp$x), FUN = function(x){
      if(length(x) != 1){
        return(as.integer(sample(x=x, size=length(x), replace = FALSE, prob = NULL)))
      } else {
        return(x)
      }}))
    names(sortlist) <- NULL
    return(sortlist)
  }

  ## inequality functions for alignment
  gini <- function(d) 2 * sum(1:k * sort(d))/(k * sum(d)) - (k + 1)/k
  clf <- function(d) sum(sqrt(d))

  switch(select,
         "MTT" = {                      # mean test statistic treshold after Kopf et al. (2013b)
           threshold <- sort(abs(colMeans(aux_tests, na.rm=TRUE)))[ceiling(0.5*k)]   # threshold definition   
           criteria  <- apply(aux_tests,2,function(x) sum(abs(x) > threshold, na.rm=TRUE))
         },

         "MPT" = {                      # mean p-value treshold after Kopf et al. (2013b)
           threshold <- sort(colMeans(aux_tests,na.rm=TRUE), decreasing = TRUE)[ceiling(0.5*k)]  # threshold definition               
           criteria  <- apply(aux_tests,2,function(x) -sum(x > threshold, na.rm=TRUE))
         },
         
         "MT" = {                       # mean test statistic after Shih and Wang (2009)
           criteria <- colMeans(abs(aux_tests), na.rm=TRUE)
         },
         
         "MP" = {                                       # mean p-value after Kopf et al. (2013b)
           criteria <- -colMeans(aux_tests, na.rm=TRUE)
         },
         
         "NST" = {                      # number of significant after Wang (2004) modified by Kopf et al. (2013a)
           criteria <- colSums(aux_tests<=0.05, na.rm=TRUE)
         },
         
         "AO" = {                       # suggestion by Woods (2009)
           criteria <- abs(diag(aux_tests)) 
           warning("The all-other selection requires strong prior knowledge that DIF is balanced.")
         },
         
         "AOP" = {                      # suggestion by Wang (2012)
           AOPtemp <- 1:k
           AOPtempnew <- which(!(abs(diag(aux_tests)) > 1.64)) ### xxx check
           
           if(length(AOPtempnew)==0){
             criteria <- abs(diag(aux_tests))
             warning("Caution: AOP-selection replaced by AO-selection (no purification necessary)!")
           } 
           
           if(all(AOPtemp %in% AOPtempnew) & all(AOPtemp %in% AOPtempnew)){            ### xxx simplify no if, only while
             dttemp <- diftests(obj1, obj2, anchor_items = AOPtempnew, adjust = "none")$test$test    
           }
           
           while(!(all(AOPtemp %in% AOPtempnew) & all(AOPtemp %in% AOPtempnew))){
             AOPtemp <- AOPtempnew
             tempdt <- diftests(obj1, obj2, anchor_items = AOPtempnew, adjust = "none")$test$test
             trep <- rep(FALSE, times=k)
             trep[-AOPtempnew[1]] <- tempdt$pvalue<=0.05                 
             AOPtempnew <- AOPtempnew[!trep[AOPtempnew]]
           }
           
           criteria <- rep(Inf, times=k)
           criteria[-AOPtempnew[1]] <- abs(tempdt$tstat)     
           
           warning("The all-other purified selection requires strong prior knowledge that DIF is balanced.")        
         },
	 
         "Gini" = {                       # alignment as suggested by Strobl et al. (2020)
           criteria <- -apply(abs(aux_tests), 2, gini)
         },

         "CLF" = {                       # alignment adaptation by Strobl et al. (2020) based on Asparouhov & Muthen (2014).
           criteria <- apply(abs(aux_tests), 2, clf)
         },
	 
         "GiniT" = {                       # modified version of Gini with t statistics rather than coefficient differences
	   diag(aux_tests) <- 0
           criteria <- -apply(abs(aux_tests), 2, gini)
         },

         "CLFT" = {                       # modified version of CLF with t statistics rather than coefficient differences
	   diag(aux_tests) <- 0
           criteria <- apply(abs(aux_tests), 2, clf)
         })

  ranking_order <- mysort(criteria)
  ranking_order <- as.integer(ranking_order)
  return(list(ranking_order = ranking_order, criteria = criteria))

}


## function anchorclass
anchorclass <- function(obj1, obj2, ranking_order, class, length, range)
{
  k <- length(coef(itempar(obj1)))
  switch(class,
         "constant" = {                 # constant anchor class, predefined number of anchor items, see, e.g. Wang (2004)
           if((abs(length - round(length)) > .Machine$double.eps^0.5) | length > k | length < 1){
             length <- 4
             warning("length argument was not correctly specified and set to default 4.")
           }
           anchor_items <- ranking_order[1:length]
         },
         "forward" = {                  # iterative forward class, see Kopf et al. (2013a)
           if((abs(length - round(length)) > .Machine$double.eps^0.5) | length > k | length < 1){
             length <- k
             warning("length argument was not correctly specified and set to default k.")
           }
           if(range[1]<0 | range[2]>1 | range[1]>range[2]){
             range <- c(0.1, 0.8)
             warning("range argument was not correctly specified and set to default c(0.1, 0.8).")				 		
           }
           if(range[1] == 0){
             indices <- 1:k 
           }else{
             indices <- c(ceiling(range[1]*k+1):k,1:(ceiling(range[1]*k)))
           }            
           index <- 0           
           anchor_cand <- ranking_order[indices]
           while(index == 0 || length(res_anchor) < min(range[2] * sum(diftesttemp$test$pvalue>0.05, na.rm=TRUE),ifelse(is.null(length),k,length))) {
             index <- index + 1             
             res_anchor  <- anchor_cand[1:index]             
             diftesttemp <- diftests(obj1, obj2, anchor_items = res_anchor, adjust = "none")$test
           }  
           anchor_items <- res_anchor
         })

  return(list(anchor_items = anchor_items))  
}


## diftest: pairwise comparison of item parameters of 
##   - two raschmodel objects (obj1, obj2)
##   - with a certain anchoring
##   - using traditional covariances (i.e., losing one
##     aliased coefficient) 
##   - the Wald tests can either be adjusted for multiple testing
##     or report just marginal Wald tests

diftests <- function(obj1, obj2, anchor_items, adjust){

  ## check for multcomp
  stopifnot(requireNamespace("multcomp"))
    
  ip1 <- itempar(obj1, ref = anchor_items)
  ip2 <- itempar(obj2, ref = anchor_items)   
	
  alias <- anchor_items[1]
  cf <- c(coef(ip1)[-alias],coef(ip2)[-alias])
	
  k <- length(cf)/2
  vc <- matrix(0, 2 * k, 2 * k)
  vc[1:k, 1:k] <- vcov(ip1)[-alias, -alias]
  vc[-(1:k), -(1:k)] <- vcov(ip2)[-alias, -alias]   
	
  names(cf) <- colnames(vc) <- rownames(vc) <- c(paste(colnames(obj1$data)[-alias], 1, sep = "_"), paste(colnames(obj2$data)[-alias], 2, sep = "_"))
	
  ## collect in trivial "model object"
  mod <- list(coefficients = cf, vcov = vc)
  class(mod) <- "raschmodel"
	
  ## new pairwise contrasts
  contr <- cbind(diag(k), -diag(k))
  colnames(contr) <- names(cf)
  rownames(contr) <- colnames(obj1$data)[-alias]
	
  ## test employed
  test <- if(adjust == "none") multcomp::univariate() else multcomp::adjusted(type = adjust)
  ftest <- summary(multcomp::glht(mod, linfct = contr), test = test)
  return(list(test = ftest, itempars = cf, vcovs = vc))

}


## allothertests: pairwise all-other comparison of item parameters of 
##   - two raschmodel objects (obj1, obj2)
##   - with all-other as anchor
##   - using traditional covariances no coefficients set aliased

allothertests <- function(obj1, obj2, adjust) {
	
  ## check for multcomp
  stopifnot(requireNamespace("multcomp"))

  ## extract estimates and covariances
  ip1 <- itempar(obj1, ref = 1)
  ip2 <- itempar(obj2, ref = 1)
  vc1 <- vcov(ip1)
  vc2 <- vcov(ip2)
  k <- length(ip1)
  
  cf <- c(ip1, ip2)
  vc <- rbind(cbind(vc1, matrix(0, nrow = k, ncol = k)),
              cbind(matrix(0, nrow = k, ncol = k), vc2))
  
  names(cf) <- colnames(vc) <- rownames(vc) <- c(
      paste(colnames(obj1$data), 1, sep = "_"),
      paste(colnames(obj2$data), 2, sep = "_"))
  
  contrnew <- matrix(data = - 1/(k-1), nrow = k, ncol = k)
  diag(contrnew) <- 1
  
  contr <- rbind(cbind(contrnew, matrix(0, nrow = k, ncol = k)),
                 cbind(matrix(0, nrow = k, ncol = k), contrnew))
  ipreturn <- c(contr %*% cf)
  
  contrnew <- cbind(contrnew, -contrnew)
  colnames(contrnew) <- names(cf)
  rownames(contrnew) <- colnames(obj1$data)
  
  ## collect in trivial "model object"	
  mod <- list(coefficients = cf, vcov = vc)
  class(mod) <- "raschmodel"
  
  ## test employed
  test <- if(adjust == "none") multcomp::univariate() else multcomp::adjusted(type = adjust)
  ftest <- summary(multcomp::glht(mod, linfct = contrnew), test = test)
  return(list(test = ftest, itempars = ipreturn, vcovs = NA))

}
