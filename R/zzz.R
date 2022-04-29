register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1L)
  stopifnot(is.character(generic), length(generic) == 1L)
  stopifnot(is.character(class), length(class) == 1L)

  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }

  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

.onLoad <- function(libname, pkgname) {
  if(getRversion() < "3.6.0") {
    register_s3_method("sandwich", "estfun", "btmodel")
    register_s3_method("sandwich", "estfun", "raschmodel")
    register_s3_method("sandwich", "estfun", "rsmodel")
    register_s3_method("sandwich", "estfun", "pcmodel")
    register_s3_method("sandwich", "estfun", "gpcmodel")
    register_s3_method("sandwich", "estfun", "nplmodel")
    register_s3_method("sandwich", "estfun", "mptmodel")
    register_s3_method("sandwich", "bread",  "raschmodel")
    register_s3_method("sandwich", "bread",  "rsmodel")
    register_s3_method("sandwich", "bread",  "pcmodel")
    register_s3_method("sandwich", "bread",  "gpcmodel")
    register_s3_method("sandwich", "bread",  "nplmodel")
  }
  invisible()
}
