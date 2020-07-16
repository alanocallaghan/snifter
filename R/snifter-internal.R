NULL

#' @importFrom basilisk BasiliskEnvironment
python_env <- BasiliskEnvironment("fitsne", pkgname="snifter",
    packages=c("opentsne==0.4.3"))
