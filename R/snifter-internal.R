#' @importFrom assertthat assert_that is.count
#' @importFrom basilisk BasiliskEnvironment basiliskStart
NULL

python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c("opentsne==0.4.3")
)
