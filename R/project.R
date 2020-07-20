#' Project new data into an existing t-SNE embedding object.
#'
#' @param x t-SNE embedding created with \code{\link{fitsne}}.
#' @param new New data to project into existing embedding
#' @return Numeric matrix of t-SNE co-ordinates resulting from embedding
#'  \code{new} into the t-SNE embedding \code{x}.
#' @param old Data used to create the original embedding.
#' @param perplexity Numeric scalar. Perplexity can be thought of
#'  as the continuous number of nearest neighbors, for which t-SNE
#'  will attempt to preserve distances. However, when projecting,
#'  we only consider neighbors in the existing embedding i.e.
#'  each data point is placed into the embedding, independently
#'  of other new data points.
#' @param initialization Character scalar specifying the method
#'  used to compute the initial point positions to be used in the
#'  embedding space. Can be "median", "weighted" or "random".
#'  In all cases, "median" or "weighted" should be preferred.
#' @param k Integer scalar specifying the number of nearest neighbors to
#'  consider when initially placing the point onto the embedding. This is
#'  different from "perplexity" because perplexity affects
#'  optimization while this only affects the initial point positions.
#' @param learning_rate The learning rate for t-SNE optimization.
#'  When \code{learning_rate="auto"} the appropriate learning rate
#'  is selected according to max(200, N / 12), as determined in
#'  Belkina et al. Otherwise, a numeric scalar.
#' @param early_exaggeration Numeric scalar; the exaggeration factor to use
#'  during the *early exaggeration* phase. Typical values range from 12 to
#'  32.
#' @param early_exaggeration_iter The number of iterations to run
#'  in the *early exaggeration* phase.
#' @param exaggeration The exaggeration factor to use during the
#'  normal optimization phase. This can be used to form more densely
#'  packed clusters and is useful for large data sets.
#' @param n_iter The number of iterations to run in the normal
#'  optimization regime.
#' @param initial_momentum The momentum to use during the
#'  *early exaggeration* phase.
#' @param final_momentum The momentum to use during the normal
#'  optimization phase.
#' @param max_grad_norm Maximum gradient norm. If the norm exceeds
#'  this value, it will be clipped. When adding points into an
#'  existing embedding, and the new points
#'  overlap with the reference points, this may lead to large
#'  gradients.
#'  This can make points "shoot off" from the embedding, causing the
#'  interpolation method to compute a very large grid, and leads to
#'  worse results.
#' @param tolerance Numeric scalar specifying the numeric tolerance
#'  used to ensure the affinities calculated on the old data match
#'  those of the original embedding.
#' @references
#'  Automated optimized parameters for T-distributed stochastic neighbor
#'  embedding improve visualization and analysis of large datasets.
#'  Belkina, A.C., Ciccolella, C.O., Anno, R. et al.
#'  Nature Communications 10, 5415 (2019).
#'  doi: \url{https://doi.org/10.1038/s41467-019-13055-y}
#' @examples
#'  set.seed(42)
#'  m <- matrix(rnorm(2000), ncol=20)
#'  out_binding <- fitsne(m[-(1:2), ], random_state = 42L)
#'  new_points <- project(out_binding, new = m[1:2, ], old = m[-(1:2), ])
#'  plot(as.matrix(out_binding), col = "black", pch = 19,
#'      xlab = "t-SNE 1", ylab = "t-SNE 2")
#'  points(new_points, col = "red", pch = 19)
#' @export
project <- function(
        x,
        new,
        old,
        perplexity = 5,
        initialization = c("median", "weighted", "random"),
        k = 25L,
        learning_rate = 0.1,
        early_exaggeration = 4,
        early_exaggeration_iter = 0L,
        exaggeration = 1.5,
        n_iter = 250L,
        initial_momentum = 0.5,
        final_momentum = 0.8,
        max_grad_norm = 0.25,
        tolerance = 1e-4
    ) {
    new <- as.matrix(new)
    old <- as.matrix(old)
    initialization <- match.arg(initialization)
    .project(x, new, old,
        perplexity = perplexity,
        initialization = initialization,
        k = k,
        learning_rate = learning_rate,
        early_exaggeration = early_exaggeration,
        early_exaggeration_iter = early_exaggeration_iter,
        exaggeration = exaggeration,
        n_iter = n_iter,
        initial_momentum = initial_momentum,
        final_momentum = final_momentum,
        max_grad_norm = max_grad_norm,
        tolerance = tolerance
    )
}

.project <- function(x, new, old, tolerance, ...) {
    .checks_project(
        x = x,
        new = new,
        old = old,
        tolerance = tolerance,
        ...
    )
    proc <- basiliskStart(python_env)
    on.exit(basiliskStop(proc))
    out <- basiliskRun(proc, .run_project,
        x = x, new = new, old = old, tolerance = tolerance,
        ...
    )
    out
}

.run_project <- function(x, new, old, tolerance, ...) {
    openTSNE <- reticulate::import("openTSNE", convert = FALSE)
    affinities <- openTSNE$affinity$PerplexityBasedNN(
        old,
        perplexity = attr(x, "perplexity"),
        method = attr(x, "neighbors"),
        metric = attr(x, "metric"),
        metric_params = attr(x, "metric_params")
    )
    check <- all.equal(
        as.matrix(reticulate::py_to_r(affinities$P)),
        as.matrix(attr(x, "affinities")),
        tolerance = tolerance
    )
    if (!check) {
        stop("Affinities differ from original affinities (tolerance =", tolerance)
    }
    args <- list(embedding = x, affinities = affinities)
    attr <- c(
        "negative_gradient_method",
        "dof",
        "theta",
        "n_interpolation_points",
        "min_num_intervals"
    )
    attr_args <- lapply(attr, function(parameter) attr(x, parameter))
    names(attr_args) <- attr
    args <- c(args, attr_args)
    x_embedding <- do.call(openTSNE$TSNEEmbedding, args)
    reticulate::py_to_r(x_embedding$transform(new, ...))
}

.checks_project <- function(
        x,
        new,
        old,
        perplexity,
        initialization,
        k,
        learning_rate,
        early_exaggeration,
        early_exaggeration_iter,
        exaggeration,
        n_iter,
        initial_momentum,
        final_momentum,
        max_grad_norm,
        tolerance
    ) {
    assert_that(
        is.numeric(x),
        inherits(x, "snifter"),
        ncol(new) == ncol(old),
        nrow(old) == nrow(x),
        length(perplexity) == 1,
        is.numeric(perplexity),
        perplexity > 0,
        length(k) == 1,
        is.integer(k),
        k > 0,
        length(learning_rate) == 1,
        is.numeric(learning_rate),
        learning_rate > 0,
        length(early_exaggeration) == 1,
        is.numeric(early_exaggeration),
        early_exaggeration > 0,
        length(early_exaggeration_iter) == 1,
        is.integer(early_exaggeration_iter),
        early_exaggeration_iter >= 0,
        length(n_iter) == 1,
        is.integer(n_iter),
        n_iter > 0,
        length(initial_momentum) == 1,
        is.numeric(initial_momentum),
        initial_momentum > 0,
        length(final_momentum) == 1,
        is.numeric(final_momentum),
        final_momentum > 0,
        length(max_grad_norm) == 1,
        is.numeric(max_grad_norm),
        max_grad_norm > 0,
        length(tolerance) == 1,
        is.numeric(tolerance),
        tolerance > 0
    )
}

