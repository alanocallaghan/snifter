#' Run FI-tSNE algorithm
#' 
#' See \href{https://opentsne.readthedocs.io/en/latest/}{the openTSNE documentation}
#' for further details on these arguments and the general usage of this 
#' algorithm.
#' @param x Input data matrix.
#' @param simplified Logical scalar. When \code{FALSE}, the function
#' returns an object of class \code{snifter}. This contains all information
#' necessary to project new data into the embedding using \code{\link{project}}
#' If \code{TRUE}, all extra attributes will be omitted, and the return value
#' is a base matrix.
#' @param n_components Number of t-SNE components to be produced.
#' @param n_jobs Integer scalar specifying the number of cores to be used.
#' @param perplexity Numeric scalar controlling the neighborhood used
#'  when estimating the embedding.
#' @param n_iter Integer scalar specifying the number of iterations to complete.
#' @param initialization Character scalar specifying the initialization
#'  to use. "pca" may preserve global distance better than other options.
#' @param neighbors Character scalar specifying the nearest neighbour
#'  algorithm to use.
#' @param negative_gradient_method Character scalar specifying the 
#'  negative gradient approximation to use. "bh", referring to Barnes-Hut,
#'  is more appropriate for smaller data sets, while "fft" referring
#'  to fast Fourier transform, is more appropriate for larger datasets.
#' @param learning_rate Numeric scalar specifying the learning rate, or the
#'  string "auto", which uses \code{max(200, N / 12)}, where \code{N} is
#'  the number of observations.
#' @param early_exaggeration Numeric scalar specifying the exaggeration factor
#'  to use during the early exaggeration phase. Typical values range from 12 to
#'  32.
#' @param early_exaggeration_iter Integer scalar specifying the number of 
#'  iterations to run in the early exaggeration phase.
#' @param exaggeration Numeric scalar specifying the exaggeration factor to use 
#'  during the normal optimization phase. This can be used to form more densely 
#'  packed clusters and is useful for large data sets.
#' @param dof Numeric scalar specifying the degrees of freedom, as described in
#'  Kobak et al. (2019).
#' @param theta Numeric scalar, only used when negative_gradient_method="bh". 
#'  This is the trade-off parameter between speed and accuracy of the tree 
#'  approximation method. Typical values range from 0.2 to 0.8. The value 0 
#'  indicates that no approximation is to be made and produces exact results 
#'  also producing longer runtime.
#' @param n_interpolation_points Integer scalar, only used when
#'  negative_gradient_method="fft". The number of 
#'  interpolation points to use within each grid cell for interpolation based 
#'  t-SNE. It is highly recommended leaving this value at the default 3.
#' @param min_num_intervals Integer scalar, only used when
#'  negative_gradient_method="fft". The minimum number of grid cells to use, 
#'  regardless of the ints_in_interval parameter. Higher values provide more 
#'  accurate gradient estimations.
#' @param ints_in_interval Numeric scalar, only used when
#'  negative_gradient_method="fft". Indicates how large a grid cell should be
#'  e.g. a value of 3 indicates a grid side length of 3. Lower values provide 
#'  more accurate gradient estimations.
#' @param metric Character scalar specifying the metric to be used to compute 
#'  affinities between points in the original space.
#' @param metric_params Named list of additional keyword arguments for the 
#'  metric function.
#' @param initial_momentum Numeric scalar specifying the momentum to use during
#'  the early exaggeration phase.
#' @param final_momentum Numeric scalar specifying the momentum to use during 
#'  the normal optimization phase.
#' @param max_grad_norm Numeric scalar specifying the maximum gradient norm. 
#'  If the norm exceeds this value, it will be clipped.
#' @param random_state Integer scalar specifying the seed used by the random
#'  number generator.
#' @param verbose Logical scalar controlling verbosity.
#' @return A matrix of t-SNE embeddings.
#'
#' @references
#'  openTSNE: a modular Python library for t-SNE dimensionality reduction and
#'  embedding
#'  Pavlin G. Poličar, Martin Stražar, Blaž Zupan
#'  bioRxiv (2019) 731877; doi: \url{https://doi.org/10.1101/731877}
#'
#'  Fast interpolation-based t-SNE for improved visualization of single-cell
#'  RNA-seq data
#'  George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger,
#'  and Yuval Kluger
#'  Nature Methods 16, 243–245 (2019)
#'  doi: \url{https://doi.org/10.1038/s41592-018-0308-4}
#'
#'  Accelerating t-SNE using Tree-Based Algorithms
#'  Laurens van der Maaten
#'  Journal of Machine Learning Research (2014)
#'  \url{http://jmlr.org/papers/v15/vandermaaten14a.html}
#'
#'  Heavy-tailed kernels reveal a finer cluster structure in t-SNE visualisations
#'  Dmitry Kobak, George Linderman, Stefan Steinerberger, Yuval Kluger and
#'  Philipp Berens
#'  arXiv (2019)
#'  doi: \url{https://doi.org/10.1007/978-3-030-46150-8_8}.
#' @examples
#'  set.seed(42)
#'  m <- matrix(rnorm(2000), ncol=20) 
#'  out <- fitsne(m, random_state = 42L)
#'  plot(out, pch = 19, xlab = "t-SNE 1", ylab = "t-SNE 2")
#' 
#'  ## openTSNE allows us to project new points into the existing
#'  ## embedding - useful for extremely large data.
#'  ## see https://opentsne.readthedocs.io/en/latest/api/index.html
#' 
#'  out_binding <- fitsne(m[-(1:2), ], random_state = 42L)
#'  new_points <- project(out_binding, new = m[1:2, ], old = m[-(1:2), ])
#'  plot(as.matrix(out_binding), col = "black", pch = 19,
#'      xlab = "t-SNE 1", ylab = "t-SNE 2")
#'  points(new_points, col = "red", pch = 19)
#' 
#' @rdname snifter
#' @export
fitsne <- function(
        x,
        simplified = FALSE,
        n_components = 2L,
        n_jobs = 1L,
        perplexity = 30,
        n_iter = 500L,
        initialization = c("pca", "spectral", "random"),
        neighbors = c("auto", "exact", "annoy", "pynndescent", "approx"),
        negative_gradient_method = c("fft", "bh"),
        learning_rate = "auto",
        early_exaggeration = 12,
        early_exaggeration_iter = 250L,
        exaggeration = NULL,
        dof = 1,
        theta = 0.5,
        n_interpolation_points = 3L,
        min_num_intervals = 50L,
        ints_in_interval = 1,
        metric = "euclidean",
        metric_params = NULL,
        initial_momentum = 0.5,
        final_momentum = 0.8,
        max_grad_norm = NULL,
        random_state = NULL,
        verbose = FALSE
    ) {
    x <- as.matrix(x)

    initialization <- match.arg(initialization)
    neighbors <- match.arg(neighbors)
    negative_gradient_method <- match.arg(negative_gradient_method)
    random_state <- as.integer(random_state)
    early_exaggeration_iter <- as.integer(early_exaggeration_iter)
    n_jobs <- as.integer(n_jobs)
    n_components <- as.integer(n_components)
    min_num_intervals <- as.integer(min_num_intervals)
    n_iter <- as.integer(n_iter)

    out <- .create_tsne(
        x = x,
        n_iter = n_iter,
        n_components = n_components,
        n_jobs = n_jobs,
        learning_rate = learning_rate,
        perplexity = perplexity,
        early_exaggeration = early_exaggeration,
        early_exaggeration_iter = early_exaggeration_iter,
        exaggeration = exaggeration,
        dof = dof,
        theta = theta,
        n_interpolation_points = n_interpolation_points,
        min_num_intervals = min_num_intervals,
        ints_in_interval = ints_in_interval,
        initialization = initialization,
        metric = metric,
        metric_params = metric_params,
        initial_momentum = initial_momentum,
        final_momentum = final_momentum,
        max_grad_norm = max_grad_norm,
        neighbors = neighbors,
        negative_gradient_method = negative_gradient_method,
        random_state = random_state,
        verbose = verbose
    )
    if (simplified) {
        return(out$x)
    }
    structure(
        .Data = out$x,
        affinities = out$affinities,
        n_iter = n_iter,
        n_components = n_components,
        n_jobs = n_jobs,
        learning_rate = learning_rate,
        perplexity = perplexity,
        early_exaggeration = early_exaggeration,
        early_exaggeration_iter = early_exaggeration_iter,
        exaggeration = exaggeration,
        dof = dof,
        theta = theta,
        n_interpolation_points = n_interpolation_points,
        min_num_intervals = min_num_intervals,
        ints_in_interval = ints_in_interval,
        initialization = initialization,
        metric = metric,
        metric_params = metric_params,
        initial_momentum = initial_momentum,
        final_momentum = final_momentum,
        max_grad_norm = max_grad_norm,
        neighbors = neighbors,
        negative_gradient_method = negative_gradient_method,
        class = c("snifter", "matrix")
    )
}


#' @export
print.snifter <- function(x, ...) {
    attributes(x) <- NULL
    print(as.matrix(x))
}

.create_tsne <- function(
        x,
        ...
    ) {
    .checks_snifter(x, ...)
    proc <- basiliskStart(python_env)
    on.exit(basiliskStop(proc))
    out <- basiliskRun(proc, .run_tsne, x = x, ...)
    out
}
.run_tsne <- function(x, ...) {
    openTSNE <- reticulate::import("openTSNE", convert = FALSE)
    obj <- openTSNE$TSNE(...)
    out <- obj$fit(x)
    list(
        x = reticulate::py_to_r(out),
        affinities = reticulate::py_to_r(out$affinities$P)
    )
}


.checks_snifter <- function(
        x,
        n_iter,
        n_components,
        n_jobs,
        learning_rate,
        perplexity,
        early_exaggeration,
        early_exaggeration_iter,
        exaggeration,
        dof,
        theta,
        n_interpolation_points,
        min_num_intervals,
        ints_in_interval,
        initialization,
        metric,
        metric_params,
        initial_momentum,
        final_momentum,
        max_grad_norm,
        neighbors,
        negative_gradient_method,
        random_state,
        verbose
    ) {
    assert_that(
        is.numeric(x),
        length(n_iter) == 1,
        is.integer(n_iter),
        n_iter > 0,
        length(n_components) == 1,
        is.integer(n_components),
        n_components > 0,
        length(n_jobs) == 1,
        is.integer(abs(n_jobs)),
        (n_jobs %in% c(-1L, -2L)) | n_jobs > 0,
        length(learning_rate) == 1,
        (is.numeric(learning_rate) & learning_rate > 0) |
            learning_rate == "auto",
        length(perplexity) == 1,
        is.numeric(perplexity),
        perplexity > 0,
        length(early_exaggeration) == 1,
        is.numeric(early_exaggeration),
        early_exaggeration > 0,
        length(early_exaggeration_iter) == 1,
        is.integer(early_exaggeration_iter),
        early_exaggeration_iter >= 0,
        is.null(exaggeration) ||
            (is.numeric(exaggeration) &
                exaggeration > 0 &
                length(exaggeration) == 1),
        length(dof) == 1,
        is.numeric(dof),
        dof > 0,
        length(theta) == 1,
        is.numeric(theta),
        theta > 0,
        length(n_interpolation_points) == 1,
        is.integer(n_interpolation_points),
        n_interpolation_points > 0,
        length(min_num_intervals) == 1,
        is.integer(min_num_intervals),
        min_num_intervals > 0,
        length(ints_in_interval) == 1,
        is.numeric(ints_in_interval),
        min_num_intervals > 0,
        is.character(initialization) | is.matrix(initialization),
        is.character(metric),
        is.null(metric_params) |
            is.character(metric_params) |
            is.list(metric_params),
        length(initial_momentum) == 1,
        is.numeric(initial_momentum),
        initial_momentum > 0,
        length(final_momentum) == 1,
        is.numeric(final_momentum),
        final_momentum > 0,
        is.null(max_grad_norm) | is.numeric(max_grad_norm),
        is.null(random_state) | is.integer(random_state),
        is.logical(verbose)
    )
}
# max_grad_norm for embedding new points
# when adding points into an existing embedding and the new points overlap 
# with the reference points, leading to large gradients. This can make points 
# “shoot off” from the embedding, causing the interpolation method to compute 
# a very large grid, and leads to worse results.
