#' Run FI-tSNE algorithm
#' 
#' See \href{https://opentsne.readthedocs.io/en/latest/}{the openTSNE documentation}
#' for further details on these arguments and the general usage of this 
#' algorithm.
#' @param x Input data matrix.
#' @param n_components Number of t-SNE components to be produced.
#' @param n_jobs Integer scalar specifying the number of corest to be used.
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
#'  to fast Fourier transform, is more appropriate.
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
#'  Kobak et al. “Heavy-tailed kernels reveal a finer cluster structure in t-SNE
#'  visualisations”, 2019.
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
#' @examples
#'  set.seed(42)
#'  m <- matrix(rnorm(2000), ncol=20) 
#'  out <- fi_tsne(m, random_state = 42L)
#'  plot(out, pch = 19, xlab = "t-SNE 1", ylab = "t-SNE 2")
#' 
#'  ## openTSNE allows us to project new points into the existing
#'  ## embedding - useful for extremely large data.
#'  ## see https://opentsne.readthedocs.io/en/latest/api/index.html
#' 
#'  out_binding <- fi_tsne(m[-(1:2), ], random_state = 42L)
#'  new_points <- project(out_binding, new = m[1:2, ], old = m[-(1:2), ])
#'  plot(as.matrix(out_binding), col = "black", pch = 19,
#'      xlab = "t-SNE 1", ylab = "t-SNE 2")
#'  points(new_points, col = "red", pch = 19)
#' 
#' @rdname snifter
#' @export
snifter <- function(
        x,
        n_components = 2L,
        n_jobs = 1L,
        perplexity = 30,
        n_iter = 500L,
        initialization = c("pca", "spectral", "random"),
        neighbors = c("auto", "exact", "annoy", "pynndescent", "approx"),
        negative_gradient_method = c("bh", "fft"),
        learning_rate = "auto",
        early_exaggeration = 250,
        early_exaggeration_iter = 12L,
        exaggeration = NULL,
        dof = 1,
        theta = 0.5,
        n_interpolation_points = 3L,
        min_num_intervals = 50,
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

#' @rdname snifter
#' @export
fi_tsne <- snifter


#' Project new data into an existing t-SNE embedding object.
#' 
#' @param x t-SNE embedding created with \code{\link{fi_tsne}}.
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
#'  out_binding <- fi_tsne(m[-(1:2), ], random_state = 42L)
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
    stopifnot(inherits(x, "snifter"))
    stopifnot(ncol(new) == ncol(old))
    stopifnot(nrow(old) == nrow(x))
    new <- as.matrix(new)
    old <- as.matrix(old)
    initialization <- match.arg(initialization)

    proc <- basiliskStart(python_env)
    basiliskRun(proc,
        function(x,
                new,
                old,
                perplexity,
                initialization,
                k,
                learning_rate,
                early_exaggeration_iter,
                early_exaggeration,
                exaggeration,
                n_iter,
                initial_momentum,
                final_momentum,
                max_grad_norm,
                tolerance
            ) {
            openTSNE <- reticulate::import("openTSNE")
            affinities <- openTSNE$affinity$PerplexityBasedNN(
                old,
                perplexity = attr(x, "perplexity"),
                method = attr(x, "neighbors"),
                metric = attr(x, "metric"),
                metric_params = attr(x, "metric_params")
            )
            affinities$P - attr(x, "affinities")
            check <- all.equal(
                as.matrix(affinities$P),
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
            x_embedding$transform(
                new,
                perplexity = perplexity,
                initialization = initialization,
                k = k,
                learning_rate = learning_rate,
                early_exaggeration = early_exaggeration,
                early_exaggeration_iter = as.integer(early_exaggeration_iter),
                exaggeration = exaggeration,
                n_iter = n_iter,
                initial_momentum = initial_momentum,
                final_momentum = final_momentum,
                max_grad_norm = max_grad_norm
            )
        },
        x = x,
        new = new,
        old = old,
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

#' @export
print.snifter <- function(x, ...) {
    attributes(x) <- NULL
    print(as.matrix(x))
}

#' @importFrom reticulate py_to_r
#' @export
py_to_r.openTSNE.tsne.TSNEEmbedding <- function(x) {
    reticulate::py_to_r_wrapper(x)
    # list(
    #     affinities = x$affinities$P,
    #     mat = NextMethod()
    # )
}

.create_tsne <- function(
        x,
        ...
    ) {
    proc <- basiliskStart(python_env)
    basiliskRun(proc,
        function(x, ...) {
            openTSNE <- reticulate::import("openTSNE")
            obj <- openTSNE$TSNE(...)
            out <- obj$fit(x)
            list(
                x = reticulate:::py_ref_to_r(out),
                affinities = out$affinities$P
            )
        },
        x = x,
        ...
    )
}
.check_params <- function(x,
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
        is.integer(n_iter),
        is.integer(n_components),
        is.integer(abs(n_jobs)),
        is.numeric(learning_rate) | is.character(learning_rate),
        is.numeric(perplexity),
        is.numeric(early_exaggeration),
        is.integer(early_exaggeration_iter),
        is.null(exaggeration) | is.numeric(exaggeration),
        is.numeric(dof),
        is.numeric(theta),
        is.integer(n_interpolation_points),
        is.integer(min_num_intervals),
        is.numeric(ints_in_interval),
        is.character(initialization) | is.matrix(initialization),
        is.character(metric),
        is.null(metric_params) | is.character(metric_params),
        is.numeric(initial_momentum),
        is.numeric(final_momentum),
        is.null(max_grad_norm) | is.numeric(max_grad_norm),
        is.character(neighbors),
        is.character(negative_gradient_method),
        is.null(random_state) | is.integer(random_state),
        is.logical(verbose)
    )
}

# max_grad_norm for embedding new points
# when adding points into an existing embedding and the new points overlap 
# with the reference points, leading to large gradients. This can make points 
# “shoot off” from the embedding, causing the interpolation method to compute 
# a very large grid, and leads to worse results.
