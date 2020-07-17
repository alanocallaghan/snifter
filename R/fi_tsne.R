#' Run FI-tSNE algorithm
#' 
#' See [the openTSNE documentation](https://opentsne.readthedocs.io/en/latest/)
#' for further details on these arguments and the general usage of this 
#' algorithm.
#' @param x Input data matrix.
#' @param n_components Number of t-SNE components to be produced.
#' @param n_jobs Integer scalar specifying the number of corest to be used.
#' @param perplexity Numeric scalar controlling the neighborhood used
#' when estimating the embedding.
#' @param n_iter Integer scalar specifying the number of iterations to complete.
#' @param initialization Character scalar specifying the initialization
#' to use. "pca" may preserve global distance better than other options.
#' @param neighbors Character scalar specifying the nearest neighbour
#' algorithm to use.
#' @param negative_gradient_method Character scalar specifying the 
#' negative gradient approximation to use. "bh", referring to Barnes-Hut,
#' is more appropriate for smaller data sets, while "fft" referring
#' to fast Fourier transform, is more appropriate.
#' @param learning_rate Numeric scalar specifying the learning rate, or the
#' string "auto", which uses \code{max(200, N / 12)}, where \code{N} is
#' the number of observations.
#' @param early_exaggeration Numeric scalar specifying the exaggeration factor
#' to use during the early exaggeration phase. Typical values range from 12 to
#' 32.
#' @param early_exaggeration_iter Integer scalar specifying the number of 
#' iterations to run in the early exaggeration phase.
#' @param exaggeration Numeric scalar specifying the exaggeration factor to use 
#' during the normal optimization phase. This can be used to form more densely 
#' packed clusters and is useful for large data sets.
#' @param dof Numeric scalar specifying the degrees of freedom, as described in
#' Kobak et al. “Heavy-tailed kernels reveal a finer cluster structure in t-SNE
#' visualisations”, 2019.
#' @param theta Numeric scalar, only used when negative_gradient_method="bh". 
#' This is the trade-off parameter between speed and accuracy of the tree 
#' approximation method. Typical values range from 0.2 to 0.8. The value 0 
#' indicates that no approximation is to be made and produces exact results 
#' also producing longer runtime.
#' @param n_interpolation_points Integer scalar, only used when
#' negative_gradient_method="fft". The number of 
#' interpolation points to use within each grid cell for interpolation based 
#' t-SNE. It is highly recommended leaving this value at the default 3.
#' @param min_num_intervals Integer scalar, only used when
#' negative_gradient_method="fft". The minimum number of grid cells to use, 
#' regardless of the ints_in_interval parameter. Higher values provide more 
#' accurate gradient estimations.
#' @param ints_in_interval Numeric scalar, only used when
#' negative_gradient_method="fft". Indicates how large a grid cell should be
#' e.g. a value of 3 indicates a grid side length of 3. Lower values provide 
#' more accurate gradient estimations.
#' @param metric Character scalar specifying the metric to be used to compute 
#' affinities between points in the original space.
#' @param metric_params Named list of additional keyword arguments for the 
#' metric function.
#' @param initial_momentum Numeric scalar specifying the momentum to use during
#' the early exaggeration phase.
#' @param final_momentum Numeric scalar specifying the momentum to use during 
#' the normal optimization phase.
#' @param max_grad_norm Numeric scalar specifying the maximum gradient norm. 
#' If the norm exceeds this value, it will be clipped.
#' This is most beneficial 
#' @param random_state Integer scalar specifying the seed used by the random
#' number generator.
#' @param verbose Logical scalar controlling verbosity.
#' @param ... Unused.
#' @return A matrix of t-SNE embeddings.
fi_tsne <- function(x, ...) {
    UseMethod("fi_tsne")
}

#' @export
#' @rdname fi_tsne
fi_tsne.matrix <- function(
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
        # max_step_norm = 5,
        random_state = NULL,
        verbose = FALSE
    ) {

    initialization <- match.arg(initialization)
    neighbors <- match.arg(neighbors)
    negative_gradient_method <- match.arg(negative_gradient_method)
    obj <- .create_tsne(
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
        # max_step_norm = max_step_norm,
        neighbors = neighbors,
        negative_gradient_method = negative_gradient_method,
        random_state = random_state,
        verbose = verbose
    )
    obj$fit(x)
}

#' @export
fi_tsne.data.frame <- function(x, ...) {
    fi_tsne(as.matrix(x), ...)
}

#' Embed new data in an existing t-SNE embedding object.
#' 
#' @param x t-SNE embedding created with \code{\link{fi_tsne_embedding}}.
#' @param new New data to project into existing embedding
#' @export
embed <- function(x, ...) {
    UseMethod("embed")
}

#' @export
#' @rdname embed
embed.openTSNE.tsne.TSNEEmbedding <- function(x, new) {
    proc <- basiliskStart(python_env)
    basiliskRun(proc,
        function(x, new) {
            x$transform(new)
        },
        x = x,
        new = new
    )
}

#' @export
print.openTSNE.tsne.TSNEEmbedding <- function(x) {
    cat("reticulate binding to openTSNE.tsne.TSNEEmbedding.\n")
    cat("Use as.matrix() for values.\n")
}

#' @export
as.matrix.openTSNE.tsne.TSNEEmbedding <- function(x) {
    reticulate:::py_to_r.default(x)
}


#' @importFrom reticulate py_to_r
#' @export
py_to_r.openTSNE.tsne.TSNEEmbedding <- function(x) {
    reticulate::py_to_r_wrapper(x)
}

.create_tsne <- function(
        x,
        ...
    ) {
    proc <- basiliskStart(python_env)
    basiliskRun(proc,
        function(...) {
            openTSNE <- reticulate::import("openTSNE")
            openTSNE$TSNE(...)
        },
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
        # is.numeric(max_step_norm),
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
