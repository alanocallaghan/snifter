#' Run FI-tSNE algorithm
#' 
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
        learning_rate = "auto",
        perplexity = 30,
        early_exaggeration = 250,
        early_exaggeration_iter = 12L,
        n_iter = 500L,
        exaggeration = NULL,
        dof = 1,
        theta = 0.5,
        n_interpolation_points = 3L,
        min_num_intervals = 50,
        ints_in_interval = 1,
        initialization = c("pca", "spectral", "random"),
        metric = "euclidean",
        metric_params = NULL,
        initial_momentum = 0.5,
        final_momentum = 0.8,
        max_grad_norm = NULL,
        # max_step_norm = 5,
        neighbors = c("auto", "exact", "annoy", "pynndescent", "approx"),
        negative_gradient_method = c("bh", "fft"),
        random_state = NULL,
        verbose = FALSE,
        ...
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
        verbose = verbose,
        ...
    )
    obj$fit(x)
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
    .check_args(x, ...)
    proc <- basiliskStart(python_env)
    basiliskRun(proc,
        function(...) {
            openTSNE <- reticulate::import("openTSNE")
            openTSNE$TSNE(...)
        },
        ...
    )
}

.check_args <- function(
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
        # max_step_norm,
        neighbors,
        negative_gradient_method,
        random_state,
        verbose,
        ...
    ) {
    assert_that(
        is.numeric(x),
        is.count(n_iter),
        is.count(n_components),
        is.count(abs(n_jobs)),
        is.numeric(learning_rate) | is.character(learning_rate),
        is.numeric(perplexity),
        is.numeric(early_exaggeration),
        is.count(early_exaggeration_iter),
        is.null(exaggeration) | is.numeric(exaggeration),
        is.numeric(dof),
        is.numeric(theta),
        is.count(n_interpolation_points),
        is.count(min_num_intervals),
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
        is.null(random_state) | is.count(random_state),
        is.logical(verbose)
    )
}
