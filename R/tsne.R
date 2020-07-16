tsne <- function(mat, n_iter = 500, perplexity = 30, n_jobs = 8L, ...) {
  oot <- reticulate::import("openTSNE")
  obj <- oot$TSNE(perplexity=30, n_jobs=8L, n_iter = n_iter, ...)
  res <- obj$fit(m)
}


# n_components=2,
# perplexity=30,
# learning_rate='auto',
# early_exaggeration_iter=250,
# early_exaggeration=12,
# n_iter=500,
# exaggeration=None,
# dof=1,
# theta=0.5,
# n_interpolation_points=3,
# min_num_intervals=50,
# ints_in_interval=1,
# initialization='pca',
# metric='euclidean',
# metric_params=None,
# initial_momentum=0.5,
# final_momentum=0.8,
# max_grad_norm=None,
# max_step_norm=5,
# n_jobs=1,
# affinities=None,
# neighbors='auto',
# negative_gradient_method='fft',
# callbacks=None,
# callbacks_every_iters=50,
# random_state=None,
# verbose=False
