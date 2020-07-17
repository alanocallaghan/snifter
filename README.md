snifter
=======

R binding for
[openTSNE](https://opentsne.readthedocs.io/en/latest/index.html) using
[basilisk](https://bioconductor.org/packages/devel/bioc/html/basilisk.html).

Example usage
-------------

Basic usage of the function returns a matrix of t-SNE coordinates.

    library("snifter")
    set.seed(42)
    m <- matrix(rnorm(20000), ncol=20) 
    snifter <- fi_tsne(m[-(1:2), ], random_state = 42L)
    plot(snifter, pch = 19, xlab = "t-SNE 1", ylab = "t-SNE 2")

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

By setting `as_matrix = FALSE`, the function will return an R binding to
the underlying python object. This allows the user to embed new points
into an existing embedding.

    snifter <- fi_tsne(m[-(1:2), ], as_matrix = FALSE, random_state = 42L)
    plot(as.matrix(snifter), pch = 19, col = "black",
        xlab = "t-SNE 1", ylab = "t-SNE 2")
    points(embed(snifter, m[1:2, ]), pch = 19, col = "red")

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
