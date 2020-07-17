snifter
=======

R binding for
[openTSNE](https://opentsne.readthedocs.io/en/latest/index.html) using
[basilisk](https://bioconductor.org/packages/devel/bioc/html/basilisk.html).

Example usage
-------------

    devtools::load_all()
    #> Loading snifter
    set.seed(42)
    m <- matrix(rnorm(20000), ncol=20) 
    out <- fi_tsne(m[-(1:2), ], random_state = 42L)
    embed(out, m[1:2, ])
    #>            [,1]      [,2]
    #> [1,]  3.3273666 -4.083796
    #> [2,] -0.7138549  3.726564

    plot(as.matrix(out), pch = 16)

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
