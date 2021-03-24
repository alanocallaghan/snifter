test_that("fitsne works", {
    set.seed(42)
    m <- matrix(rnorm(2000), ncol=20) 
    out <- fitsne(m, random_state = NULL)
})
