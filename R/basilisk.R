python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "python=3.9.13",
      "opentsne=0.6.2",
      "scikit-learn=1.1.2",
      "scipy=1.7.3",
      "numpy=1.22.0"
    )
)
