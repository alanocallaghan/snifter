python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "python=3.10.6",
      "opentsne=1.0.2",
      "scikit-learn=1.6.1",
      "scipy=1.15.2",
      "numpy=2.2.6"
    )
)
