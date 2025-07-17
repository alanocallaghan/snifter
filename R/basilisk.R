python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "python=3.12.10",
      "opentsne=1.0.2",
      "scikit-learn=1.7.0",
      "scipy=1.16.0",
      "numpy=2.3.1"
    )
)
