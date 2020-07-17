python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "opentsne==0.4.3",
      "scikit-learn=0.23.1",
      "scipy=1.5.1",
      "numpy=1.19.0"
    )
)
