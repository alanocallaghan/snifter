python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "python=3.10.6",
      "opentsne=0.6.2",
      "scikit-learn=1.1.2",
      if (basilisk.utils::isWindows()) "scipy=1.7.3" else "scipy=1.7.3",
      "numpy=1.22.0"
    )
)
