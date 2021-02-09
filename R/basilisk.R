python_env <- BasiliskEnvironment(
    "fitsne",
    pkgname = "snifter",
    packages = c(
      "opentsne=0.5.1",
      "scikit-learn=0.24.1",
      if (basilisk.utils::isWindows()) "scipy=1.6.0" else "scipy=1.6.0",
      "numpy=1.20.0"
    )
)
