# Joao Bueno 17 Jun 2026
# Script to install required packages to convert Mandyoc outputs from text files to netcdf
using Pkg

necessary_packages = ["DataFrames", "NCDatasets", "Glob", "CSV", "Printf", "StatsBase"]
Pkg.add(necessary_packages)

println("All necessary packages are installed!")