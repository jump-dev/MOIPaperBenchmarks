using PackageCompiler
using Libdl

PackageCompiler.create_sysimage(
    [:Benchmark];
    sysimage_path = "moibenchmark." * Libdl.dlext,
    precompile_execution_file = "pmedian-precompile.jl",
)
