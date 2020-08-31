using PackageCompiler, Libdl
PackageCompiler.create_sysimage(
    [:Benchmark];
    sysimage_path="moibenchmark."*Libdl.dlext,
    precompile_execution_file="pmedian-precompile.jl",
    precompile_statements_file="statements.jl"
    # see statements.jl for more details
    )