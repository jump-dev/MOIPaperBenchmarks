using PackageCompiler
PackageCompiler.create_sysimage(
    [:MathOptInterface, :GLPK, :JSON, :SCS, :TimerOutputs];
    sysimage_path="moibenchmark.so",
    precompile_execution_file="pmedian-precompile.jl")

exit(0)