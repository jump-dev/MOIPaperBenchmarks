# MOIPaperBenchmark

MOIPaperBenchmark is a repository containing scripts and data related to the
paper: _MathOptInterface: a data structure for mathematical optimization
problems_.

## Bridging experiments

To run the filesize experiment, use:
```
julia --project=. pmedian.jl
```

## MathOptFormat experiments

To run the filesize experiment, use:
```
cd MathOptFormat
julia --project=. mof_experiments.jl --conversion
```
To speed up this process, and because we are not measuring timing information,
you can use multiple threads:
```
cd MathOptFormat
export JULIA_NUM_THREADS=4; julia --project=. mof_experiments.jl --conversion
```

To run the timing experiment, use:
```
cd MathOptFormat
julia --project=. mof_experiments.jl --timing
```

Note: you must run the filesize experiment first to generate the files used in
the timing experiment.
