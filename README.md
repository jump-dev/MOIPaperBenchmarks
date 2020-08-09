# MOIPaperBenchmark

MOIPaperBenchmark is a repository containing scripts and data related to the
paper: _MathOptInterface: a data structure for mathematical optimization
problems_.

## One-time setup

First, install Julia 1.5, Python 3.8, and pipenv.

Once installed, you can initialize the environments as follows:

```
$ julia --project=. -e "import Pkg; Pkg.instantiate()"
$ julia --project=MathOptFormat -e "import Pkg; Pkg.instantiate()"
$ pipenv install
```

## Bridging experiments

Run the bridging experiment as follows:

```
$ pipenv run python pmedian.py
$ julia --project=. pmedian.jl
```

Note that you must run the Python version first, since it makes a temporary file
with the CVXPY results that is used by the Julia script to build the Latex
table.

## MathOptFormat experiments

Run the filesize experiment as follows:

```
$ cd MathOptFormat
$ julia --project=. mof_experiments.jl --conversion
```

To speed up this process, and because we are not measuring timing information,
you can use multiple threads:

```
$ cd MathOptFormat
$ export JULIA_NUM_THREADS=4; julia --project=. mof_experiments.jl --conversion
```

Run the timing experiment as follows:

```
$ cd MathOptFormat
$ julia --project=. mof_experiments.jl --timing
```

Note: you must run the filesize experiment first to generate the files used in
the timing experiment.
