# MOIPaperBenchmark

MOIPaperBenchmark is a repository containing scripts and data related to the
paper: _MathOptInterface: a data structure for mathematical optimization
problems_.

## One-time setup

First, install Julia 1.5, Python 3.8, and pipenv.

Once installed, you can initialize the environments as follows:

```
$ julia --project=. -e "import Pkg; Pkg.instantiate()"
$ julia --project=. precompile.jl
$ pipenv install
```

## Run experiment

Run the bridging experiment as follows:

```
$ pipenv run python pmedian.py
$ julia --project=. -Jmoibenchmark pmedian.jl
```

Note that you must run the Python version first, since it makes a temporary file
with the CVXPY results that is used by the Julia script to build the Latex
table.
