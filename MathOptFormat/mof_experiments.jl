import Dates
import DelimitedFiles
import JSON
import MathOptInterface
import Statistics

const MOI = MathOptInterface

const BENCH_DIRECTORY = joinpath(@__DIR__, "benchmark")
const MPS_DIRECTORY = joinpath(@__DIR__, "miplib_mps")
const MOF_DIRECTORY = joinpath(@__DIR__, "miplib_mof")

if !isdir(BENCH_DIRECTORY)
    mkdir(BENCH_DIRECTORY)
end

if !isdir(MPS_DIRECTORY)
    mkdir(MPS_DIRECTORY)
end

if !isdir(MOF_DIRECTORY)
    mkdir(MOF_DIRECTORY)
end

if !isfile("benchmark.zip")
    run(`curl https://miplib.zib.de/downloads/benchmark.zip -o benchmark.zip`)
end

if !isfile(joinpath(BENCH_DIRECTORY, "30n20b8.mps.gz"))
    run(`unzip benchmark.zip -d benchmark`)
end

###
### Conversion functions
###

"Append `message` to the file `log_filename` and write `message` to `stdout`."
function log(log_filename, message)
    open(log_filename, "a") do io
        println(io, message)
    end
    println(message)
end

function convert_file(filename, i, N)
    bench_filename = joinpath(BENCH_DIRECTORY, filename)
    mps_filename = joinpath(MPS_DIRECTORY, filename)
    mof_filename = joinpath(
        MOF_DIRECTORY,
        replace(filename, ".mps.gz" => ".mof.json.gz")
    )
    # if isfile(mof_filename)
    #     return "[SKIPPED] $(i)/$(N) $(filename)"
    # end
    try
        mps_model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MPS)
        MOI.read_from_file(mps_model, bench_filename)
        MOI.write_to_file(mps_model, mps_filename)

        mof_model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MOF)
        MOI.copy_to(mof_model, mps_model)
        MOI.write_to_file(mof_model, mof_filename)
    catch ex
        return "[FAIL   ] $(i)/$(N) $(filename) $(ex)"
    end
    return "[SUCCESS] $(i)/$(N) $(filename)"
end

function convert_files()
    log_filename = joinpath(
        @__DIR__, Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM-SS.log")
    )
    logfile_lock, i_lock = Base.ReentrantLock(), Base.ReentrantLock()
    files = collect(readdir(BENCH_DIRECTORY))
    N, i = length(files), 0
    Threads.@threads for filename in files
        ii = lock(() -> i += 1, i_lock)
        msg = convert_file(filename, ii, N)
        lock(() -> log(log_filename, msg), logfile_lock)
    end
    return
end

function summarize_filesize()
    bench = Dict{String, Int}(
        file => Base.stat(joinpath(BENCH_DIRECTORY, file)).size
        for file in readdir(BENCH_DIRECTORY)
    )
    mof = Dict{String, Int}(
        file => Base.stat(joinpath(MOF_DIRECTORY, file)).size
        for file in readdir(MOF_DIRECTORY)
    )
    mps = Dict{String, Int}(
        file => Base.stat(joinpath(MPS_DIRECTORY, file)).size
        for file in readdir(MPS_DIRECTORY)
    )
    ratios = Float64[]
    open(joinpath(@__DIR__, "filesizes.csv"), "w") do io
        println(io, "file, bench, mof, mps")
        for (key, mps_size) in mps
            bench_size = get(bench, key, "")
            filename = replace(key, ".mps.gz" => "")
            mof_size = get(mof, filename * ".mof.json.gz", "")
            println(io, "$(filename), $(bench_size), $(mof_size), $(mps_size)")
            if !isempty(bench_size) && !isempty(mof_size)
                push!(ratios, mof_size / bench_size)
            end
        end
    end
    bench_mu = Statistics.mean(collect(values(bench)))
    mps_mu = Statistics.mean(collect(values(mps)))
    mof_mu = Statistics.mean(collect(values(mof)))
    ratio = Statistics.mean(ratios)
    return """
    Bench average size : $(bench_mu / 1024)
    MPS average size   : $(mps_mu / 1024)
    MOF average size   : $(mof_mu / 1024)
    ratio              : $(ratio)
    """
end

###
### Timing functions
###

function time_load_save(filename, directory)
    @info "Timing $(filename)"
    model = MOI.FileFormats.Model(filename = filename)
    t = time()
    MOI.read_from_file(model, joinpath(directory, filename))
    read_time = time() - t
    t = time()
    MOI.write_to_file(model, joinpath(directory, "out_" * filename))
    write_time = time() - t
    rm(joinpath(directory, "out_" * filename))
    return read_time, write_time
end

function time_load_save(directory)
    times = Dict{String, Any}()
    files = readdir(directory)
    for filename in files
        if startswith(filename, "out_")
            continue
        end
        try
            read_time, write_time = time_load_save(filename, directory)
            times[filename] = Dict("read" => read_time, "write" => write_time)
        catch e
            @info "[FAILED] $(filename) $(e)"
        end
    end
    dir = replace(
        replace(split(directory, '/')[end], ".json" => ""), "miplib_" => ""
    )
    open("timing_$(dir).json", "w") do io
        write(io, JSON.json(times))
    end
end

function read_mof(filename)
    start = time()
    MathOptInterface.FileFormats.compressed_open(
        filename, "r", MOI.FileFormats.AutomaticCompression()
    ) do io
        data = JSON.parse(io; dicttype = Dict{String, Any})
    end
    return time() - start
end

function time_json()
    times = Dict{String, Float64}()
    for filename in readdir("miplib_mof")
        @info filename
        try
            times[filename] = read_mof(joinpath("miplib_mof", filename))
        catch e
            @info "Failed: $(e)"
        end
    end
    open("timing_julia.json", "w") do io
        write(io, JSON.json(times))
    end
end

function summarize_timing()
    mps_data = JSON.parsefile("timing_mps.json"; use_mmap = false)
    mof_data = JSON.parsefile("timing_mof.json"; use_mmap = false)
    julia_data = JSON.parsefile("timing_julia.json"; use_mmap = false)
    python_data = JSON.parsefile("timing_python.json"; use_mmap = false)
    A = Matrix{Any}(undef, length(mps_data), 7)
    for (i, file) in enumerate(sort(collect(keys(mps_data))))
        mof_file = replace(file, ".mps.gz" => ".mof.json.gz")
        mof = get(mof_data, mof_file, Dict("read" => "", "write" => ""))
        A[i, 1] = replace(file, ".mps.gz" => "")
        A[i, 2] = mps_data[file]["read"]
        A[i, 3] = mof["read"]
        A[i, 4] = mps_data[file]["write"]
        A[i, 5] = mof["write"]
        A[i, 6] = get(julia_data, mof_file, "")
        A[i, 7] = get(python_data, mof_file, "")
    end
    header = ["name" "MPS_read" "MOF_read" "MPS_write" "MOF_write" "JSON.parse" "json.parse"]
    DelimitedFiles.writedlm("timing.csv", vcat(header, A))
    mps_read = Statistics.mean(A[:, 2])
    mof_read = Statistics.mean(A[:, 3])
    read_ratio = Statistics.mean(A[:, 3] ./ A[:, 2])
    mps_write = Statistics.mean(A[:, 4])
    mof_write = Statistics.mean(A[:, 5])
    write_ratio = Statistics.mean(A[:, 5] ./ A[:, 4])
    julia_json = Statistics.mean(A[:, 6])
    python_json = Statistics.mean(A[:, 7])
    json_ratio = Statistics.mean(A[:, 6] ./ A[:, 7])
    return """
    Read MPS, MOF, MOF / MPS : $(mps_read), $(mof_read), $(read_ratio)
    Write MPS, MOF, MOF / MPS : $(mps_write), $(mof_write), $(write_ratio)
    JSON : $(julia_json), $(python_json), $(json_ratio)
    """
end

###
### Main
###

function print_help()
    println("""
    To run the conversion and summarization script, use:

        export JULIA_NUM_THREADS=4; julia mof_experiments.jl --conversion

    To time reading and writing, use

        julia mof_experiments.jl --timing
    """)
end

if length(ARGS) > 0
    if ARGS[1] == "--conversion"
        convert_files()
        @info "Summarizing..."
        msg = summarize_filesize()
        println(msg)
    elseif ARGS[1] == "--timing"
        @info "Timing MPS"
        time_load_save(MPS_DIRECTORY)
        @info "Timing MOF"
        time_load_save(MOF_DIRECTORY)
        @info "Timing Julia.JSON"
        time_json()
        @info "Timing Python.json"
        run(`python main.py`)
        @info "Summarizing..."
        summarize_timing()
    else
        print_help()
    end
end
