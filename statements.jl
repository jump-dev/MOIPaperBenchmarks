#=
This file was manually genererated by executing:

    julia --project=. -Jmoibenchmark.so --trace-compile=stderr pmedian.jl

The procompile statements are written in the screen. They were copied and
pasted to this file.
=#
precompile(Tuple{Base.Fix2{typeof(Base.isequal), DataType}, Type{T} where T})
precompile(Tuple{Base.Fix2{typeof(Base.isequal), UnionAll}, Type{T} where T})
precompile(Tuple{typeof(MathOptInterface.Bridges.Variable.concrete_bridge_type), DataType, Type{MathOptInterface.Nonpositives}})
precompile(Tuple{typeof(Base.Broadcast.broadcasted), Function, Function, Array{MathOptInterface.ScalarAffineTerm{Float64}, 1}})
precompile(Tuple{typeof(Base.println), Base.TTY, String, Vararg{String, N} where N})
precompile(Tuple{typeof(Base.convert), Type{Array{Float64, 1}}, Array{Float64, 1}})
precompile(Tuple{typeof(Base.show), Base.TTY, TimerOutputs.TimerOutput})