module Benchmark

import MathOptInterface
import SCS
import GLPK
import JSON
import Printf
import Random

using SparseArrays
using TimerOutputs

const MOI = MathOptInterface

struct PMedianData
    num_facilities::Int
    num_customers::Int
    num_locations::Int
    customer_locations::Vector{Float64}
end

# This is the LP relaxation.
function generate_moi_problem(model, data::PMedianData)
    NL = data.num_locations
    NC = data.num_customers

    ###
    ### 0 <= facility_variables <= 1
    ###

    facility_variables = MOI.add_variables(model, NL)

    for v in facility_variables
        MOI.add_constraint(model, MOI.SingleVariable(v), MOI.Interval(0.0, 1.0))
    end

    ###
    ### assignment_variables >= 0
    ###

    assignment_variables = reshape(MOI.add_variables(model, NC * NL), NC, NL)
    for v in assignment_variables
        MOI.add_constraint(model, MOI.SingleVariable(v), MOI.GreaterThan(0.0))
        # "Less than 1.0" constraint is redundant.
    end

    ###
    ### Objective function
    ###

    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            [
                MOI.ScalarAffineTerm(
                    abs(data.customer_locations[i] - j),
                    assignment_variables[i, j]
                )
                for i in 1:NC for j in 1:NL
            ],
            0.0,
        ),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    ###
    ### assignment_variables[i, j] <= facility_variables[j]
    ###

    for i in 1:NC, j in 1:NL
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction(
                [
                    MOI.ScalarAffineTerm(1.0, assignment_variables[i, j]),
                    MOI.ScalarAffineTerm(-1.0, facility_variables[j])
                ],
                0.0,
            ),
            MOI.LessThan(0.0),
        )
    end

    ###
    ### sum_j assignment_variables[i, j] = 1
    ###

    for i in 1:NC
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction(
                [
                    MOI.ScalarAffineTerm(1.0, assignment_variables[i, j])
                    for j in 1:NL
                ],
                0.0,
            ),
            MOI.EqualTo(1.0),
        )
    end

    ###
    ### sum_j facility_variables[j] == num_facilities
    ###

    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.(1.0, facility_variables),
            0.0,
        ),
        MOI.EqualTo{Float64}(data.num_facilities),
    )

    return assignment_variables, facility_variables
end

function generate_moi_problem_vector(model, data::PMedianData)
    NL = data.num_locations
    NC = data.num_customers

    ###
    ### facility_variables    ∈ R₊
    ### facility_variables -1 ∈ R₋
    ###

    facility_variables = MOI.add_variables(model, NL)
    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(
            [
                MOI.VectorAffineTerm(i, MOI.ScalarAffineTerm(1.0, x))
                for (i, x) in enumerate(facility_variables)
            ],
            zeros(NL),
        ),
        MOI.Nonnegatives(NL),
    )

    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(
            [
                MOI.VectorAffineTerm(i, MOI.ScalarAffineTerm(1.0, x))
                for (i, x) in enumerate(facility_variables)
            ],
            -ones(NL),
        ),
        MOI.Nonpositives(NL),
    )

    ###
    ### assignment_variables ∈ R₊
    ###

    assignment_variables = reshape(MOI.add_variables(model, NC * NL), NC, NL)
    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(
            [
                MOI.VectorAffineTerm(
                    i, MOI.ScalarAffineTerm(1.0, assignment_variables[i])
                )
                for i in eachindex(assignment_variables)
            ],
            zeros(NC * NL),
        ),
        MOI.Nonnegatives(NC * NL),
    )
    # "Less than 1.0" constraint is redundant.

    ###
    ### Objective function
    ###

    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            [
                MOI.ScalarAffineTerm(
                    abs(data.customer_locations[i] - j),
                    assignment_variables[i, j],
                )
                for i in 1:NC for j in 1:NL
            ],
            0.0,
        )
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    ###
    ### assignment_variables[i,j] <= facility_variables[j]
    ###     ⟺
    ### assignment_variables[i,j] - facility_variables[j] ∈ R₋
    ###

    terms = Vector{MOI.VectorAffineTerm{Float64}}(undef, 2 * NC * NL)
    row = 1
    for i in 1:NC, j in 1:NL
        terms[2 * row - 1] = MOI.VectorAffineTerm(
            row, MOI.ScalarAffineTerm(1.0, assignment_variables[i, j])
        )
        terms[2 * row] = MOI.VectorAffineTerm(
            row, MOI.ScalarAffineTerm(-1.0, facility_variables[j])
        )
        row += 1
    end

    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(terms, zeros(NC * NL)),
        MOI.Nonpositives(NC * NL),
    )

    ###
    ### sum_j assignment_variables[i,j] == 1
    ###     ⟺
    ### sum_j assignment_variables[i,j] - 1 ∈ 0
    ###

    terms_2 = Vector{MOI.VectorAffineTerm{Float64}}(undef, NC * NL)
    k = 1
    for i in 1:NC, j in 1:NL
        terms_2[k] = MOI.VectorAffineTerm(
            i, MOI.ScalarAffineTerm(1.0, assignment_variables[i, j])
        )
        k += 1
    end
    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(terms_2, -ones(NC)),
        MOI.Zeros(NC),
    )

    ###
    ### sum_j facility_variables[j] - num_facilities ∈ 0
    ###

    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction(
            MOI.VectorAffineTerm.(
                1, MOI.ScalarAffineTerm.(1.0, facility_variables)
            ),
            Float64[-data.num_facilities],
        ),
        MOI.Zeros(1),
    )

    return assignment_variables, facility_variables
end

function generate_glpk_problem(prob, data::PMedianData)
    NL = data.num_locations
    NC = data.num_customers

    facility_variables = [GLPK.glp_add_cols(prob, 1) for i in 1:NL]

    for j in 1:NL
        GLPK.glp_set_col_bnds(prob, facility_variables[j], GLPK.GLP_DB, 0.0, 1.0)
    end

    assignment_variables = [
        GLPK.glp_add_cols(prob, 1)
        for i in 1:data.num_customers, j in 1:data.num_locations
    ]

    for i in 1:data.num_customers, j in 1:data.num_locations
        GLPK.glp_set_col_bnds(prob, assignment_variables[i, j], GLPK.GLP_DB, 0.0, 1.0)
        GLPK.glp_set_obj_coef(prob, assignment_variables[i, j], abs(data.customer_locations[i] - j))
    end
    GLPK.glp_set_obj_dir(prob, GLPK.GLP_MIN)

    I = Cint[0] # Extra padding because GLPK starts reading from the second element
    J = Cint[0]
    V = Float64[0.0]
    for i in 1:data.num_customers, j in 1:data.num_locations
        # assignment_variables[i,j] <= facility_variables[j]
        index = GLPK.glp_add_rows(prob, 1)
        GLPK.glp_set_row_bnds(prob, index, GLPK.GLP_UP, 0.0, 0.0)
        push!(I, index)
        push!(J, assignment_variables[i, j])
        push!(V, 1.0)
        push!(I, index)
        push!(J, facility_variables[j])
        push!(V, -1.0)
    end
    for i in 1:data.num_customers
        # sum_j assignment_variables[i,j] = 1
        index = GLPK.glp_add_rows(prob, 1)
        GLPK.glp_set_row_bnds(prob, index, GLPK.GLP_FX, 1.0, 1.0)
        for j in 1:data.num_locations
            push!(I, index)
            push!(J, assignment_variables[i, j])
            push!(V, 1.0)
        end
    end
    # sum_j facility_variables[j] = num_facilities
    index = GLPK.glp_add_rows(prob, 1)
    GLPK.glp_set_row_bnds(prob, index, GLPK.GLP_FX, data.num_facilities, data.num_facilities)
    for j in 1:data.num_locations
        push!(I, index)
        push!(J, facility_variables[j])
        push!(V, 1.0)
    end
    GLPK.glp_load_matrix(prob, length(I) - 1, I, J, V)
    return
end

function solve_glpk_direct(data::PMedianData; time_limit_sec=Inf)
    @timeit "GLPK direct" begin
        prob = GLPK.glp_create_prob()
        param = GLPK.glp_smcp()
        GLPK.glp_init_smcp(param)
        param.msg_lev = GLPK.GLP_MSG_ERR
        if isfinite(time_limit_sec)
            param.tm_lim = 1000 * time_limit_sec
        end
        @timeit "generate" generate_glpk_problem(prob, data)
        @timeit "solve" GLPK.glp_simplex(prob, param)
        objval = GLPK.glp_get_obj_val(prob)
        GLPK.glp_delete_prob(prob)
    end
    return objval
end

function generate_scs_problem(data::PMedianData)
    # Equalities must be ordered first.
    num_variables = 0
    facility_variables = [num_variables += 1 for i in 1:data.num_locations]
    assignment_variables = [num_variables += 1 for i in 1:data.num_customers, j in 1:data.num_locations]

    b = Float64[]
    I = Int[]
    J = Int[]
    V = Float64[]

    num_constraints = 0
    for i in 1:data.num_customers
        # sum_j assignment_variables[i,j] = 1
        num_constraints += 1
        push!(b, 1.0)
        for j in 1:data.num_locations
            push!(I, num_constraints)
            push!(J, assignment_variables[i, j])
            push!(V, 1.0)
        end
    end
    # sum_j facility_variables[j] = num_facilities
    num_constraints += 1
    push!(b, data.num_facilities)
    for j in 1:data.num_locations
        push!(I, num_constraints)
        push!(J, facility_variables[j])
        push!(V, 1.0)
    end

    # Now inequality constraints in the form b - a'x >= 0.

    for i in 1:data.num_customers, j in 1:data.num_locations
        # assignment_variables[i,j] <= facility_variables[j]
        num_constraints += 1
        push!(b, 0.0)
        push!(I, num_constraints)
        push!(J, assignment_variables[i, j])
        push!(V, 1.0)
        push!(I, num_constraints)
        push!(J, facility_variables[j])
        push!(V, -1.0)
    end

    # All variables are nonnegative.
    for v in 1:num_variables
        num_constraints += 1
        push!(b, 0.0)
        push!(I, num_constraints)
        push!(J, v)
        push!(V, -1.0)
    end

    # facility_variables <= 1
    for j in 1:data.num_locations
        num_constraints += 1
        push!(b, 1.0)
        push!(I, num_constraints)
        push!(J, facility_variables[j])
        push!(V, 1.0)
    end

    c = zeros(num_variables)
    for i in 1:data.num_customers, j in 1:data.num_locations
        c[assignment_variables[i, j]] = abs(data.customer_locations[i] - j)
    end
    num_zero_cones = data.num_customers + 1
    num_positive_orthants = num_constraints - num_zero_cones

    A = sparse(I, J, V, num_constraints, num_variables)

    return num_constraints, num_variables, A, b, c, num_zero_cones, num_positive_orthants,
        Int[], Int[], 0, 0, Float64[]
end

function solve_scs_direct(data::PMedianData; max_iters)
    @timeit "SCS direct" begin
        @timeit "generate" scs_prob = generate_scs_problem(data)
        @timeit "solve" solution = SCS.SCS_solve(
            SCS.IndirectSolver,
            scs_prob...;
            max_iters=max_iters,
            acceleration_lookback=0,
            verbose=0,
        )
    end
    objval = scs_prob[5]' * solution.x
    return objval
end

function solve_moi(data::PMedianData, optimizer; vector_version, params)
    cache = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        optimizer()
    )
    model = MOI.Bridges.full_bridge_optimizer(cache, Float64)
    # resetting optimizer is necessary to use copy_to to transfer all the data
    # from caching optimizer to GLPK in a single batch. If that is not called
    # constraints are passed one by one during model generation phase.
    MOI.Utilities.reset_optimizer(cache)
    for (param, value) in params
        MOI.set(model, param, value)
    end
    @timeit "generate" x, y = if vector_version
        generate_moi_problem_vector(model, data)
    else
        generate_moi_problem(model, data)
    end
    @timeit "solve" MOI.optimize!(model)
    return MOI.get(model, MOI.ObjectiveValue())
end


function solve_glpk_moi(data::PMedianData; vector_version, time_limit_sec=Inf)
    params = []
    if isfinite(time_limit_sec)
        push!(params, (MOI.TimeLimitSec(), time_limit_sec))
    end
    s_type = vector_version ? "vector" : "scalar"
    @timeit(
        "GLPK MOI $(s_type)",
        solve_moi(
            data, GLPK.Optimizer; vector_version=vector_version, params=params
        )
    )
end

function solve_scs_moi(data::PMedianData; vector_version, max_iters::Int)
    params = Tuple{MOI.RawParameter, Int}[
        (MOI.RawParameter("max_iters"), max_iters),
        (MOI.RawParameter("verbose"), 0),
        (MOI.RawParameter("acceleration_lookback"), 0)
    ]
    s_type = vector_version ? "vector" : "scalar"
    @timeit(
        "SCS MOI $(s_type)",
        solve_moi(
            data, SCS.Optimizer; vector_version=vector_version, params=params
        )
    )
end

function run_benchmark(;
    num_facilities, num_customers, num_locations, time_limit_sec, max_iters
)
    Random.seed!(10)
    reset_timer!()
    data = PMedianData(num_facilities, num_customers, num_locations, rand(num_customers) .* num_locations)
    GC.gc()
    glpk_moi_vector_obj = solve_glpk_moi(data, vector_version=true, time_limit_sec=time_limit_sec)
    @show glpk_moi_vector_obj
    GC.gc()
    glpk_moi_scalar_obj = solve_glpk_moi(data, vector_version=false, time_limit_sec=time_limit_sec)
    @show glpk_moi_scalar_obj
    GC.gc()
    glpk_direct_obj = solve_glpk_direct(data, time_limit_sec=time_limit_sec)
    @show glpk_direct_obj
    GC.gc()
    scs_moi_vector_obj = solve_scs_moi(data, vector_version=true, max_iters=max_iters)
    @show scs_moi_vector_obj
    GC.gc()
    scs_moi_scalar_obj = solve_scs_moi(data, vector_version=false, max_iters=max_iters)
    @show scs_moi_scalar_obj
    GC.gc()
    scs_direct_obj = solve_scs_direct(data, max_iters=max_iters)
    @show scs_direct_obj

    print_timer()
    println()
end

function run_paper_benchmark(L)
    to = TimerOutputs.get_defaulttimer()
    s = ""
    if isfile("cvxpy.json")
        cvxpy = JSON.parsefile("cvxpy.json")
    else
        cvxpy = nothing
        println("WARNING: cvxpy.json not found.")
    end
    for l in L
        row_starts = Dict(
            "generate" => "\t\t\\multirow{3}{*}{$(l)} & generate",
            "solve" => "\t\t& load",
            "total" => "\t\t& total",
        )
        row_computation = Dict(
            "generate" => (to, key) -> round(
                TimerOutputs.time(to[key]["generate"]) / 1e9; digits = 2
            ),
            "solve" => (to, key) -> round(
                TimerOutputs.time(to[key]["solve"]) / 1e9; digits = 2
            ),
            "total" => (to, key) -> round(
                (
                    TimerOutputs.time(to[key]["generate"]) +
                    TimerOutputs.time(to[key]["solve"])
                ) / 1e9;
                digits = 2
            ),
        )
        run_benchmark(
            num_facilities = 100,
            num_customers = 100,
            num_locations = l,
            time_limit_sec = 0.001,
            max_iters = 1,
        )
        for row in ["generate", "solve", "total"]
            s *= row_starts[row]
            for key in [
                "GLPK MOI scalar",
                "GLPK MOI vector",
                "GLPK direct",
                "GLPK CVXPY",
                "SCS MOI scalar",
                "SCS MOI vector",
                "SCS direct",
                "SCS CVXPY",
            ]
                t = if key == "GLPK CVXPY"
                    cvxpy !== nothing ? cvxpy["glpk"]["$(l)"][row] : NaN
                elseif key == "SCS CVXPY"
                    cvxpy !== nothing ? cvxpy["scs"]["$(l)"][row] : NaN
                else
                    row_computation[row](to, key)
                end
                s *= Printf.@sprintf " & %.2f" t
            end
            s *= " \\\\\n"
        end
    end
    println(s)
    return s
end

# run_paper_benchmark([1_00, 5_00, 10_00, 50_00])

end # module
