import MathOptInterface
import SCS
import GLPK
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
    facility_variables = MOI.add_variables(model, data.num_locations)
    for v in facility_variables
        MOI.add_constraint(model, v, MOI.GreaterThan(0.0))
        MOI.add_constraint(model, v, MOI.LessThan(1.0))
    end
    assignment_variables = [MOI.add_variable(model) for i in 1:data.num_customers, j in 1:data.num_locations]
    for v in assignment_variables
        MOI.add_constraint(model, v, MOI.GreaterThan(0.0))
        # "Less than 1.0" constraint is redundant.
    end
    objective = MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(abs(data.customer_locations[i] - j), assignment_variables[i,j])
        for i in 1:data.num_customers for j in 1:data.num_locations],
            0.0)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objective)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    for i in 1:data.num_customers, j in 1:data.num_locations
        # assignment_variables[i,j] <= facility_variables[j]
        MOI.add_constraint(model, MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, assignment_variables[i,j]),
            MOI.ScalarAffineTerm(-1.0, facility_variables[j])], 0.0),
            MOI.LessThan(0.0))
    end
    for i in 1:data.num_customers
        # sum_j assignment_variables[i,j] = 1
        MOI.add_constraint(model, MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, assignment_variables[i,j]) for j in 1:data.num_locations],
            0.0),
            MOI.EqualTo(1.0)
        )
    end
    # sum_j facility_variables[j] = num_facilities
    MOI.add_constraint(model, MOI.ScalarAffineFunction(
        MOI.ScalarAffineTerm.(1.0, facility_variables), 0.0),
        MOI.EqualTo{Float64}(data.num_facilities)
    )
    return assignment_variables, facility_variables
end

function generate_glpk_problem(prob, data::PMedianData)
    facility_variables = [GLPK.glp_add_cols(prob, 1) for i in 1:data.num_locations]
    assignment_variables = [GLPK.glp_add_cols(prob, 1) for i in 1:data.num_customers, j in 1:data.num_locations]
    for j in 1:data.num_locations
        GLPK.glp_set_col_bnds(prob, facility_variables[j], GLPK.GLP_DB, 0.0, 1.0)
    end
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
        @timeit "solve" solution = SCS.SCS_solve(SCS.IndirectSolver, 
                                    scs_prob...; max_iters=max_iters,
                                    acceleration_lookback=0, verbose=0)
    end
    objval = scs_prob[5]'*solution.x
    return objval
end

function solve_moi(data::PMedianData, optimizer; params)
    cached = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        optimizer())
    model = MOI.Bridges.full_bridge_optimizer(cached, Float64)
    MOI.Utilities.reset_optimizer(cached)
    for (param, value) in params
        MOI.set(model, param, value)
    end
    @timeit "generate" x, y = generate_moi_problem(model, data)
    @timeit "solve" MOI.optimize!(model)
    return MOI.get(model, MOI.ObjectiveValue())
end


function solve_glpk_moi(data::PMedianData; time_limit_sec=Inf)
    params = []
    if isfinite(time_limit_sec)
        push!(params, (MOI.TimeLimitSec(), time_limit_sec))
    end
    @timeit "GLPK MOI" solve_moi(data, GLPK.Optimizer, params=params)
end

function solve_scs_moi(data::PMedianData; max_iters)
    params = [(MOI.RawParameter("max_iters"), max_iters),
              (MOI.Silent(), true),
              (MOI.RawParameter("acceleration_lookback"), 0)]
    @timeit "SCS MOI" solve_moi(data, SCS.Optimizer, params=params)
end



function run_benchmark(;num_facilities, num_customers, num_locations,
        time_limit_sec, max_iters)
    Random.seed!(10)
    reset_timer!()
    data = PMedianData(num_facilities, num_customers, num_locations, rand(num_customers) .* num_locations)
    
    glpk_moi_obj = solve_glpk_moi(data, time_limit_sec=time_limit_sec)
    @show glpk_moi_obj

    glpk_direct_obj = solve_glpk_direct(data, time_limit_sec=time_limit_sec)
    @show glpk_direct_obj

    scs_moi_obj = solve_scs_moi(data, max_iters=max_iters)
    @show scs_moi_obj

    scs_direct_obj = solve_scs_direct(data, max_iters=max_iters)
    @show scs_direct_obj
    
    print_timer()
    println()
end

# JIT warm-up
run_benchmark(num_facilities=5, num_customers=20, num_locations=10,
    time_limit_sec=Inf, max_iters=10000)

run_benchmark(num_facilities=5, num_customers=20, num_locations=10,
    time_limit_sec=Inf, max_iters=10000)

run_benchmark(num_facilities=10, num_customers=2000, num_locations=1000,
    time_limit_sec=1, max_iters=10)
