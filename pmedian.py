import cvxpy as cp
import json
import numpy as np
import scipy.sparse
import time

def generate_problem(
    num_facilities, num_customers, num_locations, customer_locations
):
    facility_variables = cp.Variable(num_locations)
    assignment_variables = cp.Variable((num_customers, num_locations))
    facility_locations = np.arange(1, num_locations + 1)
    objective = cp.sum(
        cp.multiply(
            np.abs(customer_locations[:, np.newaxis] - facility_locations),
            assignment_variables
        )
    )
    constraints = [
        # sum_j assignment_variables[i, j] = 1 for all i
        cp.sum(assignment_variables, axis = 1) == 1,
        # sum_j facility_variables[j] = num_facilities
        cp.sum(facility_variables) == num_facilities,
        facility_variables >= 0,
        facility_variables <= 1,
        assignment_variables >= 0
    ]
    # assignment_variables[i,j] <= facility_variables[j] for all i, j
    #
    # Version 1: This verison is very slow.
    #
    # assignment_on_off = [
    #     assignment_variables[i, j] <= facility_variables[j]
    #     for i in range(num_customers) for j in range(num_locations)
    # ]
    #
    # Version 2: A bit less slow.
    #
    # assignment_on_off = list()
    # for j in range(num_locations):
    #     assignment_on_off.append(
    #         assignment_variables[:, j] <= facility_variables[j]
    #     )
    #
    # Version 3: Faster, with a steep readability loss.
    facility_mat = scipy.sparse.block_diag(
        num_locations * [np.ones((num_customers, 1))]
    )
    assignment_on_off = [
        cp.vec(assignment_variables) <= facility_mat @ facility_variables
    ]
    constraints.extend(assignment_on_off)
    prob = cp.Problem(cp.Minimize(objective), constraints)
    return prob, facility_variables, assignment_variables

def solve(
    solver,
    num_facilities,
    num_customers,
    num_locations,
    customer_locations,
    scs_max_iters = 10_000,
    glpk_tm_lim = 5_000,
):
    start_generate = time.time()
    prob, f, x = generate_problem(
        num_facilities, num_customers, num_locations, customer_locations
    )
    data, chain, inverse_data = prob.get_problem_data(solver)

    start_solve = time.time()
    try:
        if solver == cp.GLPK:
            soln = chain.solve_via_data(
                prob,
                data,
                solver_opts={
                    "glpk": {
                        "tm_lim": glpk_tm_lim,
                        "msg_lev": "GLP_MSG_OFF"
                    }
                }
            )
        elif solver == cp.SCS:
            # Specifying acceleration_lookback is needed to avoid cvxpy trying
            # to call SCS twice.
            soln = chain.solve_via_data(
                prob,
                data,
                solver_opts={
                    "max_iters": scs_max_iters,
                    "use_indirect": True,
                    "acceleration_lookback": 0
                }
            )
        else:
            raise Exception("Unrecognized solver")
        prob.unpack_results(soln, chain, inverse_data)
    except cp.SolverError:
        print("Solve didn't finish (possibly hit the limit)")

    solve_time = time.time() - start_solve
    gen_time = start_solve - start_generate
    print("Generation time (sec):", gen_time)
    print("Solve time (sec):", solve_time)
    print("Total time (sec):", gen_time + solve_time)
    print("Obj:", prob.value)
    return {
        'generate': round(gen_time, 2),
        'solve': round(solve_time, 2),
        'total': round(gen_time + solve_time, 2)
    }

# Customer locations are fixed for the small problem so that the objectives can
# be compared with the Julia code.
test_customer_locations = np.array([
    1.1258244478647295,
    3.6831406658084287,
    3.444540231363058,
    0.5664544616214151,
    1.2078054506961555,
    1.7957407667101322,
    3.8181274408522614,
    8.151038332483566,
    2.422083248151139,
    8.1977797040088,
    6.6993139516121625,
    4.530583319523316,
    8.440074795625907,
    6.791898821149229,
    7.252498905219804,
    9.236760031564318,
    0.6609803322813423,
    9.991722180201814,
    1.7165450700728413,
    4.204179066639919,
])

print("SCS small problem")
solve(cp.SCS, num_facilities=5, num_customers=20, num_locations=10,
      customer_locations=test_customer_locations)
print("GLPK small problem")
solve(cp.GLPK, num_facilities=5, num_customers=20, num_locations=10,
      customer_locations=test_customer_locations)

customer_locations = np.random.rand(100) * 1000
timing = {'glpk': {}, 'scs': {}}
for l in [1_000, 5_000, 10_000, 50_000]:
    print("GLPK %s" % l)
    timing['glpk'][l] = solve(
        cp.GLPK,
        num_facilities = 100,
        num_customers = 100,
        num_locations = l,
        customer_locations = customer_locations,
        glpk_tm_lim = 1,
    )
    print()
    print("SCS %s" % l)
    timing['scs'][l] = solve(
        cp.SCS,
        num_facilities = 100,
        num_customers = 100,
        num_locations = l,
        customer_locations = customer_locations,
        scs_max_iters = 1
    )
    print()

with open('cvxpy.json', 'w') as io:
    json.dump(timing, io)
