using AssignmentProblems
using JuMP
using HiGHS
using Knapsacks

function solveBinaryModel(data::Knapsack)
    I = collect(1:ni(data))
    W = data.capacity
    w = data.weights
    p = data.profits

    model = Model(HiGHS.Optimizer)

    @variable(model, x[j in I], Bin)
    @objective(model, Max, sum(p[j] * x[j] for j in I))
    @constraint(model, sum(w[j] * x[j] for j in I) <= W)

    set_silent(model)

    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        result = [ j for j in I if value(x[j]) > 1e-7 ]
        return round(Int64, objective_value(model)), result
    end
end

function solvePricingSubproblem(Q, q, c, pi_dual, delta_dual, strategy)
    n = length(q)
    profits = zeros(Int64, n)
    
    for j in 1:n
        profit_val = pi_dual[j] - c[j]
        if profit_val > 0
            profits[j] = round(Int64, 1e6 * profit_val)
        end
    end

    knap = Knapsack(Q, q, profits)

    if strategy == :BinaryModel
        value, subset = solveBinaryModel(knap)
    else
        value, subset = solveKnapsack(knap, strategy)
    end

    cost = sum(c[j] for j in subset; init = 0.0)
    sum_pi = sum(pi_dual[j] for j in subset; init = 0.0)
    
    red_cost = cost - sum_pi - delta_dual

    return red_cost, subset, cost
end

function updateVectors!(subsets, xs, ys, costs, subset, agent, cost, n_jobs, n_agents)
    push!(subsets, subset)
    push!(costs, cost)

    x = zeros(Int64, n_jobs)
    for elem in subset
        x[elem] = 1
    end
    push!(xs, x)

    y = zeros(Int64, n_agents)
    y[agent] = 1
    push!(ys, y)
end

function generateColumn!(model, lambdas, pi, delta, subset, agent, cost, new_id)
    var = @variable(model, lower_bound = 0.0)
    set_name(var, "λ[$new_id]")

    set_objective_coefficient(model, var, cost)
    
    for elem in subset
        set_normalized_coefficient(pi[elem], var, 1.0)
    end
    set_normalized_coefficient(delta[agent], var, 1.0)

    push!(lambdas, var)
end

function solve(data, time_limit = -1)
    n_agents = na(data)
    n_jobs   = nj(data)
    Agents = 1:n_agents
    Jobs = 1:n_jobs

    subsets = []
    xs = []
    ys = []
    costs = []
    lambdas = []

    agent_cap = copy(data.capacities)
    init_jobs = [Int64[] for _ in Agents]
    for j in Jobs
        assigned = false
        for i in Agents
            if data.consumptions[i, j] <= agent_cap[i]
                push!(init_jobs[i], j)
                agent_cap[i] -= data.consumptions[i, j]
                assigned = true
                break
            end
        end
        if !assigned
            error("Failed to initialize")
            return nothing, nothing, nothing, nothing, nothing, nothing
        end
    end

    for i in Agents
        cost = sum(data.costs[i, j] for j in init_jobs[i]; init = 0.0)
        updateVectors!(subsets, xs, ys, costs, init_jobs[i], i, cost, n_jobs, n_agents)
    end

    omegas = 1:length(subsets)

    model = Model(HiGHS.Optimizer)

    @variable(model, lambda[omegas] >= 0)
    push!(lambdas, lambda)

    @objective(model, Min, sum(costs[s] * lambda[s] for s in omegas))
    
    @constraint(model, pi[j in Jobs], sum(xs[s][j] * lambda[s] for s in omegas) == 1)
    @constraint(model, delta[i in Agents], sum(ys[s][i] * lambda[s] for s in omegas) == 1)

    set_silent(model)
    set_time_limit_sec(model, time_limit)

    start_time = time()

    strategy = :Heuristic
    while true
        optimize!(model)

        # Early stop after optimizing but before changing model
        curr_time = time()
        if ((time_limit > 0) && (curr_time - start_time > time_limit)) break end

        pi_dual = dual.(pi)
        delta_dual = dual.(delta)
        
        best_rc = 0.0
        best_agent = -1
        best_subset = Int64[]
        best_cost = 0.0

        for i in Agents
            Q = data.capacities[i]
            q = convert(Vector{Int64}, data.consumptions[i, :])
            c = data.costs[i, :]

            red_cost, subset, cost = solvePricingSubproblem(Q, q, c, pi_dual, delta_dual[i], strategy)

            if red_cost < best_rc
                best_rc = red_cost
                best_agent = i
                best_subset = subset
                best_cost = cost
            end
        end
        
        if best_rc >= -1e-5
            if strategy == :Heuristic
                strategy = :BinaryModel
            else
                println("Column generation finished")
                break 
            end
        else
            updateVectors!(subsets, xs, ys, costs, best_subset, best_agent, best_cost, n_jobs, n_agents)

            generateColumn!(model, lambdas, pi, delta, best_subset, best_agent, best_cost, length(subsets))

            strategy = :Heuristic
        end

    end

    result = objective_value(model)
    model_bound = objective_bound(model)
    real_bound = data.ub

    model_gap = result - model_bound
    real_gap = result - real_bound
    
    jumps_gap = relative_gap(model)
    return result, model_bound, real_bound, model_gap, real_gap, jumps_gap
end

#data = loadAssignmentProblem(:a05100)
#println(solve(data, 60))

function evaluate_all_assignments()
    result = Dict()
    x = 0.0
    for case in names(AssignmentProblems)
        println(case)
        try
            x += 1
            data = loadAssignmentProblem(case)
            ret = @timed solve(data, 60)
            result[case] = [ret.time, ret.value...]
        catch
            # Invalid problem code
        end
    end
    return result
end

res = evaluate_all_assignments()
println(res)

df = DataFrame( [(Case = k, Time=v[1], Value=v[2], model_bound = v[3], real_bound = v[4], model_gap=v[5], real_gap=v[6], relative_gap =v[7]) for (k,v) in res])
CSV.write("resultados.csv", df)