using AssignmentProblems

using JuMP
using HiGHS
using CSV
using DataFrames

function gap_solver(data) 
    num_agents = 1:na(data)
    num_jobs = 1:nj(data)
    
    model = Model(HiGHS.Optimizer)

    @variable(model, x[num_agents,num_jobs], Bin)

    @objective(model, Min, sum(x[i,j]data.costs[i,j] for i in num_agents, j in num_jobs))

    @constraint(model, [i in num_agents], sum(x[i,j]data.consumptions[i,j] for j in num_jobs) <= data.capacities[i])
    @constraint(model, [j in num_jobs], sum(x[i,j] for i in num_agents) == 1)
    
    set_time_limit_sec(model, 60)
    set_silent(model)
    optimize!(model)

    result = objective_value(model)
    model_bound = objective_bound(model)
    real_bound = data.lb

    model_gap = result - model_bound
    real_gap = result - real_bound
    
    jumps_gap = relative_gap(model)
    return result, model_bound, real_bound, model_gap, real_gap, jumps_gap
end

#data = loadAssignmentProblem(:a05100)
#println(gap_solver(data))

function evaluate_all_assignments()
    result = Dict()
    x = 0.0
    for case in names(AssignmentProblems)
        println(case)
        try
            x += 1
            data = loadAssignmentProblem(case)
            ret = @timed gap_solver(data)
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