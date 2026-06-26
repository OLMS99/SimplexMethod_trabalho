using AssignmentProblems

using JuMP
using HiGHS
using CSV
using DataFrames

function find_cover(x_val, a, b, tol = 1e-4)
    n, m = size(a) # agents x jobs

    violated_covers = Tuple{Int, Vector{Int}}[]

    for i in 1:n
        submodel = Model(HiGHS.Optimizer)

        @variable(submodel, z[1:m], Bin)

        # Minimize sum_j (1 - x*) * z_j
        @objective(submodel, Min, sum((1.0 - x_val[i, j]) * z[j] for j in 1:m))
        @constraint(submodel, sum(a[i, j] * z[j] for j in 1:m) >= b[i] + tol)

        set_silent(submodel)
        optimize!(submodel)

        if termination_status(submodel) == MOI.OPTIMAL
            ## verify candidate cover is violated
            if objective_value(submodel) < (1.0 - tol) 
                z_sol = value.(z)
                C = findall(val -> val > 0.5, z_sol)
                push!(violated_covers, (i, C))
            end
        end
    end

    #println("Found covers:", violated_covers)
    return violated_covers
end

function gap_solver_cover(data, time_limit = -1) 
    num_agents = 1:na(data)
    num_jobs = 1:nj(data)
    
    model = Model(HiGHS.Optimizer)

    @variable(model, 0 <= x[num_agents,num_jobs] <= 1)


    @objective(model, Min, sum(x[i,j]data.costs[i,j] for i in num_agents, j in num_jobs))

    @constraint(model, [i in num_agents], sum(x[i,j]data.consumptions[i,j] for j in num_jobs) <= data.capacities[i])
    @constraint(model, [j in num_jobs], sum(x[i,j] for i in num_agents) == 1)
    
    set_silent(model)
    set_time_limit_sec(model, time_limit)

    start_time = time()
    
    optimize!(model)

    while true
        xval = value.(x)

        remains = max(5, time() - start_time)
        covers = find_cover(xval, data.consumptions, data.capacities, remains)

        if length(covers) < 1 
            println("No new cover cut found!")
            break 
        end

        for (i, C) in covers
            #println("Adding cover cut for Agent $i (size $(length(C)))")
            @constraint(model, sum(x[i, j] for j in C) <= length(C) - 1)
        end

        optimize!(model)

        curr_time = time()
        if ((time_limit > 0) && (curr_time - start_time > time_limit)) break end
    end

    for i in num_agents, j in num_jobs
        set_binary(x[i, j])
    end
    #unset_silent(model)
    remains = max(5, time() - start_time)
    set_time_limit_sec(model, remains)
    optimize!(model)

    result = objective_value(model)
    model_bound = objective_bound(model)
    real_bound = data.lb

    model_gap = result - model_bound
    real_gap = result - real_bound
    
    jumps_gap = relative_gap(model)
    return result, model_bound, real_bound, model_gap, real_gap, jumps_gap
end

#data = loadAssignmentProblem(:c0520_5)
#println(@timed gap_solver_cover(data, 5))

function evaluate_all_assignments()
    result = Dict()
    x = 0.0
    for case in names(AssignmentProblems)
        println(case)
        try
            x += 1
            data = loadAssignmentProblem(case)
            ret = @timed gap_solver_cover(data, 60)
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