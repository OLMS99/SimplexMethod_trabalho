module GomoryCuts
    using AssignmentProblems

    using CSV
    using DataFrames
    using LinearAlgebra

    include("../SimplexBase/simplexOrigem.jl")
    using .SimplexMethod

    export solve_gap
    
    function build_lp(data)
        m, n = size(data.costs) ## agents x tasks
        n_vars = m*n

        idx(i, j) = (i - 1) * n + j

        ## Restrições de assignment - 1 agente por tarefa
        A_eq = zeros(n, n_vars) # tasks x n_vars
        b_eq = ones(n) # tasks
        for j in 1:n
            for i in 1:m
                A_eq[j, idx(i,j)] = 1.0
            end
        end
        sign_eq = fill("==", n)

        ## Restrições de capacidade - limite de capacidade por agente
        A_cap = zeros(m, n_vars) # agents x n_vars
        b_cap = copy(data.capacities) # agents
        for i in 1:m
            for j in 1:n
                A_cap[i, idx(i,j)] = data.consumptions[i,j]
            end
        end
        sign_cap = fill("<=", m)

        ## Restrições de limite superior para Assignments
        A_sup = Matrix{Float64}(I, n_vars, n_vars) # n_vars x n_vars
        b_sup = ones(n_vars)
        sign_sup = fill("<=", n_vars)

        A_full = vcat(A_eq, A_cap, A_sup)
        b_full = vcat(b_eq, b_cap, b_sup)
        sign_full = vcat(sign_eq, sign_cap, sign_sup)
        
        c = vec(data.costs)

        println(size(A_full))
        
        println(size(b_full))
        
        println(size(sign_full))
        return c, A_full, b_full, sign_full
    end

    function find_fractions(base, tol = 1e-5)
        return findall(v -> tol < v < 1.0 - tol, x)
    end

    function find_cut_idx(b, frac_idx)
        cut_idx = -1
        cut_val = typemax(Float64)
        for i in frac_idx
            if !isapprox(b[i], round(b[i]), 1e-5) && b[i] < cut_val
                cut_idx = i
                cut_val = x
            end
        end
        return cut_idx
    end

    function cut_row!(A, b, k)
        new_row = floor.(A[k,:]) - A[k,:]
        new_rhs = floor.(b[k,:]) - b[k,:]

        #push!(A, new_row')
        #push!(b, new_rhs)
        
        A = vcat(A, new_row')
        b = vcat(b, new_rhs)
        return A, b
    end

    function gomory_cuts(data, time_limit = 60)
        c, A, b, sign = build_lp(data)

        start_time = time()
        result = 0

        obj_value, tableau, status = SimplexMethod.simplex_method(c, A, b, sign, "MIN")
        while true
            global result = obj_value

            c = tableau.z_c
            A = tableau.Y
            b = tableau.x_B
            base = tableau.base_idx

            #if status != "optimal"
            #    print("Simplex failed to converge")
            #    break
            #end

            if time() - start_time > time_limit
                print("Time limit exceeded")
                break
            end

            frac_idx = find_fractions(base)
            cut_idx = find_cut_idx(b, frac_idx)

            if cut_idx == -1
                print("No valid cut found")
                break
            end
            
            A, b = cut_row(A, b, cut_idx)

            eqs = fill("<=", length(b))
            obj_value, tableau, status = SimplexMethod.simplex_method(c, A, b, eqs, "MIN")
        end

        elapsed_time = time() - start_time
        return result, elapsed_time
    end

end

using AssignmentProblems

data = loadAssignmentProblem(:a05100)

c, A, b, eqs = GomoryCuts.gomory_cuts(data)
