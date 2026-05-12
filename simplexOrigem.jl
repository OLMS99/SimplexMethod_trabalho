module SimplexMethod

  using LinearAlgebra, Printf

  export simplex_method, canonize_simplex

  mutable struct SimplexTableau
    z_c     ::Array{Float64} # z_j - c_j
    Y       ::Array{Float64} # inv(B) * A
    x_B     ::Array{Float64} # inv(B) * b
    obj     ::Float64        # c_B * x_B
    b_idx   ::Array{Int64}   # indices for basic variables x_B
  end

  function is_nonnegative(x::Vector)
    return length( x[ x .< 0] ) == 0
  end

  function is_nonpositive(z::Array)
    return length( z[ z .> 0] ) == 0
  end

  function initial_BFS(A, b)
    m, n = size(A)

    comb = collect(combinations(1:n, m))
    for i in length(comb):-1:1
      b_idx = comb[i]
      B = A[:, b_idx]
      x_B = inv(B) * b
      if is_nonnegative(x_B)
        return b_idx, x_B, B
      end
    end

    error("Infeasible")
  end

  function initial_Origin(A, b)
    m, n = size(A)

    b_idx = n-m+1:n
    B = A[:, b_idx]
    x_B = inv(B) * b
    if is_nonnegative(x_B)
      return b_idx, x_B, B
    end

    error("Infeasible")
  end

  function print_tableau(t::SimplexTableau)
    m, n = size(t.Y)

    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7*n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    @printf("%6s|", "")
    for j in 1:length(t.z_c)
      @printf("%6.2f ", t.z_c[j])
    end
    @printf("| %6.2f\n", t.obj)

    println(hline)

    for i in 1:m
      @printf("x[%2d] |", t.b_idx[i])
      for j in 1:n
        @printf("%6.2f ", t.Y[i,j])
      end
      @printf("| %6.2f\n", t.x_B[i])
    end

    println(hline)
  end

  function pivoting!(t::SimplexTableau)
    m, n = size(t.Y)

    entering, exiting = pivot_point(t)
    println("Pivoting: entering = x_$entering, exiting = x_$(t.b_idx[exiting])")

    # Pivoting: exiting-row, entering-column
    # updating exiting-row
    coef = t.Y[exiting, entering]
    t.Y[exiting, :] /= coef
    t.x_B[exiting] /= coef

    # updating other rows of Y
    for i in setdiff(1:m, exiting)
      coef = t.Y[i, entering]
      t.Y[i, :] -= coef * t.Y[exiting, :]
      t.x_B[i] -= coef * t.x_B[exiting]
    end

    # updating the row for the reduced costs
    coef = t.z_c[entering]
    t.z_c -= coef * t.Y[exiting, :]'
    t.obj -= coef * t.x_B[exiting]

    # Updating b_idx
    t.b_idx[ findfirst(t.b_idx .== t.b_idx[exiting]) ] = entering
  end

  function pivot_point(t::SimplexTableau)
    # Finding the entering variable index
    entering = findfirst( t.z_c .> 0)[2]
    if entering == 0
      error("Optimal")
    end

    # min ratio test / finding the exiting variable index
    pos_idx = findall( t.Y[:, entering] .> 0 )
    if length(pos_idx) == 0
      error("Unbounded")
    end
    exiting = pos_idx[ argmin( t.x_B[pos_idx] ./ t.Y[pos_idx, entering] ) ]

    return entering, exiting
  end

  function initialize(c, A, b)
    c = Array{Float64}(c)
    A = Array{Float64}(A)
    b = Array{Float64}(b)

    m, n = size(A)

    # Finding an initial BFS
    #b_idx, x_B, B = initial_BFS(A,b)
    b_idx, x_B, B = initial_Origin(A,b)

    Y = inv(B) * A
    c_B = c[b_idx]
    obj = dot(c_B, x_B)

    # z_c is a row vector
    z_c = zeros(1,n)
    n_idx = setdiff(1:n, b_idx)
    z_c[n_idx] = c_B' * inv(B) * A[:,n_idx] - c[n_idx]'

    return SimplexTableau(z_c, Y, x_B, obj, b_idx)
  end

  function is_optimal(t::SimplexTableau)
    return is_nonpositive(t.z_c)
  end

  function simplex_method(c, A, b)
    tableau = initialize(c, A, b)
    print_tableau(tableau)

    while !is_optimal(tableau)
      pivoting!(tableau)
      print_tableau(tableau)
    end

    opt_x = zeros(length(c))
    opt_x[tableau.b_idx] = tableau.x_B

    return opt_x, tableau.obj
  end

  function calc_artificial_goal(A, flags)
    W = zeros(size(A, 2))
    num_naturais = size(A, 2) - size(flags, 1)

    println("W cost initialize as:",W)
    println("AA ", A)
    base_idx = []
    for (k,flag) in enumerate(flags)
      if flag[1] == 2
        W[num_naturais + k] -= 1
        push!(base_idx, num_naturais + k)

        j = flag[2]
        for i in 1:size(A, 2)
          W[i] += A[j,i]
        end
      end
    end
    println("W cost is:",W, flags)
    return  W, base_idx
  end

  function add_excess_or_slack_variables!(A::Array, c::Array, flags::Array)

    println(A)
    numRestrictions = size(A, 1)
    numExtraVariables = size(flags, 1)
    newMatrixPart = zeros(numRestrictions,numExtraVariables)
    println(numRestrictions, ":::", numExtraVariables)


    for (idx, flag) in enumerate(flags)
      newMatrixPart[flag[2],idx] = clamp(1.0*flag[1],-1,1)
    end

    return c, LinearAlgebra.hcat(A, newMatrixPart)
  end

  function make_variable_list(B, constriction_list)
    result = []
    artificial = 0
    for (i, inequation) in enumerate(constriction_list)
      if inequation[2] == "<="
        push!(result,(1,i))
        if B[i] < 0
          push!(result, (2,i)) # variavel artificial
          artificial = 1
        end
      end
      if inequation[2] == ">="
        push!(result,(-1,i))
        if B[i] > 0
          push!(result, (2,i)) # variavel artificial
          artificial = 1
        end
      end
      if inequation[2] == "=="
        push!(result,(2,i)) # variavel artificial
        artificial = 1
      end
    end
    return result, artificial
  end

  function canonize_simplex(c, A, b, dir="MIN")
    if uppercase(dir) == "MAX"
      c = c * -1
    end

    extraVariables, artificial = make_variable_list(b, A)
    println("List complete", extraVariables)

    constraint_matrix = matrix_construction(A)
    println("Constraint matrix built", constraint_matrix)

    c, canonized_A = add_excess_or_slack_variables!(constraint_matrix, c, extraVariables)
    println("Constraint matrix cannonized", canonized_A)

    if artificial == 1
      print("Creating artificial objective")
      W, idx = calc_artificial_goal(canonized_A, extraVariables)
      simplex_method(c, canonized_A, b)
    end

    return canonized_A
    simplex_method(c, canonized_A, b)
  end

  function matrix_construction(constraints)
    m = length(constraints)
    n = length(constraints[1][1])

    constraint_matrix = zeros(m,n)
    for (idx,constraint) in enumerate(constraints)
      for (idy,coef) in enumerate(constraint[1])
        constraint_matrix[idx,idy] = coef
      end
    end
    return constraint_matrix
  end
end
