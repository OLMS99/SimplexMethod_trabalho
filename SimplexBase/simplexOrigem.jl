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

  function update(c, tableau, artificial)
    artificial = Array{Int}(artificial)
    tableau.Y[:, artificial] .= 0 ## zerando coeficientes de variaveis artificial (mais barato que tirar da matriz)

    if length(intersect(artificial, tableau.b_idx)) > 0
      error("Infeasible")
    end
    
    c = Array{Float64}(c)
    obj = 0
    for idx in tableau.b_idx
      coef = c[idx]
      row = 1
      for (row, x) in enumerate(tableau.Y[:, idx])
        if x == 1.0
          break
        end
      end

      obj += coef*tableau.x_B[row]
      for (j,x) in enumerate(tableau.Y[row, :])
        c[j] -= coef*x
      end
    end

    return c, tableau.Y, tableau.x_B, tableau.b_idx, obj
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

  function initial_Origin(A, b, base_idx)
    m, n = size(A)

    #b_idx = n-m+1:n
    b_idx = base_idx
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

  function pivoting!(t::SimplexTableau, verbose::Bool)
    m, n = size(t.Y)

    entering, exiting = pivot_point(t)
    verbose && println("Pivoting: entering = x_$entering, exiting = x_$(t.b_idx[exiting])")

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

  function initialize(c, A, b, base_idx)
    c = Array{Float64}(c)
    A = Array{Float64}(A)
    b = Array{Float64}(b)

    base_idx = Array{Int}(base_idx)
    m, n = size(A)

    # Finding an initial BFS
    #b_idx, x_B, B = initial_BFS(A,b)
    b_idx, x_B, B = initial_Origin(A,b,base_idx)

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

  function solve_simplex(c, A, b, base_idx, verbose, obj = 0.0)
    tableau = initialize(c, A, b, base_idx)
    tableau.obj += obj
    verbose && print_tableau(tableau)

    while !is_optimal(tableau)
      pivoting!(tableau, verbose)
      verbose && print_tableau(tableau)
    end

    opt_x = zeros(length(c))
    opt_x[tableau.b_idx] = tableau.x_B

    return tableau
  end

  function calc_artificial_goal(A, artificial)
    W = zeros(size(A, 2))
    for idx in artificial
      W[idx] = 1
    end
    return  W
  end

  function add_excess_or_slack_variables!(A::Array, c::Array, flags::Array)

    numRestrictions = size(A, 1)
    numExtraVariables = size(flags, 1)
    newMatrixPart = zeros(numRestrictions,numExtraVariables)

    for (idx, flag) in enumerate(flags)
      push!(c, 0)
      newMatrixPart[flag[2],idx] = clamp(1.0*flag[1],-1,1)
    end

    return c, LinearAlgebra.hcat(A, newMatrixPart)
  end

  function make_variable_list(B, eqs, num_var)
    result = []
    artificial = []

    base_idx = []
    for (i, inequation) in enumerate(eqs)
      if inequation == "<="
        push!(result,(1,i))
        if B[i] < 0
          push!(result, (2,i)) # variavel artificial
          push!(artificial, num_var + length(result))
        end
      end
      if inequation == ">="
        push!(result,(-1,i))
        if B[i] > 0
          push!(result, (2,i)) # variavel artificial
          push!(artificial, num_var + length(result))
        end
      end
      if inequation == "=="
        push!(result,(2,i)) # variavel artificial
        push!(artificial, num_var + length(result))
      end
      push!(base_idx, num_var + length(result))
    end

    return result, artificial, base_idx
  end

  function canonize_simplex(c, A, b, eqs)
    extraVariables, artificial, base_idx = make_variable_list(b, eqs, length(c))

    #constraint_matrix = matrix_construction(A)
    constraint_matrix = A

    c, canonized_A = add_excess_or_slack_variables!(constraint_matrix, c, extraVariables)

    return c, canonized_A, artificial, base_idx
  end

  function simplex_method(c, A, b, eqs, dir = "MAX", verbose = true, canonize = true, artificial = [], base_idx = [])
    direction = 1
    if uppercase(dir) == "MAX"
      global direction = -1
    end

    c = c * direction

    if canonize == true
      verbose && println("Transformando o simplex para versão canônica")
      c, A, artificial, base_idx = canonize_simplex(c, A, b, eqs)
    end

    obj = 0
    if length(artificial) > 0
      W = calc_artificial_goal(A, artificial)
      verbose && println("Simplex de duas fases, minimizando obtejivo artificial W:", W)
      verbose && println("Variaveis artificiais: ", artificial)
      verbose && println("Base Inicial: ", base_idx)
      tableau = solve_simplex(W, A, b, base_idx, verbose)
      verbose && println("Primeira fase resolvida")
      c, A, b, base_idx, obj = update(c, tableau, artificial)
      verbose && println("Iniciando segunda fase do simplex de duas fases")
    end

    tableau = solve_simplex(c, A, b, base_idx, verbose, obj)
    return tableau.obj * direction
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
