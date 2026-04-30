module SimplexMethod

  using LinearAlgebra, Combinatorics, Printf

  export simplex_method

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

  function add_excess_or_slack_variables(A::Array, b::Vector)
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

  function initialize(c, A, b)
    c = Array{Float64}(c)
    A = Array{Float64}(A)
    b = Array{Float64}(b)

    m, n = size(A)

    # Finding an initial BFS
    b_idx, x_B, B = initial_BFS(A,b)

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

end