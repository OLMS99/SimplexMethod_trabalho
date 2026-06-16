c = [1; 1]
b = [12; 1; 4]
A = [[2 3];[4 3];[3 2]]
eqs = ["<=", ">=", ">="]

#println(size(c),size())
include("simplexOrigem.jl")
using Main.SimplexMethod
SimplexMethod.simplex_method(c, A, b, eqs, "MAX")
