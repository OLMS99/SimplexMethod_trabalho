c = [1; 1]
b = [12; 1; 44]
A= [([2 3],"<="),([4 3],">="),([3 2],">=")]

include("simplexOrigem.jl")
using Main.SimplexMethod
SimplexMethod.simplex_method(c, A, b,"MAX",false)
