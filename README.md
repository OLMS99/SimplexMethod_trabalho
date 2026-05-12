Implementação da classe Simplex Methood baseada no livro Julia Programming for Operations Research 2/e por Changhyun Kwon,
com adições pedidas nos trabalhos da matéria dada na matéria Otimização Combinatória na PUC-Rio no primeiro semestre 2026

para reproduzir o exemplo que o capítulo usa, utilize este bloco de código:

```
c = [-3; -2; -1; -5;]

b = [7; 3; 8]

A= [([7 3 4 1],"<="),([2 1 1 5],"<="),([1 4 5 2],"<=")]

include("simplexOrigem.jl")
using Main.SimplexMethod
SimplexMethod.canonize_simplex(c, A, b,"MIN")

```
