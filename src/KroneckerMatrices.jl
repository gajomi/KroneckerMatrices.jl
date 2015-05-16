module KroneckerMatrices

using Iterators

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,
            svdvals,eigfact,eigvals

export KronProdMat,terms,sizes,
       kronproduct,kronpower

include("utils.jl")
include("fullkronfuns.jl")
include("generics.jl")
include("products.jl")
include("factorization.jl")

end
