module KroneckerMatrices

using Iterators

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,*,
            svdfact,svd,svdvals,eigfact,eigvals,chol,qrfact,lufact

export AbstractKronMat,AbstractKronVec,AbstractKronMatOrVec,
       KronProdMat,KronProdVec,KronSumMat,terms,sizes,full,
       kronprod,kronprod,kronsum

include("utils.jl")
include("fullkronfuns.jl")
include("abstract.jl")
include("products.jl")
include("sums.jl")
include("factorization.jl")

end
