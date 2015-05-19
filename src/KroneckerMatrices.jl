module KroneckerMatrices

using Iterators

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,
            svdfact,svd,svdvals,eigfact,eigvals,chol,qrfact,lufact

export AbstractKronMat,AbstractKronVec,AbstractKronMatOrVec,
       KronProdMat,KronProdVec,terms,sizes,full,
       kronprod,kronpow

include("utils.jl")
include("fullkronfuns.jl")
include("abstract.jl")
include("products.jl")
include("factorization.jl")

end
