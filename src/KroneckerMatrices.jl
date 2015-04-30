module KroneckerMatrices

using Iterators

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,
            svdvals,eigvals

export KroneckerProduct,⊗,terms,sizes,
       kronproduct

issquare(A::AbstractMatrix) = ==(size(A)...)

include("fullkronfuns.jl")
include("generics.jl")
include("products.jl")

end
