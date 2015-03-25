module KroneckerMatrices

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,
            svdvals,eigvals

export KroneckerProduct,⊗,
       eltype,convert,size,getindex,
       transpose,ctranspose,
       inv,
       rank,det,trace,
       svdvals,eigvals


issquare(A::AbstractMatrix) = ==(size(A)...)

include("products.jl")

end
