module KroneckerMatrices

import Base:eltype,convert,size,length,
            rank,det,trace,transpose,ctranspose,inv,
            svdvals,eigvals

export KroneckerProduct,⊗

issquare(A::AbstractMatrix) = ==(size(A)...)

include("products.jl")

end
