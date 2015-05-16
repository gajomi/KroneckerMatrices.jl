#Kronecker Products
type KronProdMat{T} <: KronMat{T}
    terms::Vector{AbstractMatrix{T}}
end

function KronProdMat(terms...)
    T = promote_type([eltype(term) for term in terms]...)
    return KronProdMat(AbstractMatrix{T}[terms...])
end

⊗(terms...) =  KronProdMat(terms...)

function terms(C::KronProdMat)
    return C.terms
end

function getindex(C::KronProdMat, i::Integer,j::Integer)
    Ms = terms(C)
    hk,hl = kronindexes(C,i,j)
    return prod([M[k,l] for (M,k,l) in zip(Ms,hk,hl)])
end

#unary operations
for f = (:ctranspose,:tranpose,:inv)
    @eval ($f)(C::KronProdMat) = KronProdMat(Any[($f)(term) for term in terms(C)]...)
end

for f = (:trace,:rank)
    @eval ($f)(C::KronProdMat) = prod([($f)(term) for term in terms(C)])
end

function det{T}(C::KronProdMat{T})
    Ms,Ns = [[S...] for S in zip(sizes(C)...)]
    bigM,bigN = prod(Ms),prod(Ns)
    bigM == bigN || throw(DimensionMismatch("matrix is not square"))
    if all(Ms .== Ns)
        return prod([det(term)^div(bigM,M) for (term,M) in zip(terms(C),Ms)])
    else
        return zero(T)
    end
end

^{T<:Integer}(C::KronProdMat,n::T) = KronProdMat(Any[term^n for term in terms(C)]...)
