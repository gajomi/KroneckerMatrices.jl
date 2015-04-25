#Kronecker Products
type KroneckerProduct{T} <: KroneckerMatrix{T}
    outer::AbstractMatrix{T}
    inner::AbstractMatrix{T}
end

⊗(outer,inner) =  KroneckerProduct(outer,inner)

function terms(C::KroneckerProduct)
    return (C.outer,C.inner)
end

function getindex(C::KroneckerProduct, i::Integer,j::Integer)
    Ms = terms(C)
    hk,hl = kronindexes(C,i,j)
    return prod([M[k,l] for (M,k,l) in zip(Ms,hk,hl)])
end

#unary operations
for f = (:ctranspose,:tranpose,:inv)
    @eval ($f)(C::KroneckerProduct) = KroneckerProduct(($f)(C.outer),($f)(C.inner))
end

#a "checksquare" macro could be helpful here
trace(C::KroneckerProduct) = trace(C.outer)*trace(C.inner)
rank(C::KroneckerProduct) = rank(C.outer)*rank(C.inner)

function det(C::KroneckerProduct)
    (K,L),(M,N) = sizes(C)
    K*M==L*N || throw(DimensionMismatch("matrix is not square"))
    #if either inner or outer not square C is rank deficient
    return (K==L && M==N) ? det(C.outer)^M*det(C.inner)^K : 0.0
end

^{T<:Integer}(C::KroneckerProduct,n::T) = KroneckerProduct(C.outer^n,C.inner^n)

#factorizations (eig,svd,etc)
function eigvals(C::KroneckerProduct)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    issquare(C.outer) && issquare(C.inner) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    return vec([λ*μ for λ=eigvals(C.outer), μ=eigvals(C.inner)])
end

svdvals{T}(C::KroneckerProduct{T}) = vec([σ*τ for σ=svdvals(C.outer), τ=svdvals(C.inner)])
