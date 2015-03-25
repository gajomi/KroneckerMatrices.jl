#Kronecker Products
type KroneckerProduct{T} <: AbstractMatrix{T}
    outer::AbstractMatrix{T}
    inner::AbstractMatrix{T}
end

⊗(outer,inner) =  KroneckerProduct(outer,inner)

#sizes, indexing, etc
sizes(C::KroneckerProduct) = (size(C.outer),size(C.inner))
size(C::KroneckerProduct) = map(*,sizes(C)...)

function hasmixedproductproperty(C::KroneckerProduct,D::KroneckerProduct)
    (MCO,NCO),(MCI,NCI) = sizes(C)
    (MDO,NDO),(MDI,NDI) = sizes(D)
    return NCO==MDO && NCI==MDI
end
hasmpp = hasmixedproductproperty

function outerindex(C::KroneckerProduct,i,j)
    M,N = sizes(C)[2]
    return (div(i-1,M)+1,div(j-1,N)+1)
end

function innerindex(C::KroneckerProduct,s,t)
    M,N = sizes(C)[2]
    return (rem(s-1,M)+1,rem(t-1,N)+1)
end

function getindex(C::KroneckerProduct, i::Integer,j::Integer)
    os,ot = outerindex(C,i,j)
    is,it = innerindex(C,i,j)
    return C.outer[os,ot]*C.inner[is,it]
end

function convert{T}(::Type{AbstractMatrix{T}},C::KroneckerProduct)
    M,N = size(C)
    return [C[i,j] for i=1:M,j=1:N]
end

#equality comparisona and friends (this could get tricky)

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

#decompositions (eig,svd, etc)
function eigvals(C::KroneckerProduct)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    issquare(C.outer) && issquare(C.inner) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    return vec([λ*μ for λ=eigvals(C.outer), μ=eigvals(C.inner)])
end

function svdvals{T}(C::KroneckerProduct{T})
    M = minimum(size(C))#probably don;t want to do this nonsense actually
    vals = zeros(T,M)
    nonzerovals = vec([σ*τ for σ=svdvals(C.outer), τ=svdvals(C.inner)])
    vals[1:length(nonzerovals)] = nonzerovals
    return vals
end
