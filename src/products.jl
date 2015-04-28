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
trace(C::KroneckerProduct) = prod([trace(term) for term in terms(C)])
rank(C::KroneckerProduct) = prod([rank(term) for term in terms(C)])

function det(C::KroneckerProduct)
    Ms,Ns = [[S...] for S in zip(sizes(C)...)]
    bigM,bigN = prod(Ms),prod(Ns)
    bigM == bigN || throw(DimensionMismatch("matrix is not square"))
    if all(Ms .== Ns)
        return prod([det(term)^div(bigM,M) for (term,M) in zip(terms(C),Ms)])
    else
        return 0.
    end
end

^{T<:Integer}(C::KroneckerProduct,n::T) = KroneckerProduct(C.outer^n,C.inner^n)

#factorizations (eig,svd,etc)
function eigvals(C::KroneckerProduct)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    all(issquare,terms(C)) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    return [prod([λs...]) for λs in product([eigvals(term) for term in terms(C)]...)]
end

svdvals{T}(C::KroneckerProduct{T}) =[prod([σs...]) for σs in product([svdvals(term) for term in terms(C)]...)]
