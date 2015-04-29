#Kronecker Products
type KroneckerProduct{T} <: KroneckerMatrix{T}
    terms::Vector{AbstractMatrix{T}}
end

function KroneckerProduct(terms...)
    return KroneckerProduct(AbstractMatrix{Float64}[terms...])
end

⊗(outer,inner) =  KroneckerProduct(outer,inner)

function terms(C::KroneckerProduct)
    return C.terms
end

function getindex(C::KroneckerProduct, i::Integer,j::Integer)
    Ms = terms(C)
    hk,hl = kronindexes(C,i,j)
    return prod([M[k,l] for (M,k,l) in zip(Ms,hk,hl)])
end

#unary operations
for f = (:ctranspose,:tranpose,:inv)
    @eval ($f)(C::KroneckerProduct) = KroneckerProduct(Any[($f)(term) for term in terms(C)]...)
end

for f = (:trace,:rank)
    @eval ($f)(C::KroneckerProduct) = prod([($f)(term) for term in terms(C)])
end

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

^{T<:Integer}(C::KroneckerProduct,n::T) = KroneckerProduct(Any[term^n for term in terms(C)]...)

#factorizations (eig,svd,etc)
function eigvals(C::KroneckerProduct)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    all(issquare,terms(C)) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    return [prod([λs...]) for λs in product([eigvals(term) for term in terms(C)]...)]
end

svdvals{T}(C::KroneckerProduct{T}) =[prod([σs...]) for σs in product([svdvals(term) for term in terms(C)]...)]
