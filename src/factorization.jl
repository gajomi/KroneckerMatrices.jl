function eigfact(C::KroneckerProduct,permute::Bool=true, scale::Bool=true)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    all(issquare,terms(C)) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    U = KroneckerProduct([eigvecs(term) for term in terms(C)]...)
    Λ = KroneckerProduct([mat(eigvals(term)) for term in terms(C)]...)
    return [:vectors=>U,:values=>Λ]
end

eigvals(C::KroneckerProduct) = eigfact(C)[:values]

svdvals{T}(C::KroneckerProduct{T}) = KroneckerProduct([mat(svdvals(term)) for term in terms(C)]...)
