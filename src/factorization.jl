function eigfact(C::KronProdMat,permute::Bool=true, scale::Bool=true)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    all(issquare,terms(C)) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    U = KronProdMat([eigvecs(term) for term in terms(C)]...)
    Λ = KronProdMat([mat(eigvals(term)) for term in terms(C)]...)
    return [:vectors=>U,:values=>Λ]
end

eigvals(C::KronProdMat) = eigfact(C)[:values]

svdvals{T}(C::KronProdMat{T}) = KronProdMat([mat(svdvals(term)) for term in terms(C)]...)
