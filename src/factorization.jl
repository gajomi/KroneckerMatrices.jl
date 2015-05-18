##########################
# Cholesky Factorization #
##########################
chol(C::KronProdMat) = KronProdMat([chol(term) for term in terms(C)])


###################
#QR Factorization #
###################
function qrfact(C::KronProdMat)
    qs,rs = zip([qr(term) for term in terms(C)])
    return [:Q=>KronProdMat(qs...),:R=>KronProdMat(rs...)]
end


######################
# Eigendecomposition #
######################
function eigfact(C::KronProdMat)
    issquare(C) || throw(DimensionMismatch("matrix not square"))
    all(issquare,terms(C)) || throw(DimensionMismatch("currently only supported for square inner and outer"))
    λs,us = zip([eig(term) for term in terms(C)]...)
    return [:vectors=>KronProdMat(us...),:values=>KronProdVec(λs...)]
end

eigvals(C::KronProdMat) = eigfact(C)[:values]


#######
# SVD #
#######
function svdfact(C::KronProdMat)
    us,σs,vs = zip([svd(term) for term in terms(C)]...)
    return [:U=>KronProdMat(us...),:S=>KronProdMat(σs...),:V=>KronProdMat(vs...)]
end

function svd(C::KronProdMat)
    F = svdfact(C)
    return (F[:U],F[:S],F[:V])
end

svdvals{T}(C::KronProdMat{T}) = KronProdVec([svdvals(term) for term in terms(C)]...)
