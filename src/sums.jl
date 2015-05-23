type KronSumMat{T} <: AbstractKronMat{T}
    terms::Vector{AbstractMatrix{T}}
end

function KronSumMat(terms...)
    T = promote_type([eltype(term) for term in terms]...)
    KronSumMat(AbstractMatrix{T}[terms...])
end


#ok, evidently there is going to be too much dry violation here.. should correct before proceeding too far.
