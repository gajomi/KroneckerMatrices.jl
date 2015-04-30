abstract KroneckerMatrix{T} <: AbstractMatrix{T}

#sizes, indexing, etc
sizes(C::KroneckerMatrix) = [size(M) for M in terms(C)]
size(C::KroneckerMatrix) = map(*,sizes(C)...)

#function hasmixedproductproperty(C::KroneckerMatrix,D::KroneckerMatrix)
#    (MCO,NCO),(MCI,NCI) = sizes(C)
#    (MDO,NDO),(MDI,NDI) = sizes(D)
#    return NCO==MDO && NCI==MDI
#end
#hasmpp = hasmixedproductproperty

function hindexes(i::Int64,widths)
    return [(mod(div(i-1,u),w)+1)::Int64 for (u,w) in zip(cumprod([1, widths...]),widths)]
end

function kronindexes(C::KroneckerMatrix,i::Int64,j::Int64)
    iwidths,jwidths = zip(sizes(C)...)
    hi,hj = hindexes(i,iwidths),hindexes(j,jwidths)
    return (hi,hj)
end

function convert{T}(::Type{Matrix{T}},C::KroneckerMatrix)
    return kronproduct(terms(C)...)
end
