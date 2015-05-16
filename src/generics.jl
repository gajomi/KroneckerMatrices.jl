abstract KronMat{T} <: AbstractMatrix{T}

#sizes, indexing, etc
sizes(C::KronMat) = [size(M) for M in terms(C)]
size(C::KronMat) = map(*,sizes(C)...)

#function hasmixedproductproperty(C::KronMat,D::KronMat)
#    (MCO,NCO),(MCI,NCI) = sizes(C)
#    (MDO,NDO),(MDI,NDI) = sizes(D)
#    return NCO==MDO && NCI==MDI
#end
#hasmpp = hasmixedproductproperty

#do I need to merely reverse something?
function hindexes(i::Int64,widths)
    widths = reverse(widths)
    return reverse([(mod(div(i-1,u),w)+1)::Int64 for (u,w) in zip(cumprod([1, widths...]),widths)])
end

function kronindexes(C::KronMat,i::Int64,j::Int64)
    iwidths,jwidths = zip(sizes(C)...)
    hi,hj = hindexes(i,iwidths),hindexes(j,jwidths)
    return (hi,hj)
end

function convert{T}(::Type{Matrix{T}},C::KronMat)
    return kronproduct(terms(C)...)
end
