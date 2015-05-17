abstract AbstractKronMat{T} <: AbstractMatrix{T}
abstract AbstractKronVec{T} <: AbstractVector{T}
typealias AbstractKronVecOrMat Union(AbstractKronMat,AbstractKronVec)

#sizes, indexing, etc
sizes(C::AbstractKronVecOrMat) = [size(M) for M in terms(C)]
size(C::AbstractKronVecOrMat) = map(*,sizes(C)...)

lengths(C::AbstractKronVec) = [length(term) for term in terms(C)]
length(C::AbstractKronVec) = prod(lengths(C))

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

function kronindexes(C::AbstractKronVec,i::Int64)
    widths = lengths(C)
    return hindexes(i,widths)
end

function kronindexes(C::AbstractKronMat,i::Int64,j::Int64)
    iwidths,jwidths = zip(sizes(C)...)
    hi,hj = hindexes(i,iwidths),hindexes(j,jwidths)
    return (hi,hj)
end

function convert{T}(::Type{Matrix{T}},C::AbstractKronMat)
    return kronprod(terms(C)...)
end

function full{T}(C::AbstractKronMat{T})
    M,N = size(C)
    return T[C[i,j] for i=1:M, j = 1:N]
end

function full{T}(C::AbstractKronVec{T})
    N = length(C)
    return T[C[i] for i=1:N]
end
