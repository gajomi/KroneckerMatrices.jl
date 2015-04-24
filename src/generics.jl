abstract KroneckerMatrix{T} <: AbstractMatrix{T}

#sizes, indexing, etc
sizes(C::KroneckerMatrix) = (size(C.outer),size(C.inner))#[size(M) for M in C.factors]
size(C::KroneckerMatrix) = map(*,sizes(C)...)#something totally different

function hasmixedproductproperty(C::KroneckerMatrix,D::KroneckerMatrix)
    (MCO,NCO),(MCI,NCI) = sizes(C)
    (MDO,NDO),(MDI,NDI) = sizes(D)
    return NCO==MDO && NCI==MDI
end
hasmpp = hasmixedproductproperty

function outerindex(C::KroneckerMatrix,i,j)
    M,N = sizes(C)[2]
    return (div(i-1,M)+1,div(j-1,N)+1)
end

function innerindex(C::KroneckerMatrix,s,t)
    M,N = sizes(C)[2]
    return (rem(s-1,M)+1,rem(t-1,N)+1)
end

function convert{T}(::Type{AbstractMatrix{T}},C::KroneckerMatrix)
    M,N = size(C)
    return [C[i,j] for i=1:M,j=1:N]
end
