kronprod = kron

kronpow(A::Matrix,n::Integer) = n==1? A : kronprod([A for _ in 1:n]...)

function kronsum(terms...)
    Ms,Ns = zip([size(term) for term in terms]...)
    Ms == Ns || throw(DimensionMismatch("term not square"))
    D = [eye(M) for M in Ms]
    return sum([kronprod(setindex!(copy(D),term,i)...) for (i,term) in enumerate(terms)])
end

#function kronrepeatedsum(A::Matrix,n::Integer)
#    B = A
#    for k=1:n
#        B = kronsum(B,A)
#    end
#    return B
#end
