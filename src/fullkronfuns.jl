function kronproduct(As...)
    n = length(As)
    B = As[1]
    for k=2:n
        B = kron(B,As[k])
    end
    return B
end

function kronpower(A::Matrix,n::Integer)
    B = A
    for k=1:n
        B = kron(B,A)
    end
    return B
end

#function kronsum(A::Matrix,B::Matrix)
#    K,L = size(A)
#    M,N = size(B)
#    K==L || error("arguments must be square matrices")
#    return kronproduct(A,diagm(ones(M)))+kronproduct(diagm(ones(K)),B)
#end

#function kronrepeatedsum(A::Matrix,n::Integer)
#    B = A
#    for k=1:n
#        B = kronsum(B,A)
#    end
#    return B
#end
