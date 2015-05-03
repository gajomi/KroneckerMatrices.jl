kronproduct = kron

kronpower(A::Matrix,n::Integer) = n==1? A :kronproduct([A for _ in 1:n]...)


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
