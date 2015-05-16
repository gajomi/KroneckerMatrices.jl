issquare(A::AbstractMatrix) = ==(size(A)...)

mat(a::Vector) = reshape(a,length(a),1)
