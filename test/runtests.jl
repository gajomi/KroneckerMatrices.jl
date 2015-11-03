using KroneckerMatrices
using Base.Test

#square matrix nonsquare factors
Q,R = 3,4
S,T = 4,3
Aouter = rand(Q,R)
Ainner = ones(S,T)

#square matrix, square factors
Q,R = 3,3
S,T = 4,4
Bouter = rand(Q,R)
Binner = ones(S,T)

#tests for square Kronecker Products
for (outer,inner) in zip(Any[Aouter Bouter],Any[Ainner Binner])
    M = KronProdMat(outer,inner)
    Mfull = kronprod(outer,inner)
    @test size(M) == size(Mfull)
    K,L = size(M)
    for i = 1:K
        for j = 1:L
            @test M[i,j] == Mfull[i,j]
        end
    end
    @test full(M') == Mfull'
    @test rank(M) == rank(Mfull)
    @test det(M) == det(Mfull)

    fullsvdvals = sort(svdvals(Mfull),by = -)
    kronsvdvalstrunc = sort(full(svdvals(M)),by = -)
    kronsvdvals = zeros(eltype(M),length(fullsvdvals))
    kronsvdvals[1:length(kronsvdvalstrunc)] = kronsvdvalstrunc
    for (σ,τ) in zip(kronsvdvals,fullsvdvals)
        @test_approx_eq_eps(σ,τ,100*eps())
    end
end

#tests for square Kronecker products with square factors
Cimpl = KronProdMat(Bouter,Binner)
Cexpl = kronprod(Bouter,Binner)

Λimpl,Vimpl = eig(Cimpl)
Λimpl,Vimpl = full(Λimpl),full(Vimpl)
Λexpl,Vexpl = eig(Cexpl)
Λimplflat = [abs(λ)>10*eps() ? 1. : 0. for λ in Λimpl]
Λexplflat = [abs(λ)>10*eps() ? 1. : 0. for λ in Λexpl]
#note: this is kind of a wierd test. Needs explantion
err = maximum(abs(Vimpl*diagm(Λimplflat)*Vimpl'-Vexpl*diagm(Λexplflat)*Vexpl'))
@test_approx_eq_eps(0.,err,10*eps())


#for (λ1,λ2) in zip(sort(vec(convert(Matrix{eltype(B)},eigvals(B))),by=real),sort(eigvals(Bfull),by=real))
#    @test_approx_eq_eps(λ1,λ2,10*eps())
#end
