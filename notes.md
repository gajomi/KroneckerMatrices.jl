TODO:
+Fix indexing in kroneckerproducts, adding tests with explicit index calls
+Rewrite eigvals/svdvals to be kronceker product intstances, test explicit convert to matrix types
introduce factorizations separatel from definitions, mirror struction of base
Rename to KronMat
Introduce distibction between vectors and matrices: abstract KronMat and abstract KronVec


Example notebook with: serprinski carpet, hadamard ldu picutres, and maybe a thrid, just for fun
Example notebook: Algebraic graph theory- efficient computation of traces, eigenvalues
Example notebook: Poisson equation with different source terms
Exaple notebook: DNA and/or epigenetic sequence mutation
Example notebook: performance characteristic of slected functions

Note that implementation of eigvals and svdvals may be scerewd up at the moment. In reality these should just be
constructed with kron operation. This is efficeint and will ensure that correct order is kept so as to obtain
the right results along pending implementaion of eigvecs et al.
