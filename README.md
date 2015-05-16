# KroneckerMatrices
Implementation of types and relevant methods for implicit Kronecker products and Kronecker sums. These differ from matrices constructed from matrices constructed with from kron in base in that the full matrix is not stored, but still allows for values to be computed on demand through the standard call to getindex. Additional standard linear algebra methods are also available when they can be computed effciently.


[![Build Status](https://travis-ci.org/gajomi/KroneckerMatrices.jl.svg?branch=master)](https://travis-ci.org/gajomi/KroneckerMatrices.jl)
