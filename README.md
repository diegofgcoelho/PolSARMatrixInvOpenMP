# PolSARMatrixInvOpenMP
Fast matrix inversion in the context of PolSAR image processing and classification in C++ using OpenMP.

The inversion of the matrices of interest in performed in place using minimal number of temporary variables. 
The fast inversion (as oposite to Cholesky based approach) is based on the computation of the matrix of cofactors.
