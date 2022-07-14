Code to reproduce the figures in "Fast global spectral methods for three-dimensional partial differential equations" by C. Stroessner and D. Kressner.
(IMA J. Numer. Anal. 2022. DOI: https://doi.org/10.1093/imanum/drac030.
Arxiv: https://arxiv.org/abs/2111.04585.)

To perform the experiments you need to additionally download
- Tensor Recursive (https://www.epfl.ch/labs/anchp/wp-content/uploads/2019/05/tensor_recursive.tar.gz) for the blocked recursive solver by Chen and Kressner.
- Fast Poisson Solver (https://github.com/danfortunato/fast-poisson-solvers) for comparing to NADIM by Fortunato and Townsend.
- Chebfun (https://www.chebfun.org/) required for comparing to NADIM.
- Tensorlab (https://www.tensorlab.net/) required for computing CP decompositions efficiently.

The folders spectral_method_3D, tensor_recursive, fast-poisson-solvers and chebfun need to be added to the Matlab path. The fast-poisson-solver additionally requires running the contained setup script by hand.

The files Plot..., Experiment... and RuntimeComparison.m reproduce the plots and data in the paper.

Note: The function ExperimentHelmholtzNonConst requires more than 8GB of RAM.
