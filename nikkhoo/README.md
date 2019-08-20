# pCDM: point Compound Dislocation Model


Nikkhoo et al. [2017] model calculates analytical solution for surface displacements due to a combination of 3 mutually orthogonal tensile displocations (one horizontal and two vertical) freely oriented in space (3 angles of rotation) in an elastic half-space.

This model can be used to simulate inflation or deflation of a volumetric source observed in the far field, like magmatic instrusion under a volcano and GPS surface stations measurements. The model is able to approximate any shape: isotropic, sill, dyke, pipe, or any ellipsoid in any orientations in space.

![](pcdm_ab.png)

## pcdmv.m
The proposed Matlab script is a literal transcription of the Nikkhoo's equations from original script pCDM.m (see [www.volcanodeformation.com](http://www.volcanodeformation.com)), except for:

- the source coordinates is set to (0,0) so observation points (x,y) are relative to it;
- the volume potencies have been redefined as a total volume variation *dVtot* and two dimensionless shape parameters *A* and *B*, between 0 and 1, defined as follows (see `pcdmdesc.m` below for examples):
	- *A = dVz / dVtot* is horizontal over total volume variation ratio,
	- *B = dVy / (dVx + dVy)* is vertical volume variation ratio;
- the equations have been vectorized for all input parameters excepted Poisson's ratio *nu* which must be a scalar. So massive computation can be achieved on large vectors, matrix or even N-D matrix of inputs.

Type "doc pcdm" for help, syntax and example, and see script comments for details.

## pcdmdesc.m
A small script to transform source shape parameters *A* and *B* to a human-readable string. It uses a default 10% tolerancy. Examples:

	>> pcdmdesc(1,0)
	    'sill'
	>> pcdmdesc(0,0.5)
	    'vertical pipe'
	>> pcdmdesc(1/3,0.5)
	    'isotropic'
	>> pcdmdesc(0,1)
	    'vertical NS dyke'

Type "doc pcdmdesc" for help, syntax and examples.

## plotpcdm.m
A first tentative script to represent a pCDM source on a graph, in 3-D or 2-D projections.

##pcdmv.c
This is a transcription of `pcdm.m` in C language that includes complementary subfunction `mexFunction()` to be compiled as a MEX file (Matlab/Octave executable). To make the binary for your computer architecture, you must install a compiler first then type at the Matlab/Octave command line:

	>> mex pcdmv.c

## Note on vectorization
The three codes have almost the same input parameters but different behaviors with input vectors and matrix:

### original pCDM.m by Nikkhoo (2016)
It accepts vectors or matrix only for X and Y observation points coordinates. Other input parameters X0, Y0, DEPTH, OMEGAX, OMEGAY, OMEGAZ, DVX, DVY, DVZ, and NU must be scalars.

### pcdmv.m
Accepts vectors or matrix for input parameters X, Y, DEPTH, OMEGAX, OMEGAY, OMEGAZ, DV, A, and B and any of them can be also a scalar. Optional NU input parameter must be a scalar. Since output arguments will be set to the same size as input X, if X is a scalar and other input arguments are vectors or matrix, use `repmat(X,...)` to make X also a vector or matrix. For other input arguments, any mixing between matrix and scalars are acceptable (if all matrix have the same size, of course).

### pcdmv.c
Accepts all input arguments as scalar, vector or matrix but all of the same number of arguments or size. Mixing with scalar is prohibited. Output arguments will be set to the size of input arguments.

### performance
This a basic comparison of computational time for one million different random models using a 2.7GHz Intel Core i7 computer.

|code|time|
|:----|--------|
|pCDM.m (original)|4.0 s|
|pcdm.m|1.0 s|
|pcdmv.c (compiled)|0.4 s|
