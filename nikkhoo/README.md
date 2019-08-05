# pCDM: point Compound Dislocation Model


Nikkhoo et al. [2016] model calculates analytical solution for surface displacements due to a combination of 3 mutually orthogonal tensile displocations in an elastic half-space. This model can be used to simulate inflation or deflation of a volumetric source in the far field, like magmatic instrusion under a volcano. The model is able to approximate any shape: isotropic, sill, dyke, pipe, or any ellipsoid in all orientations in space.

## pcdm.m
The proposed Matlab script is a literal transcription of the Nikkhoo's equations, except for the volume potencies that have been rededined as a total dV and two dimensionless shape parameters A and B.

The equations are also vectorized for (x,y) coordinates and source depth.

Type "doc pcdm" for help, syntax and example, and see script comments for details.

## pcdmdesc.m
A small script to transform source shape parameters A and B to a human-readable string. It uses a default 10% tolerancy. Examples:

```
>> pcdmdesc(1,0)
    'sill'
>> pcdmdesc(0,0.5)
    'vertical pipe'
>> pcdmdesc(1/3,0.5)
    'isotropic'
>> pcdmdesc(0,1)
    'vertical NS dyke'
```
Type "doc pcdmdesc" for help, syntax and examples.

## plotpcdm.m
A first tentative script to represent a pCDM source on a graph, in 3-D or 2-D projections.