# C isoline
One-file library to create iso-surfaces.

## Function
To make a iso-surface, 3D function is being used:
`f(x,y,z, ...)`
Isosurface is a contour around aream, where `f() >= [scalar value]`.

## Method
The key idea of the method is to divide all the space into voxels, and for each voxel, try to build a small area of iso-surface by solving certain number of equations.

## About mesh 
Surface is a 3d model (mesh). Mesh consists of a set of vertex that from a set of triangles.
Even though this lib is made for 3D graphics, it does not depend on any specific technology (Like OpenGL); But it works well with GLM, because of the same vector implementation. 
