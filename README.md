# Voxelsolve 
Single-file library (and algorythm itself) to create iso-surfaces in a form of 3D mesh.

## Function
To make a iso-surface, 3D function is being used:
`f(x,y,z, ...)`
Isosurface is a contour around aream, where `f() >= [scalar value]`.

## Method
The key idea of the method is to divide all the prcessed space into voxels.
For each voxel, try to build a small area of iso-surface by solving certain number of equations.
In general, equation is being solved for each cube edge (12 equations in general). By analising soluitions (3D points), build a small part of the mesh based on them.
There aren't likely to be more than 4 dots for each vertex, but there's possibility to triangulate even more points. 

## About mesh 
Surface is a 3d model (mesh). Mesh consists of a set of vertex that from a set of triangles.
Even though this lib is made for 3D graphics, it does not depend on any specific technology (Like OpenGL); But it works well with GLM, because of the same vector implementation. 

## Triangulation issue
The last step of voxelsolve is one of the most complex steps at all. It could be resource-consuming, so maybe there's a possibility to replace it with some lighter verison. Basically, we just cut a cube and try to build the mesh based on these cuts. Somehow, it got pretty painful to programm.
