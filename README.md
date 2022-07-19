# fcc-bz-mesh


## About

This repo allows to build a mesh of face cubic centered BZ.


## Dependency
Dependencies are the following:
```python
matplotlib==3.5.2
numpy==1.22.4
gmsh=4.10.3=hc719622_0
# python-gmsh=4.10.3 (for conda environment only)
```


## What is the BZ

With the file __BZ-pyplot-allBZ-point.py__ you can visualize the convex full of the full BZ we wish to mesh.

Typing
```python
python3 BZ-pyplot-allBZ-point.py
```
produces
![what](rsc/convexhull.png)


## Creating meshes of the BZ

The file __BZ.py__ allows to build a mesh of face cubic centered BZ. Though mesh_type you can
select either 1, 1/2, 1/8 or 1/48 of the BZ. mesh_name controls the name of
the written mesh. You can visualize the geometry and the mesh though -vm
and -gm options. Then you can add refienements box with the --refine_L and
--refine_delta options.

To see all options, just type,
```python
python3 BZ.py --help
```
and read.


For instance,
```python
python3 BZ.py -t 1 -n test-mesh-bz -lmin 3e-2 -lmax 10e-2 -vm --refine_delta 0.1 0.125 0.1 --refine_L 0.1 0.125 0.3
```
produce the following mesh:
![what](rsc/refineLdelta_ex.png)


## What's next ?

You can now use 


## A parser
Just a python parser example in order to export vertices and tetrahedra informations from the .mesh created.
