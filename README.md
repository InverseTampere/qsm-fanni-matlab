# QSM-FaNNI

[![DOI](https://zenodo.org/badge/92806841.svg)](https://zenodo.org/badge/latestdoi/92806841)

Quantitative structure models - Foliage and needles naive insertion algorithm MATLAB implementation. Written and tested with Matlab version R2016b.

![Test result](https://github.com/InverseTampere/qsm-fanni-matlab/raw/master/src/test_result.png)

## Description

QSM-FaNNI can be used to generate a leaf cover for a quantitative structure model (QSM) of a tree. The shape of the leaves, as well as, their location, orientation and size can be controlled by the user. All leaves have identical basis geometry - currently the geometry can consist of any number of triangular faces - that is then manipulated by scaling, translation and rotation. 

The distribution of leaf material on the QSM is controlled by the leaf area density distribution (LADD), leaf size by the leaf size distribution (LSD) and leaf orientation by the leaf orientation distribution (LOD). Furthermore, the position, orientation and length of the petioles connecting the leaves to the QSM are controlled by additional distributions. All of the mentioned distributions are Matlab functions with a fixed interface (certain inputs and outputs). Default functions are included by the user is encouraged to create their own.

The program generates candidate leaves that are accepted to the final leaf cover if they do not intersect the with QSM geometry or any accepted leaves. If an intersection occurs the program tries to modify the leaf parameters (location, scale and/or orientation) with any number of user-customizable transformations, before finally discarding the candidate if intersections still persist.

## Basic usage

The program relies on a number of classes that have to be added to the current path before running QSM-FaNNI.

```
addpath('classes/');
```

The main program the function `qsm_fanni` defined in *qsm_fanni.m*. The function receives the QSM as a `QSMB` object (or a subclass object), an initialized `LeafModel` object (containing the leaf basis geometry), and total leaf area to be distributed. The leaf area parameter can have two components; one for the initial leaf area to be generated, and one for the target leaf area. This can be used to increase the probability that the target area is reached, even if some of the generated leaves are discarded due to unavoidable intersections.

The main function call has the following formats:
```
Leaves = qsm_fanni(QSM,Leaves,LeafArea);
[Leaves, NAccepted] = qsm_fanni(QSM,Leaves,LeafArea);
[Leaves, NAccepted] = qsm_fanni(QSM,Leaves,LeafArea,...);
```
where the inputs are self-explanatory and the output variable `Leaves` is a updated version of the leaf model, containing the resulting leaf cover, and `NAccepted` is the number of accepted leaves in the model. 

### Optional arguments

Multiple optional arguments exist and they are given in name-value pairs ('ValueName',Value) as explained below. The names are not case-sensitive. In the descriptions NLeaf is the number of candidate leaves.


**'TransformConfig'**

Matrix that contains the scaling, rotation and translation configurations to apply one-by-one when a leaf intersects with other leaves or blocks. Each row of the matrix is a single configuration, columns are as follow:

1. scale width  (x, side dir.)
2. scale length (y, front dir.)
3. scale height (z, normal dir.)
4. rotate around x (side)
5. rotate around y (front)
6. rotate around z (normal)
7. translate x (side)
8. translate y (front)
9. translate z (normal)

By default the matrix is set to:

```
x    y    z     x     y    z    x    y    z
-------------------------------------------
1    1    1     0  pi/6    0    0    0    0
1    1    1     0 -pi/6    0    0    0    0
1    1    1     0  pi/2    0    0    0    0
1    1    1  pi/6     0    0    0    0    0
1    1    1 -pi/6     0    0    0    0    0
1    1    1 -pi/6  pi/2    0    0    0    0
1    1    1  pi/6  pi/2    0    0    0    0
 ```

**'RandomTransform'**

When set to true transformation configuration are tried in random order rather than row-by-row. Default: `false`

**'AreaFunction'**

A handle to a function that distributes available total volume to blocks. The given function must have the following API.

Inputs:
1. Number of blocks.
2. Parameters of each block (struct).
3. Optional extra parameters.

Outputs:
1. Vector of relative leaf area in each block. Should sum up to one.

**'SizeFunction'**

A handle to a function that determines the number of leaves per block, and the size of each leaf, i.e., transforms area per block to leaf sizes. The function should have the following inputs and outputs:

Inputs:
1. Vector of area for each block.
2. Parameters of each block (struct).
3. Area of leaf base.
4. Dimensions of leaf base.
5. Optional extra parameters.

Outputs:
1. Matrix of leaf dimensions (NLeaf x 3).
2. Vector of leaf parents, i.e., block ids.
3. Maximum absolute dimension of any leaf, used for voxelization.

**'AngleFunction'**

A handle to a function that determines the direction and normal of each leaf. The function should have the following inputs and outputs:

Inputs:
1. Parameters of each block (struct).
2. Parent index of each leaf (block id).
3. Direction of the petiole of each leaf.
4. Optional extra parameters.

Outputs:
1. Leaf directions (NLeaf x 3).
2. Leaf normals (NLeaf x 3).

**'AreaFunctionParameters'**

Cell array of extra parameters to pass to the LADD function.

**'SizeFunctionParameters'**

Cell array of extra parameters to pass to the LSD function.

**'AngleFunctionParameters'**

Cell array of extra parameters to pass to the LOD function.

**'Seed'**

Set the seed for the random number generator before the leaf generation. The values is used in `reset(RandStream.getGlobalStream,<value>)`;

**'Verbose'**

Boolean that controls printing information on the process status.

## Supporting classes

### QSMB and QSMBCylindrical

The `QSMB` class contains the minimum API any type of QSM must provide. The class is abstract, and thus, for computations the user must initialize an `QSMBCylindrical` object for cylinder-based models, or create their own subclass of `QSMB` for any other type of models. The class constructor has only one input `ModelData` that can have two different formats. The struct-based format is compatible with the current version of the [TreeQSM](https://github.com/InverseTampere/TreeQSM) result, and thus, the documentation of that project can provide further details. Currently, the following, self-explanatory structure fields are utilized during the construction process (array dimensions in parenthesis, with `N` as the cylinder count and `M` as the branch count):

```
ModelData.cylinder.start            (N x 3)
ModelData.cylinder.axis             (N x 3)
ModelData.cylinder.length           (N x 1)
ModelData.cylinder.radius           (N x 1)
ModelData.cylinder.parent           (N x 1)
ModelData.cylinder.extension        (N x 1)
ModelData.cylinder.branch           (N x 1)
ModelData.cylinder.BranchOrder      (N x 1)
ModelData.cylinder.PositionInBranch (N x 1)
ModelData.treedata.NumberBranches   (1 x 1)
ModelData.branch.order              (M x 1)
ModelData.branch.parent             (M x 1)
ModelData.branch.volume             (M x 1)
ModelData.branch.length             (M x 1)
ModelData.branch.angle              (M x 1)
ModelData.branch.height             (M x 1)
```

The format of the constructor is the following:
```
QSM = QSMBCylindrical(ModelData);
```

### LeafModel and LeafModelTriangle

The `LeafModel` class contains the minimum API any type of leaf model must provide. The class is abstract, and thus, for computations the user must initialize an `LeafModelTriangle` object for leaves with the basis geometry consisting of triangles, or create their own subclass of `LeafModel` for any other type of models. The class constructor can take the following formats:
```
Leaves = LeafModelTriangle(Vert, Faces);
Leaves = LeafModelTriangle(Vert, Faces, Ngons);
Leaves = LeafModelTriangle(Vert, Faces, Ngons, NInit);
```
where `Vert` is a Nx3 matrix with the vertices as rows, `Faces` is a Mx3 matrix with vertex indices of a single triangle on each row. The optional argument `Ngons` is a cell array, there each element contains a vector of vertex indices forming *n*-sided polygons, used for optimization, e.g., when exporting the leaf model. `NInit` is a scalar number that can be used to initialize the objects internal matrices to optimize leaf insertion. After leaf insertion has been completed, the extra rows can be removed with the `trim_slack()` method:
```
Leaves.trim_slack();
```

Internally the `LeafModelTriangle` object stores a single copy of the basis geometry and a set of transformation parameters (translation, scaling, rotation) per each accepted leaf. Should the exact geometry of individual leaves be required, the class offers two methods for computing geometry:
```
Tris =  Leaves.triangles(origin, dir, normal, scale);
[Vertices, Faces] =  Leaves.triangles(origin, dir, normal, scale);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon, Filter);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon, origin, dir, normal, scale);
```
where `origin`, `dir`, `normal`, `scale` are single transformation parameters, `fNgon` is a flag determining whether to use n-sided polygon rather than triangles, when available. Furthermore, `Filter` is a logical vector with a length equal to the number of accepted leaves, and where the value `true` filters in the corresponding leaf. `Tris` is a N x 9 matrix with the following structure:
```
[v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z]
```
where `v1`, `v2` and `v3` are the three vertices of a triangle.
