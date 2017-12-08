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
where `origin`, `dir`, `normal`, `scale` are single transformation parameters, `fNgon` is a flag determining whether to use *n*-sided polygon rather than triangles, when available. Furthermore, `Filter` is a logical vector with a length equal to the number of accepted leaves, and where the value `true` filters in the corresponding leaf. `Tris` is a N x 9 matrix with the following structure:
```
[v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z]
```
where `v1`, `v2` and `v3` are the three vertices of a triangle.

## The test file

The *src* folder contains a simple test file, *test.m* for trying the leaf insertion procedure with example data. A detailed description of the command included in the compact demonstration is given below, with the relevant commands listed before the respective description.

```
addpath('classes/');
```

The command includes the required class definitions into Matlab's path. Those classes include the four classes described above, as well as, a helper class `CubeVoxelization` that is used internally by the main program to describe objects, such as triangles or cylinders, in three-dimensional voxelized spaces.

```
CylData = [
% Rad   Len                  Sta                  Axe CPar CExt BI BO IIB Add
0.300 1.000  0.000  0.000  0.000  0.000  0.000  1.000    0    1  1  0   1   0;
0.200 1.414  0.000  0.000  1.000  0.707  0.000  0.707    1    2  1  0   2   0;
0.100 1.000  1.000  0.000  2.000  1.000  0.000  0.000    2    0  1  1   3   0;
0.200 1.414  0.000  0.000  1.000  0.000  0.707  0.707    1    4  2  1   1   0;
0.100 1.000  0.000  1.000  2.000  0.000  1.000  0.000    4    0  2  1   2   0;
0.200 1.414  0.000  0.000  1.000  0.000 -0.707  0.707    1    6  3  1   1   0;
0.100 1.000  0.000 -1.000  2.000  0.000 -1.000  0.000    6    0  3  1   2   0;
0.200 1.414  0.000  0.000  1.000 -0.707  0.000  0.707    1    8  4  1   1   0;
0.100 1.000 -1.000  0.000  2.000 -1.000  0.000  0.000    8    0  4  1   2   0;
];

BranchData = [
% BOrd   BPar   BVol   BLen   BAng   BHei
     0      0 0.4918 3.4140 0.0000 1.3333;
     1      1 0.2091 2.4140 0.7854 1.7499;
     1      1 0.2091 2.4140 0.7854 1.7499;
     1      1 0.2091 2.4140 0.7854 1.7499;
];

TreeData = [
 1.1192; % Total volume of the tree
 0.4918; % Volume of the trunk
 0.6273; % Total volume of all the branches
 2.0000; % Total height of the tree
 3.4140; % Length of the trunk
10.6560; % Total length of all the branches
 4.0000; % Total number of branches
 1.0000; % Maximum branch order
11.5058; % Total area of cylinders
 0.0000; % DBH = Diameter at breast height, from the QSM
 0.0000; % DBH from cylinder fitted to right place
 0.0000; % DBH from triangulation
];

ModelData = {CylData, BranchData, TreeData};
```

Properties of the cylinders, branches and the complete *tree*, to be contained in the example QSM. The data are defined in three separate matrices in the test file. However, the user may also use the struct-based variable format described above, in the *QSMB and QSMBCylindrical* subsection.

```
QSMsimple = QSMBCylindrical(ModelData);
```

The class constructor initializes an object of the respective `QSMBCylindrical` class with the data included in the `ModelData` variable. The benefit of using an object, rather than the cell array or a struct, is that a class object also has methods and not just properties. Those methods are utilized inside the leaf insertion function, but also later in the *test.m* file, in the form of the `QSMBCylindrical.export_blender()` method.

```
vertices = [
    -0.3  0.0  0.0;
    -0.3  1.0  0.0;
     0.3  1.0  0.0;
     0.3  0.0  0.0
];

tris = [
     1,  2,  3;...
     1,  3,  4
];
```

These two matrices define the basis geometry of the leaf model. Each row of the `vertices` matrix contains the *x*, *y* and *z* coordinates of a single vertex. Each row of the `tris` matrix contains three indices referencing the rows of the `vertices` matrix, *i.e.*, the three vertices of the respective triangle. Notice that the same vertices can be referenced multiple times if necessary. In this example case there are four vertices and two triangles, as the geometry consists of a rectangle which is 1.0 units tall (*y*-direction) and 0.6 units wide (*x*-direction).

```
LeafArea = [10,20];
```

The two-element vector defines the target leaf area in the first element, and the area of leaf candidates to generate initially in the second element. In this example 20 m<sup>2</sup> of leaf area is generated initially in the form of leaf candidates. Those candidates are accepted one-by-one if they do not intersect with other geometry. When 10 m<sup>2</sup> of leaves are accepted, the process stops. It is also possible that the process is unable to reach the desired leaf area with the generated candidates, thus resulting in a lower final leaf area. If this is the case, there are two main ways to remedy the outcome: 1) increase the candidate leaf area, *i.e.*, the second element of the vector, or 2) increase/modify the transformations performed on leaf candidates when intersections occur, *i.e*, the `'TransformConfig'` parameter of the `qsm_fanni()` function. Note that both solutions are likely to increase the computational time.

```
Leaves = LeafModelTriangle(vertices, tris, {[1 2 3 4]});
```

The variable is initialized as a `LeafModelTriangle` object by the class constructor. The basis geometry defined above is given as the two first inputs. The third, optional input is a cell array that tell the object that the basis geometry can also be described as a single rectangle (*n*-sided polygon), rather than two triangles. This information is not required, but it can be used for optimization, *e.g.*, when exporting resulting leaf covers. The parameter is a cell array, where each cell is a vector of indices of the `vertices` matrix, defining the *n* vertices of the respective polygon.

```
[Leaves, NAccepted] = qsm_fanni(QSMsimple,...
                                Leaves,...
                                LeafArea,...
                                'Seed',1,...
                                'SizeFunctionParameters', {[0.25 0.30]},...
                                'Verbose',true);
```

The main function is run with the variables initialized above (the first three inputs). Note that the second input is the same as the first output, which means that the `Leaves` variable is updated to include the resulting leaf cover parameters. The last three lines of input parameters are optional. As there are random processes in the procedure, fixing the `'Seed'` parameter allows the result of the test run to always be the same. The `'SizeFunctionParameters'` option is used to pass the two-element vector, `[0.25 0.30]`, as an optional argument to the default leaf size distribution function. The argument instructs the function to sample leaf length uniformly inside that interval. Setting the `'Verbose'` option to `true` enables the function to print leaf insertion details during the process. By default verbose printing is off, and the function should not output anything in the console.

When inspecting the `Leaves` variable in the Matlab console, the following printout is given:

```
Leaves = 

  LeafModelTriangle with properties:

       base_vertices: [4×3 double]
     base_dimensions: [0.6000 1 0]
      base_triangles: [2×3 double]
           base_area: 0.6000
          base_ngons: {[1 2 3 4]}
      triangle_count: 2
          leaf_count: 210
           leaf_area: 10.0164
    leaf_start_point: [210×3 double]
          leaf_scale: [210×3 double]
      leaf_direction: [210×3 double]
         leaf_normal: [210×3 double]
         leaf_parent: [1×210 double]
    twig_start_point: [210×3 double]
```

The object has 210 accepted leaves as stated by the `leaf_count` property, that cumulate to a total leaf area of about 10.02 m<sup>2</sup>. As mentioned above, the object does not store the exact geometry of the individual leaves, but rather the transformation parameters that when applied to the basis geometry, give the exact geometry. The user can compute the exact triangle faces and vertices with the `compute_geometry()` method of the `LeafModelTriangle` class: 

```
[Vertices, Faces] = Leaves.compute_geometry(false);
```

The resulting matrices have the same format as the input basis geometry variables `vertices`  and `tris` above. The input parameter is a flag that defines whether to use the *n*-sided polygon rather than triangles in the output faces. The value `false` results in triangles being used.

### Plotting results in Matlab.

```
QSMsimple.plot_cylinders();
hold on;
Leaves.plot_leaves();
hold off;
axis equal;
```

The results are also visualized in Matlab, again by using the provided class methods. Please note that rendering large QSMs and leaf models in Matlab can take a long time.

### Exporting results

```
Leaves.export_geometry('OBJ',true,'test_leaves_export.obj',4);    
QSMsimple.export_blender('test_qsm_export.txt',4);
```

Both the leaf and QSM geometry are exported by class methods provided by the respective classes. Exporting the results to text files is not required. However, this is the process which was used to produce the *test_result.png* image shown also in this document. The exported text files were used to import the geometry into [Blender](https://www.blender.org/) for rendering the image, utilizing the [QSM-blender-addons](https://github.com/InverseTampere/qsm-blender-addons).