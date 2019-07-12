# QSM-FaNNI

[![DOI](https://zenodo.org/badge/92806841.svg)](https://zenodo.org/badge/latestdoi/92806841)

Quantitative structure models - Foliage and needles naive insertion algorithm MATLAB implementation. Initially written and tested with Matlab version R2016b and later updated and tested with version R2018a. See some of the customization possibilities in action, together with example results in an [animation](https://www.youtube.com/watch?v=urPDwcEf02A).

![Test result](https://github.com/InverseTampere/qsm-fanni-matlab/raw/master/src/test_result.png)

## Description

QSM-FaNNI can be used to generate a leaf cover for a quantitative structure model (QSM) of a tree. The shape of the leaves, as well as, their location, orientation and size can be controlled by the user. All leaves have identical basis geometry - currently the geometry can consist of any number of triangular faces - that is then manipulated by scaling, translation and rotation. 

The distribution of leaf material on the QSM is controlled by the leaf area density distribution (LADD), leaf size by the leaf size distribution (LSD) and leaf orientation by the leaf orientation distribution (LOD). Furthermore, the position, orientation and length of the petioles connecting the leaves to the QSM are controlled by additional distributions. All of the mentioned distributions are Matlab functions with a fixed interface (certain inputs and outputs). Default functions are included by the user is encouraged to create their own.

The program generates candidate leaves that are accepted to the final leaf cover if they do not intersect the with QSM geometry or any accepted leaves. If an intersection occurs the program tries to modify the leaf parameters (location, scale and/or orientation) with any number of user-customizable transformations, before finally discarding the candidate if intersections still persist.

## Dependencies

### QSMB and QSMBCylindrical

Quantitative structure models (QSMs) are expected to be defined as objects of QSMB subclasses, such as QSMBCylindrical. These Matlab classes are part of the [QSM-Blocks](https://github.com/InverseTampere/qsm-blocks-matlab) repository.

### CubeVoxelization

Internally the leaf generation procedure relies on partitioning the space in to voxels. The CubeVoxelization class handles partitioning and related tasks. The CubeVoxelization class is part of the [QSM-Blocks](https://github.com/InverseTampere/qsm-blocks-matlab) repository.

## Basic usage

The program relies on a number of classes that have to be added to the current path before running QSM-FaNNI.

```Matlab
addpath('classes/');
```

The main program the function `qsm_fanni` defined in *qsm_fanni.m*. The function receives the QSM as a `QSMB` object (or a subclass object), an initialized `LeafModel` object (containing the leaf basis geometry), and total leaf area to be distributed. The leaf area parameter can have two components; one for the initial leaf area to be generated, and one for the target leaf area. This can be used to increase the probability that the target area is reached, even if some of the generated leaves are discarded due to unavoidable intersections.

The main function call has the following formats:

```Matlab
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

```Matlab
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

The abstract class QSMB and the class QSMBCylindrical are used to hold geometric and topological quantitative structure model data. As of version 1.3.0 these are part of a separate [QSM-Blocks](https://github.com/InverseTampere/qsm-blocks-matlab) repository, but are required dependencies.

### LeafModel and LeafModelTriangle

The `LeafModel` class contains the minimum API any type of leaf model must provide. The class is abstract, and thus, for computations the user must initialize an `LeafModelTriangle` object for leaves with the basis geometry consisting of triangles, or create their own subclass of `LeafModel` for any other type of models. The class constructor can take the following formats:

```Matlab
Leaves = LeafModelTriangle(Vert, Faces);
Leaves = LeafModelTriangle(Vert, Faces, Ngons);
Leaves = LeafModelTriangle(Vert, Faces, Ngons, NInit);
```

where `Vert` is a Nx3 matrix with the vertices as rows, `Faces` is a Mx3 matrix with vertex indices of a single triangle on each row. The optional argument `Ngons` is a cell array, there each element contains a vector of vertex indices forming *n*-sided polygons, used for optimization, e.g., when exporting the leaf model. `NInit` is a scalar number that can be used to initialize the objects internal matrices to optimize leaf insertion. After leaf insertion has been completed, the extra rows can be removed with the `trim_slack()` method (that is called by `qsm_fanni` when leaf insertion ends):

```Matlab
Leaves.trim_slack();
```

Internally the `LeafModelTriangle` object stores a single copy of the basis geometry and a set of transformation parameters (translation, scaling, rotation) per each accepted leaf. Should the exact geometry of individual leaves be required, the class offers two methods for computing geometry:

```Matlab
Tris =  Leaves.triangles(origin, dir, normal, scale);
[Vertices, Faces] =  Leaves.triangles(origin, dir, normal, scale);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon, Filter);
[Vertices, Faces] =  Leaves.compute_geometry(fNgon, origin, dir, normal, scale);
```

where `origin`, `dir`, `normal`, `scale` are single transformation parameters, `fNgon` is a flag determining whether to use *n*-sided polygon rather than triangles, when available. Furthermore, `Filter` is a logical vector with a length equal to the number of accepted leaves, and where the value `true` filters in the corresponding leaf. `Tris` is a N x 9 matrix with the following structure:

```Matlab
[v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z]
```

where `v1`, `v2` and `v3` are the three vertices of a triangle.

## The test file

The *src* folder contains a simple test file, *test.m* for trying the leaf insertion procedure with example data. A detailed description of the command included in the compact demonstration is given below, with the relevant commands listed before the respective description.

```Matlab
addpath('classes/');
```

The command includes the required class definitions into Matlab's path. These classes include `LeafModel` and `LeafModelTriangle`. As of QSM-FaNNI version 1.3.0 you have to manually also add the path of the class definitions in [QSM-Blocks](https://github.com/InverseTampere/qsm-blocks-matlab).

In total the dependencies include the four classes described above, as well as, a helper class `CubeVoxelization` that is used internally by the main program to describe objects, such as triangles or cylinders, in three-dimensional voxelized spaces.

```Matlab
QSMsimple = QSMBCylindrical('example');
```

The class constructor initializes an object of the respective `QSMBCylindrical` class with example data. For additional means for initializing a QSMBCylindrical object, see the QSM-Blocks *README.md*. The benefit of using an object, rather than the cell array or a struct, is that a class object also has methods and not just properties. Those methods are utilized inside the leaf insertion function, but also later in the *test.m* file, in the form of the `QSMBCylindrical.export_blender()` method.

```Matlab
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

```Matlab
LeafArea = [10,50];
```

The two-element vector defines the target leaf area in the first element, and the area of leaf candidates to generate initially in the second element. In this example 50 m<sup>2</sup> of leaf area is generated initially in the form of leaf candidates. Those candidates are accepted one-by-one if they do not intersect with other geometry. When 10 m<sup>2</sup> of leaves are accepted, the process stops. It is also possible that the process is unable to reach the desired leaf area with the generated candidates, thus resulting in a lower final leaf area. If this is the case, there are two main ways to remedy the outcome: 1) increase the candidate leaf area, *i.e.*, the second element of the vector, or 2) increase/modify the transformations performed on leaf candidates when intersections occur, *i.e*, the `'TransformConfig'` parameter of the `qsm_fanni()` function. Note that both solutions are likely to increase the computational time.

```Matlab
Leaves = LeafModelTriangle(vertices, tris, {[1 2 3 4]});
```

The variable is initialized as a `LeafModelTriangle` object by the class constructor. The basis geometry defined above is given as the two first inputs. The third, optional input is a cell array that tell the object that the basis geometry can also be described as a single rectangle (*n*-sided polygon), rather than two triangles. This information is not required, but it can be used for optimization, *e.g.*, when exporting resulting leaf covers. The parameter is a cell array, where each cell is a vector of indices of the `vertices` matrix, defining the *n* vertices of the respective polygon.

```Matlab
[Leaves, NAccepted] = qsm_fanni( ...
    QSMsimple,...
    Leaves,...
    LeafArea,...
    'Seed',1,...
    'SizeFunctionParameters', {[0.25 0.30]},...
    'Verbose',true ...
);
```

The main function is run with the variables initialized above (the first three inputs). Note that the second input is the same as the first output, which means that the `Leaves` variable is updated to include the resulting leaf cover parameters. The last three lines of input parameters are optional. As there are random processes in the procedure, fixing the `'Seed'` parameter allows the result of the test run to always be the same. The `'SizeFunctionParameters'` option is used to pass the two-element vector, `[0.25 0.30]`, as an optional argument to the default leaf size distribution function. The argument instructs the function to sample leaf length uniformly inside that interval. Setting the `'Verbose'` option to `true` enables the function to print leaf insertion details during the process. By default verbose printing is off, and the function should not output anything in the console.

When inspecting the `Leaves` variable in the Matlab console, the following printout is given:

```Matlab
Leaves = 

  LeafModelTriangle with properties:

       base_vertices: [4×3 double]
     base_dimensions: [0.6000 1 0]
      base_triangles: [2×3 double]
           base_area: 0.6000
          base_ngons: {[1 2 3 4]}
      triangle_count: 2
          leaf_count: 209
           leaf_area: 10.0173
    leaf_start_point: [209×3 double]
          leaf_scale: [209×3 double]
      leaf_direction: [209×3 double]
         leaf_normal: [209×3 double]
         leaf_parent: [209×1 double]
    twig_start_point: [209×3 double]
```

The object has 210 accepted leaves as stated by the `leaf_count` property, that cumulate to a total leaf area of about 10.02 m<sup>2</sup>. As mentioned above, the object does not store the exact geometry of the individual leaves, but rather the transformation parameters that when applied to the basis geometry, give the exact geometry. The user can compute the exact triangle faces and vertices with the `compute_geometry()` method of the `LeafModelTriangle` class: 

```Matlab
[Vertices, Faces] = Leaves.compute_geometry(false);
```

The resulting matrices have the same format as the input basis geometry variables `vertices`  and `tris` above. The input parameter is a flag that defines whether to use the *n*-sided polygon rather than triangles in the output faces. The value `false` results in triangles being used.

### Plotting results in Matlab.

```Matlab
hQSM = QSMsimple.plot_cylinders();
set(hQSM,'FaceColor',[150,100,50]./255,'EdgeColor',[0 0 0]);
hold on;
hLeaf = Leaves.plot_leaves();
set(hLeaf,'FaceColor',[120,150,80]./255);
hold off;
axis equal;
```

The results are also visualized in Matlab, again by using the provided class methods. Please note that rendering large QSMs and leaf models in Matlab can take a long time.

### Exporting results

```Matlab
% Use ngons when exporting leaves.
fUseNgon = true;

% Export in OBJ-format with individual leaf vertices and faces.
Leaves.export_geometry( ...
    'OBJ', ...
    fUseNgon, ...
    'test_leaves_export.obj', ...
    4 ...
);

% Export in custom extended OBJ-format with basis leaf geometry 
% and individual leaf transformation parameters.
Leaves.export_geometry( ...
    'EXT_OBJ', ...
    fUseNgon, ...
    'test_leaves_export_extended.obj', ...
    4 ...
);

% Export QSM parameters to a text file.
QSMsimple.export( ...
    'blender', ...
    'test_qsm_export.txt', ...
    'Precision',4 ...
);
```

Both the leaf and QSM geometry are exported by class methods provided by the respective classes. Exporting the results to text files is not required. However, this is the process which was used to produce the *test_result.png* image shown also in this document. The exported text files were used to import the geometry into [Blender](https://www.blender.org/) for rendering the image, utilizing the [QSM-blender-addons](https://github.com/InverseTampere/qsm-blender-addons). 

As of QSM-FaNNI version 1.2.0, there are tow different export formats for leaf model exportation. With the Wavefront OBJ-format (`'OBJ'`) the geometry of individual leaves, *i.e.* vertices and faces, is exported. With the custom *Extended OBJ*-format (`'EXT_OBJ'`) the leaf basis geometry is exported in standard OBJ-format, but in addition individual leaf transformation parameters are exported as well. The leaf definition lines will have the following parameters:
1. (0) 'L' (line type definition)
2. (1-3) Twig start point
3. (4-6) Leaf start point
4. (7-9) Leaf direction (start point to leaf tip)
5. (10-12) Leaf normal
6. (13-15) Leaf scale
7. (16-N) Optional parameters, such as leaf color

The benefit of the extended format is to allow growth animation and leaf level coloring when importing leaf data into Blender with the [QSM-blender-addons](https://github.com/InverseTampere/qsm-blender-addons). 
