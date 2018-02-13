%% Include class definitions.
addpath('classes/');

%% QSM definition.


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

% Cell array format of tree model data.
ModelData = {CylData, BranchData, TreeData};

% Initialize QSM object.
QSMsimple = QSMBCylindrical(ModelData);

%% Define leaf shape.

% Vertices of the leaf basis geometry.
vertices = [
    -0.3  0.0  0.0;
    -0.3  1.0  0.0;
     0.3  1.0  0.0;
     0.3  0.0  0.0
];

% Triangles of the leaf basis geometry.
tris = [
     1,  2,  3;...
     1,  3,  4
];

%% Leaf insertion.

% Genereate 20 m2 of leaf candidates,
% stop if 10 m2 of leaf area is accepted.
LeafArea = [10,20];

% Initialize the leaf model with the basis geometry.
Leaves = LeafModelTriangle(vertices, tris, {[1 2 3 4]});

% Generate leaves.
[Leaves, NAccepted] = qsm_fanni(QSMsimple,...
                                Leaves,...
                                LeafArea,...
                                'Seed',1,...
                                'SizeFunctionParameters', {[0.25 0.30]},...
                                'Verbose',true);
%-

% Optional steps to follow:

%% Plot results.

% Plot QSM.
hQSM = QSMsimple.plot_cylinders();
% Set bark color.
set(hQSM,'FaceColor',[150,100,50]./255);

hold on;

% Plot leaves.
hLeaf = Leaves.plot_leaves();
% Set leaf color.
set(hLeaf,'FaceColor',[120,150,80]./255);

hold off;
axis equal;

%% Compute geometry of accepted leaves.

[Vertices, Faces] = Leaves.compute_geometry(false);

%% Export result.

% Use ngons when exporting leaves.
fUseNgon = true;

% Export in OBJ-format with individual leaf vertices and faces.
Leaves.export_geometry('OBJ',fUseNgon,'test_leaves_export.obj',4);

% Export in custom extended OBJ-format with basis leaf geometry 
% and individual leaf transformation parameters.
Leaves.export_geometry('EXT_OBJ',fUseNgon,'test_leaves_export_extended.obj',4);

% Export QSM parameters to a text file.
QSMsimple.export_blender('test_qsm_export.txt',4);
