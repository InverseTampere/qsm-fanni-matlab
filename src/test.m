%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add QSM-Blocks classes from where you have them.  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Include class definitions.
addpath('classes/');

%% QSM definition.

% Initialize QSM object.
QSMsimple = QSMBCylindrical('example');

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

% Genereate 50 m2 of leaf candidates,
% stop if 10 m2 of leaf area is accepted.
LeafArea = [10,50];

% Initialize the leaf model with the basis geometry.
Leaves = LeafModelTriangle(vertices, tris, {[1 2 3 4]});

% Generate leaves.
[Leaves, NAccepted] = qsm_fanni( ...
    QSMsimple,...
    Leaves,...
    LeafArea,...
    'Seed',1,...
    'SizeFunctionParameters', {[0.25 0.30]},...
    'Verbose',true ...
);

% Optional steps to follow:

%% Plot results.

% Plot QSM.
hQSM = QSMsimple.plot_model();
% Set bark color.
set(hQSM,'FaceColor',[150,100,50]./255,'EdgeColor',[0 0 0]);

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
