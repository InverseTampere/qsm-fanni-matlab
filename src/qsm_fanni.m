% Generate random leaves for a given cylinder model. <QSM> contains
% the structure model, <Leaves> is the leaf model to be used (shape etc.),
% and <Area> is the total leaf area to be distributed to the QSM.
%
% [Leaves,Twigs,N,Parent] = qsm_fanni(QSM,Leaves,Area);
%                  Leaves = qsm_fanni(QSM,Leaves,Area,...);
%
% =================== Outputs =============================================
%
% Leaves              Class LeafModel object that contains the added
%                     leaves, and corresponding twig start points.
%
% NAccepted           Final number of generated leafs.
%
% 'Others'            The other optional outputs are for debuging.
%
% =================== Inputs ==============================================
%
% QSM                 Object that contains structure information.
%
% Leaves              Object that defines leaf geometry and will store leaf
%                     parameters after insertion.
%
% Area                Target leaf area. If the variable has two components,
%                     then Area(1) is the target area, and Area(2) is the
%                     initial area to be generated.
%
% Several other parameters can also be given in name-value pairs
% as explained below. The names are not case-sensitive.
%
% 'TransformConfig'   Matrix that contains the scaling, rotation and
%                     translation configurations to apply one-by-one when a
%                     leaf intersects with other leaves or blocks. Each
%                     row of the matrix is a single configuration, columns
%                     are as follow: 
%                          1) scale width  (x, side dir.)
%                          2) scale length (y, front dir.)
%                          3) scale height (z, normal dir.)
%                          4) rotate around x (side)
%                          5) rotate around y (front)
%                          6) rotate around z (normal)
%                          7) translate x (side)
%                          8) translate y (front)
%                          9) translate z (normal)
%
%                     Default value:
%                          x    y    z     x     y    z    x    y    z
%                          -------------------------------------------
%                          1    1    1     0  pi/6    0    0    0    0
%                          1    1    1     0 -pi/6    0    0    0    0
%                          1    1    1     0  pi/2    0    0    0    0
%                          1    1    1  pi/6     0    0    0    0    0
%                          1    1    1 -pi/6     0    0    0    0    0
%                          1    1    1 -pi/6  pi/2    0    0    0    0
%                          1    1    1  pi/6  pi/2    0    0    0    0
%
% 'RandomTransform'   When set to true transformation configuration are
%                     tried in random order rather than row-by-row.
%                     Default: false
%
% 'AreaFunction'      A handle to a function that distributes available
%                     total volume to blocks.
%
%                     Inputs:
%                          1) Number of blocks.
%                          2) Parameters of each block (struct).
%                          3) Optional extra parameters.
%
%                     Outputs:
%                          1) Vector of relative leaf area in each block.
%                             Should sum up to one.
%
% 'SizeFunction'      A handle to a function that determines the number of
%                     leaves per block, and the size of each leaf, i.e.,
%                     transforms area per block to leaf sizes. The function
%                     should have the following inputs and outputs:
%
%                     Inputs:
%                          1) Vector of area for each block.
%                          2) Parameters of each block (struct).
%                          3) Area of leaf base.
%                          4) Dimensions of leaf base.
%                          5) Optional extra parameters.
%
%                     Outputs:
%                          1) Matrix of leaf dimensions (NLeaf x 3).
%                          2) Vector of leaf parents, i.e., block ids.
%                          3) Maximum absolute dimension of any leaf, used
%                             for voxelization.
%
% 'AngleFunction '    A handle to a function that determines the direction
%                     and normal of each leaf. 
%
%                     Inputs:
%                          1) Parameters of each block (struct).
%                          2) Parent index of each leaf (block id).
%                          3) Direction of the twig of each leaf.
%                          4) Optional extra parameters.
%
%                     Outputs:
%                          1) Leaf directions (NLeaf x 3).
%                          2) Leaf normals (NLeaf x 3).
%
% 'AreaFunctionParameters'     Cell array of extra parameters to pass to
%                              the AreaFunction.
%
% 'SizeFunctionParameters'     Cell array of extra parameters to pass to
%                              the SizeFunction.
%
% 'AngleFunctionParameters'    Cell array of extra parameters to pass to
%                              the AngleFunction.
%
% 'Seed'              Set the seed for the random number generator before
%                     the leaf generation. The values is used in 
%                     reset(RandStream.getGlobalStream,<value>);
%
% 'Verbose'           Boolean that controls printing information on 
%                     the process status.


% This file is part of QSM-FaNNI.
% 
% QSM-FaNNI is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QSM-FaNNI is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QSM-FaNNI.  If not, see <http://www.gnu.org/licenses/>.

function [Leaves, NAccepted, IAccepted, NConfigsTried, NNeighbour] ...
                           = qsm_fanni(QSM,Leaves,Area,varargin)

%% Inputs and default values.

% Leaf target and possible over-sampling area.
if length(Area) == 1
    AreaInsert = Area;
    AreaTarget = Area;
else
    AreaInsert = Area(2);
    AreaTarget = Area(1);
end

if AreaInsert < AreaTarget
    warning( ...
        'Leaf candidate area (%g) is less than leaf target area (%g).', ...
        AreaInsert, ...
        AreaTarget ...
    );
end

% Matrix of transform configurations.
% 1: scale X
% 2: scale Y
% 3: scale Z
% 4: rotate around X (side)
% 5: rotate around Y (dir
% 6: rotate around Z (normal)
% 7: translate X (side)
% 8: translate Y (dir)
% 9: translate Z (normal)
%TransformConfig = zeros(0,9);
TransformConfig = [...
%   x    y    z     x     y    z    x    y    z
    1    1    1     0  pi/6    0    0    0    0;
    1    1    1     0 -pi/6    0    0    0    0;
    1    1    1     0  pi/2    0    0    0    0;
    1    1    1  pi/6     0    0    0    0    0;
    1    1    1 -pi/6     0    0    0    0    0;
    1    1    1 -pi/6  pi/2    0    0    0    0;
    1    1    1  pi/6  pi/2    0    0    0    0;
];

% Flag: randomize transform order.
RandomizeTransform = false;

%%% Default distribution functions.

% By default all cylinder receive equal portition of leaf area.
fun_area = @default_fun_area;
fun_size = @default_fun_size;
fun_angle = @default_fun_angle;

%%% Default extra parameters for distribution functions.
extra_area = cell(0);
extra_size = cell(0);
extra_angle = cell(0);


% Printing on/off.
verbose = false;

%% Read inputs.

% Check additional parameters.
i = 1;
NArg = numel(varargin);
while i <= NArg

    if ischar(varargin{i})

        switch lower(varargin{i})
            
            case 'transformconfig'
                assert(i < NArg && ...
                       isnumeric(varargin{i+1}) && ...
                       size(varargin{i+1},2) == 9, ...
                       'Argument following ''TransformConfig'' should be a (N x 9) matrix.');
                TransformConfig = varargin{i+1};
                i = i + 1;
                
            case 'randomtransform'
                assert(i < NArg && ...
                       islogical(varargin{i+1}), ...
                       'Argument following ''RandomTransform'' should be boolean.');
                RandomizeTransform = varargin{i+1};
                i = i + 1;
                
            case 'areafunction'
                if i == NArg || not(isa(varargin{i+1}, 'function_handle'))
                    error('Argument following ''AreaFunction'' should be a function handle.');
                else
                    fun_area = varargin{i+1};
                    i = i + 1;
                end
                
            case 'sizefunction'
                if i == NArg || not(isa(varargin{i+1}, 'function_handle'))
                    error('Argument following ''SizeFunction'' should be a function handle.');
                else
                    fun_size = varargin{i+1};
                    i = i + 1;
                end
                
            case 'anglefunction'
                if i == NArg || not(isa(varargin{i+1}, 'function_handle'))
                    error('Argument following ''AngleFunction'' should be a function handle.');
                else
                    fun_angle = varargin{i+1};
                    i = i + 1;
                end
                
            case 'areafunctionparameters'
                assert(i < NArg && ...
                       iscell(varargin{i+1}), ...
                       'Argument following ''AreaFunctionParameters'' should be a cell array.');
                extra_area = varargin{i+1};
                i = i + 1;
                
            case 'sizefunctionparameters'
                assert(i < NArg && ...
                       iscell(varargin{i+1}), ...
                       'Argument following ''SizeFunctionParameters'' should be a cell array.');
                extra_size = varargin{i+1};
                i = i + 1;
                
            case 'anglefunctionparameters'
                assert(i < NArg && ...
                       iscell(varargin{i+1}), ...
                       'Argument following ''AngleFunctionParameters'' should be a cell array.');
                extra_angle = varargin{i+1};
                i = i + 1;
            
            case 'seed'
                assert(i < NArg && ...
                       isnumeric(varargin{i+1}), ...
                       isscalar(varargin{i+1}), ...
                       'Argument following ''Seed'' should be an integer.');
                %reset(RandStream.getGlobalStream,varargin{i+1});
                rng(varargin{i+1},'twister');
                
            case 'verbose'
                assert(i < NArg && ...
                       islogical(varargin{i+1}), ...
                       'Argument following ''Verbose'' should be boolean.');
                verbose = varargin{i+1};
                i = i + 1;
                
            otherwise
                warning(['Skipping unknown parameters: ''' varargin{i} '''']);

        end
    end
    i = i + 1;
end

%% Initialize rest of parameters.

% Number of transforms.
NTransform = size(TransformConfig,1);

%% Cylinder parameters.

% Paramenter labels.
ParLabels = {
                'relative_height',...
                'relative_position',...
                'branch_order',...
                'radius',...
                'is_last'
            };
%-

BlocksParameters = QSM.get_block_properties(ParLabels);


%% Sample leaf parameters: dimensions and orientation.

% Compute percentages for each cylinder.
RelativeBlockArea = fun_area(QSM.block_count,...
                             BlocksParameters,...
                             extra_area{:});
%-

% Scale by total area to get area per block.
BlockArea = AreaInsert.*RelativeBlockArea;

% Sample size distribution to get count and dimensions of leaves
% attached to the block.
[LeafDimensions, LeafParent, MaxLeafSize] = fun_size(BlockArea,...
                                                 BlocksParameters,...
                                                 Leaves.base_area,...
                                                 Leaves.base_dimensions,...
                                                 extra_size{:});
%-

% Compute total number of candidate leaves.
NLeafCandidate = length(LeafParent);


% Attach leaves to blocks.
% Compute twigs, i.e., connection points and lengths of twigs.
[TwigStart, TwigEnd] = QSM.generate_twigs(LeafParent);

TwigDir = TwigEnd - TwigStart;
TwigLen = sqrt(sum(TwigDir.^2,2));
MaxTwigLen = max(TwigLen);
TwigDir = bsxfun(@times,TwigDir,1./TwigLen);

% Sample leaf angle distribution.
[LeafDir, LeafNormal] = fun_angle(BlocksParameters,...
                                  LeafParent,...
                                  TwigDir,...
                                  extra_angle{:});
%-

if verbose
    disp('Leaf insertion started.')
    fprintf('\tGenerated %d leaves.\n',NLeafCandidate);
end
        

%% Voxelization parameters.

% Extreme values of QSM.
TreeBox = QSM.tree_limits;

MinPoint = TreeBox(1,:) - MaxLeafSize - MaxTwigLen;
MaxPoint = TreeBox(2,:) + MaxLeafSize + MaxTwigLen;

% Compute voxelization of QSM.
QSMVoxelization = QSM.toVoxels(MaxLeafSize, MinPoint, MaxPoint);

% Initialize voxelization of accepted leaves.
LeafVoxelization = CubeVoxelization(MaxLeafSize, MinPoint, MaxPoint);


%% Randomize order.

% Generate random order.
JRand = randperm(NLeafCandidate);

% Reorder matrices' rows.
TwigStart      =      TwigStart(JRand,:);
TwigEnd        =        TwigEnd(JRand,:);
LeafDir        =        LeafDir(JRand,:);
LeafNormal     =     LeafNormal(JRand,:);
LeafDimensions = LeafDimensions(JRand,:);
LeafParent     =     LeafParent(JRand);

% Acceptance vector.
IAccepted = false(NLeafCandidate,1);

% Number of configs tried.
NConfigsTried = nan(NLeafCandidate,1);

% Number of neighbour on initial run.
% 1: Block, 2: Leaf
NNeighbour = nan(NLeafCandidate,2);

% Number of accepted leaves.
NAccepted = 0;

% Total area of accepted leaves.
AreaAccepted = 0;

JTransform = 1:NTransform;

fTargetReached = false;

% Start adding leaves to the leaf model one-by-one.
timerVal = tic;
for iLeaf = 1:NLeafCandidate
    elapsedTime = toc(timerVal);
    progress = round(AreaAccepted/AreaTarget*100,1);
    if iLeaf == 1
        f = waitbar(0,'Starting leaf insertion');
    else 
        waitbar(progress/100,f,sprintf(['Leaf insertion %3.1f%% done,' ...
        ' elapsed time: %3.1f min\nTime remaining estimate: %3.1f min'], ...
        [progress,elapsedTime/60,100/progress*elapsedTime/60-elapsedTime/60]))
    end         
    % If target area reached stop adding leaves.
    if AreaAccepted >= AreaTarget
        fTargetReached = true;
        break;
    end
    
    % Randomize transform config order, if neccessary.
    if RandomizeTransform
        JTransform = randperm(NTransform);
    end
    
    % Cube coordinates of leaf on last try.
    LastCC = [];
    
    % Compute leaf neighbour at least once.
    ComputeLeafNeigbours = true;
    
    for iTransform = 1:NTransform+1
        
        % Paramenters of the current leaf.
        origin = TwigEnd(iLeaf,:);
        dir    = LeafDir(iLeaf,:);
        normal = LeafNormal(iLeaf,:);
        scale  = LeafDimensions(iLeaf,:);
        
        % Third axis (x) of leaf coordinate system.
        side = cross(dir,normal);
        
        % Flag: neighbour index computation required.
        ComputeNeighbours = false;
        
        % Do transformation.
        if iTransform > 1
            
            % Transformation configuration.
            Config = TransformConfig(JTransform(iTransform-1),:);
            
            % Do transformation.
            
            % Scaling: multiply current dimensions.
            if any(Config(1:3) ~= ones(1,3))
                scale = scale.*Config(1:3);
            end
            
            % Elevation / around side axis.
            if Config(4)
                Rx = rotation_matrix(side,Config(4));
                dir = (Rx*dir')';
                normal = (Rx*normal')';
            end
            
            % Rotation / around direction.
            if Config(5)
                Ry = rotation_matrix(dir,Config(5));
                %side = (Ry*side')';
                normal = (Ry*normal')';
            end
            
            % Azimuth / around normal.
            if Config(6)
                Rz = rotation_matrix(normal,Config(6));
                dir = (Rz*dir')';
                %side = (Rz*side')';
            end
            
            % Translation.
            if any(Config(7:9) ~= zeros(1,3))
                origin = origin + Config(7:9);
            end
        end
        
        % Compute leaf center point.
        cen = origin + scale(2)*dir;
        
        % Compute leaf triangles for intersection computations.
        LeafTris = Leaves.triangles(origin, dir, normal, scale);
        
        % Convert center point to cube voxelization coordinates.
        LeafCC = LeafVoxelization.get_coordinates(cen);
    
        % On first run neighbour computations are always required.
        if iTransform == 1
            ComputeNeighbours = true;
            LastCC = LeafCC;
            
        % Otherwise check if cube coordinate of leaf has changed.
        elseif any(LeafCC ~= LastCC)
            ComputeNeighbours = true;
            LastCC = LeafCC;
        end
    
        % Find neighbour blocks and their count only if changes have
        % occured.
        if ComputeNeighbours
            % Get neighbour block indices. Exclude parent block.
            BlockNeighbour = QSMVoxelization.get_neighbor_objects(LeafCC);

            % Number of block neighbours.
            NBlockNei = length(BlockNeighbour);
            %%%disp(['Block neighbors: ' num2str(NBlockNei)]);
            
            if isnan(NNeighbour(iLeaf,1))
                NNeighbour(iLeaf,1) = NBlockNei;
            end
        end

        % If any overlap, try next configuration.
        if NBlockNei && QSM.block_triangle_intersection(BlockNeighbour,...
                                                        LeafTris)
            %-
            
            %%%disp(['Discarded leaf (block inter.): ' num2str(iLeaf)]);
            continue;
        end
        
        % Find neighbour blocks and their count only if changes have
        % occured.
        if ComputeNeighbours || ComputeLeafNeigbours
            % Get neighbour block indices. Exclude parent block.
            LeafNeighbour = LeafVoxelization.get_neighbor_objects(LeafCC);
            
            % After first compute, respect <ComputeNeighbours>.
            ComputeLeafNeigbours = false;

            % Number of block neighbours.
            NLeafNei = length(LeafNeighbour);
            %%%disp(['Leaf neighbors: ' num2str(NLeafNei)]);
            
            if isnan(NNeighbour(iLeaf,2))
                NNeighbour(iLeaf,2) = NLeafNei;
            end
        end

        % If any overlap, try next configuration.
        if NLeafNei && Leaves.leaf_intersect(LeafNeighbour,LeafTris)
            %%%disp(['Discarded leaf (leaf inter.): ' num2str(iLeaf)]);
            continue;
        end
    
        % Parent of accepted leaf.
        parent = LeafParent(iLeaf);
        
        % Twig start point.
        twig = TwigStart(iLeaf,:);
        
        % No intersections: leaf is accepted and inserted into model.
        [LeafIndex,LeafArea] = Leaves.add_leaf(origin,...
                                               dir,...
                                               normal,...
                                               scale,...
                                               parent,...
                                               twig,...
                                               LeafTris);
        %-
        
        % Add leaf to voxelization.
        LeafVoxelization.add_object_by_cc(LeafCC, LeafIndex);
        
        % Mark success.
        IAccepted(iLeaf) = true;
        
        % Increase count and total area.
        NAccepted = NAccepted + 1;
        AreaAccepted = AreaAccepted + LeafArea;
        
        % Skip remaining transform configurations.
        break;
        
    end
    
    % Store applied transformation config id.
    NConfigsTried(iLeaf) = (iTransform - 1);
    
end

if fTargetReached
    iLeaf = iLeaf - 1;
end

if verbose
    if fTargetReached
        fprintf('\tTarget leaf area reached.\n');
    else
        fprintf('\tTarget area missed by %.2f.\n',...
                AreaTarget - AreaAccepted);
        %-
    end
    fprintf('\tAccepted %d leaves.\n',NAccepted);
    fprintf('\tAccepted %.2f%% of tried leaves.\n',...
            100*nnz(IAccepted(1:iLeaf))/iLeaf);
    %-
    fprintf('Leaf insertion finished.\n')
end

% Prepare accepted flags.
if nargout > 2

    IAccepted = IAccepted(1:iLeaf);

end

% Prepare transformation config info.
if nargout > 3

    NConfigsTried = NConfigsTried(1:iLeaf);

end

% Prepare neighbor info.
if nargout > 4

    NNeighbour = NNeighbour(1:iLeaf,:);

end

% Clear possible empty rows.
Leaves.trim_slack();

end


function R = rotation_matrix(u,k)
% Compute rotation around axis <u> by <k> radians.

    R = zeros(3,3);
    c = cos(k); 
    s = sin(k);
    R(1,:) = [u(1)^2+(1-u(1)^2)*c,    ...
              u(1)*u(2)*(1-c)-u(3)*s, ...
              u(1)*u(3)*(1-c)+u(2)*s];
    %-
    R(2,:) = [u(1)*u(2)*(1-c)+u(3)*s, ...
              u(2)^2+(1-u(2)^2)*c,    ...
              u(2)*u(3)*(1-c)-u(1)*s];
    %-
    R(3,:) = [u(1)*u(3)*(1-c)-u(2)*s, ...
              u(2)*u(3)*(1-c)+u(1)*s, ...
              u(3)^2+(1-u(3)^2)*c];
    %-

end

% Default function for leaf area density distribution.
function RelativeBlockArea = default_fun_area(NBlock, ...
                                              BlocksParameters,...
                                              varargin)
    %-
    
    % Flag: true if block is the last in its branch, false otherwise.
    is_last = BlocksParameters.is_last;
    
    
    if any(is_last)
        % Distribute evenly to all branch tip blocks.
        RelativeBlockArea = zeros(NBlock,1);
        RelativeBlockArea(is_last) = 1/nnz(is_last);
    else
        % Fallback, if none of the blocks are marked to be last. Distribute
        % evenly to all blocks.
        RelativeBlockArea = ones(NBlock,1)/NBlock;
    end
    

end

% Default function for leaf size distribution.
function [LeafDimensions, LeafParent, MaxLeafSize] = ...
                              default_fun_size(BlockAreaTarget,...
                                               BlocksParameters,...
                                               BaseArea,...
                                               ~,... % BaseDim
                                               varargin)
    %-

    % Relative height of blocks.
    relheight = BlocksParameters.relative_height;
    
    % Height limits.
    minh = min(relheight);
    maxh = max(relheight);
    
    % How much size scaling based on vertical positioning.
    HeightMaxVar = 0.25;
    
    % If all cylinders are on the same vertical level.
    if minh == maxh
        % Constant scaling.
        HeightPoly = polyfit(0,1,0);
    else
        % Linear scaling based on height.
        HeightPoly = polyfit([minh,maxh],...
                             [1-HeightMaxVar, 1+HeightMaxVar],...
                             1);
        %-
    end
    
    % Number of blocks.
    NBlock = length(BlockAreaTarget);
    
    % Generate random order.
    JRand = randperm(NBlock);

    % Sampling limits for leaf Y-direction.
    LeafLengthLimits = [0.04 0.06];
    
    if numel(varargin)
        LeafLengthLimits = varargin{1};
    end

    % Average leaf length.
    MeanLeafLength = mean(LeafLengthLimits);

    % Lower limit for block area. If block area difference is below
    % this value, sample another leaf.
    AreaLimit = BaseArea*MeanLeafLength*0.5;
    
    % Area distributed to each block.
    BlockArea = zeros(NBlock,1);

    % Estimate the number of leaves to be sampled.
    NLeafEstimate = sum(BlockArea)/(BaseArea*MeanLeafLength);
    NLeafEstimate = floor(1.2*NLeafEstimate);

    % Initialize outputs.
    LeafDimensions = zeros(NLeafEstimate,3);
    LeafParent = zeros(NLeafEstimate,1);

    % Number of generated leaves.
    NLeaf = 0;

    for iBlock = 1:NBlock
        
        % How much total leaf area should be distributed.
        CurrentLeafAreaTarget = sum(BlockAreaTarget(JRand(1:iBlock)));
        
        % How much total leaf area has been distributed.
        TotalLeafArea = sum(BlockArea(JRand(1:iBlock)));

        % While the block in question either:
        % 1) has personal capacity, or
        % 2) total difference is negative,
        % sample another leaf.
        while BlockAreaTarget(JRand(iBlock)) - BlockArea(JRand(iBlock)) ...
              >= AreaLimit || ...
              CurrentLeafAreaTarget - TotalLeafArea >= AreaLimit
            %-
            
            % Increase leaf count.
            NLeaf = NLeaf + 1;

            % Sample the length of the leaf.
            LeafLength = LeafLengthLimits(1) ...
                       + (LeafLengthLimits(2) - LeafLengthLimits(1)) ...
                       * rand(1) * polyval(HeightPoly,...
                                           relheight(JRand(iBlock)));
            %-

            % Area of new leaf is computed by scaling the area of the base.
            LeafArea = BaseArea*(LeafLength^2);

            % Update area difference in block.
            BlockArea(JRand(iBlock)) = BlockArea(JRand(iBlock)) + LeafArea;
            
            % Update cumulative area difference.
            TotalLeafArea = TotalLeafArea + LeafArea;

            % Store dimensions that are computed based on sampled length.
            LeafDimensions(NLeaf,:) = repmat(LeafLength,1,3);
            
            % Store parent info.
            LeafParent(NLeaf) = JRand(iBlock);

        end

    end
    
    % Remove empty rows from the end, if neccessary.
    if NLeaf < NLeafEstimate
        LeafDimensions = LeafDimensions(1:NLeaf,:);
        LeafParent      = LeafParent(1:NLeaf,:);
    end

    % Find maximum leaf size.
    MaxLeafSize = max(max(LeafDimensions));

end

% Default function for leaf orientation distribution.
function [LeafDir, LeafNormal] = default_fun_angle(~, ...
                                                   ~,...
                                                   TwigDir,...
                                                   varargin)
    %-
    
    % Params: relative_height, relative_position, branch_order, radius
    
    % Direction where leaves try to turn, i.e., reference normal.
    RefDir = [0 0 1];
    
    % Angle limit in degrees. If the initial normal differs from the
    % reference direction less than this value, the reference direction is
    % used as the normal. Otherwise further computations are made.
    AngLim = 20;
    
    % Number of leaves.
    NLeaf = size(TwigDir,1);
    
    % Initialize variables that will hold the individual directions and
    % normals.
    LeafDir    = zeros(NLeaf,3);
    LeafNormal = zeros(NLeaf,3);
    
    for iLeaf = 1:NLeaf
        
        % Initial direction of the leaf is the direction of the twig.
        Dir = TwigDir(iLeaf,:);
        
        % If the initial direction is the same as the reference direction,
        % the initial direction is updated by rotation.
        if all(RefDir == Dir)
           
            rotcor = 20/180*pi;
            R = rotation_matrix([1 0 0], rotcor);
            
            Dir = (R*Dir')';
            
        end
        
        % Find a unit vector perpendicular to the reference and initial
        % direction.
        Side = cross(Dir, RefDir);
        Side = Side./norm(Side);
        
        % Compute initial normal penpendicular to the direction.
        Normal = cross(Side,Dir);
        
        % Compute angle between normal and reference direction.
        Alpha = atan2d(norm(cross(Normal,RefDir)),dot(Normal,RefDir));
        
        % If angle is less than given limit.
        if Alpha < AngLim
            
            % Set normal as reference direction.
            LeafNormal(iLeaf,:) = RefDir;
            
            % Find leaf direction perpendicular to the normal.
            Dir = cross(RefDir,Side);
            Dir = Dir./norm(Dir);
            LeafDir(iLeaf,:) = Dir;
            
        % If initial normal is too far from reference direction.
        else
            
            % Rotate initial normal and directions by angle limit towards
            % reference direction.
            R = rotation_matrix(Side,-pi*(Alpha - AngLim)/180);
            LeafNormal(iLeaf,:) = (R*Normal')';
            LeafDir(iLeaf,:) = (R*Dir')';
            
        end
        
    end

end
