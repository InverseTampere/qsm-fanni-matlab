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

classdef QSMBCylindrical < QSMB
% Container for cylindrical quantitative structure model block information.
% A subclass of the abstract class QSMB.

    properties

        % Cylinder-level properties.
        cylinder_start_point = [];
        cylinder_axis        = [];
        cylinder_length      = [];
        cylinder_radius      = [];
        cylinder_parent      = [];
        cylinder_extension   = [];

        cylinder_branch_index    = [];
        cylinder_branch_order    = [];
        cylinder_index_in_branch = [];

        cylinder_end_point = [];
        cylinder_mid_point = [];
        cylinder_is_last   = [];

        % Number of branches.
        branch_count = 0;

        % Branch-level properties.
        branch_order  = [];
        branch_parent = [];
        branch_volume = [];
        branch_length = [];
        branch_angle  = [];
        branch_height = [];
    end
    
    properties(Access=private)
        
        % Coordinate system change matrix for each block.
        % Used for optimization in intersection checking.
        cylinder_coordinate_system;
        
    end
    
    properties(Constant)
        
        % Debug flag for easier development.
        debug = false;
    end


    methods

        function ob = QSMBCylindrical(ModelData)
        % Constructor of class object.

            % Set number of blocks.
            ob.block_count = size(ModelData{1},1);

            % Set cylinder-level properties.
            ob.cylinder_start_point = ModelData{1}(:,3:5);
            ob.cylinder_axis        = ModelData{1}(:,6:8);
            ob.cylinder_length      = ModelData{1}(:,2);
            ob.cylinder_radius      = ModelData{1}(:,1);
            ob.cylinder_parent      = ModelData{1}(:,9);
            ob.cylinder_extension   = ModelData{1}(:,10);

            ob.cylinder_end_point = ob.cylinder_start_point ...
                                  + bsxfun(@times,ob.cylinder_axis,...
                                                  ob.cylinder_length);
            %-
            
            ob.cylinder_mid_point = ob.cylinder_start_point ...
                                  + bsxfun(@times,ob.cylinder_axis,...
                                                  ob.cylinder_length/2);
            %-
            
            % Find extreme points, add and substract radius to each
            % cylinder start and end point to be sure.
            Ep_pr = bsxfun(@plus,ob.cylinder_end_point,...
                                 ob.cylinder_radius);
            %-
            Ep_mr = bsxfun(@plus,ob.cylinder_end_point,...
                                -ob.cylinder_radius);
            %-
            Sp_pr = bsxfun(@plus,ob.cylinder_start_point,...
                                 ob.cylinder_radius);
            %-
            Sp_mr = bsxfun(@plus,ob.cylinder_start_point,...
                                -ob.cylinder_radius);
            %-
            
            % Upper and lower limits for the tree in the z-axis.
            ob.tree_limits = [min([Ep_pr; Ep_mr; Sp_pr; Sp_mr],[],1); ...
                              max([Ep_pr; Ep_mr; Sp_pr; Sp_mr],[],1)];
            %-

            ob.cylinder_branch_index    = ModelData{1}(:,11);
            ob.cylinder_branch_order    = ModelData{1}(:,12);
            ob.cylinder_index_in_branch = ModelData{1}(:,13);

            % Set number of branches.
            ob.branch_count = size(ModelData{2},1);

            % Set branch-level properties.
            ob.branch_order  = ModelData{2}(:,1);
            ob.branch_parent = ModelData{2}(:,2);
            ob.branch_volume = ModelData{2}(:,3);
            ob.branch_length = ModelData{2}(:,4);
            ob.branch_angle  = ModelData{2}(:,5);

            % Set height of the branches if included in the input.
            fManualHeight = true;
            if size(ModelData{2},2) > 5
                ob.branch_height = ModelData{2}(:,6);
                fManualHeight = false;
            else
                % Otherwise initialize as zeros and compute later.
                ob.branch_height = zeros(ob.branch_count,1);
            end


            % Find the last cylinder in each branch.

            % Temp. storage.
            IsLast = false(ob.block_count,1);

            % Iterate over all branch numbers.
            for iBranch = 1:ob.branch_count
                
                % Indices of cylinder in given branch number.
                JBranch = find(ob.cylinder_branch_index == iBranch);

                % If no cylinders in branch, skip.
                if isempty(JBranch)
                    continue;
                end

                % If branch height was not given as input, compute from
                % cylinders as a mean.
                if fManualHeight
                    meanp = bsxfun(@times,ob.cylinder_axis(JBranch,:),...
                                          ob.cylinder_length(JBranch)/2)...
                          + ob.cylinder_start_point(JBranch,:);
                    %-
                    ob.branch_height(iBranch) = mean(meanp(:,3));
                end
                
                % Sort, from start to finnish.
                [~,jLast] = max(ob.cylinder_index_in_branch(JBranch));

                % Set flag for last cylinder.
                IsLast(JBranch(jLast)) = true;

            end

            % Store temporary information.
            ob.cylinder_is_last = IsLast;
            
            % Set default values for twig distribution.
            ob.fun_twig_distribution = @default_twig_param_dist;
            ob.twig_length_limits = [0.02, 0.05];
            
            % Initialize cylinder coordinate system variable. Used during
            % triangle intersection detection. If matrix contains NaNs, the
            % respective matrix is computed during execution.
            ob.cylinder_coordinate_system = nan(3,3,ob.block_count);

        end
       
        % Function for plotting cylinders in model. Implementation is
        % external.
        h = plot_cylinders(ob,varargin)

        function Props = get_block_properties(ob,PropNames,vargin)
        % Function to read (or compute, in some cases) properties of
        % cylindrical blocks. Property names are given in a cell array of
        % strings. It is also possible to restrict the property
        % computations to certain cylinders with the variable <indices>.
        % The property values are returned as fields of a struct.
        %
        % Examples:
        %
        % props = ob.get_block_properties(properties)
        % props = ob.get_block_properties(properties, indices)

            % Check if filtering by cylinder ID is required.
            if nargin > 2
                JCyl = vargin{1};
            else
                JCyl = 1:ob.block_count;
            end

            % Number of properties to compute.
            NProp = numel(PropNames);

            % Properties to return.
            Props = [];

            for iProp = 1:NProp

                % Current property label.
                Label = lower(PropNames{iProp});

                switch Label

                    % Basic properties that are already stored in the
                    % object. Return appropriate matrices.
                    case {'start_point', 'axis', 'length',...
                          'radius', 'parent', 'extension',...
                          'branch_index','branch_order',...
                          'index_in_branch','is_last'}
                        %-

                        Props.(Label) = ob.(['cylinder_' Label])(JCyl,:);

                    % Relative height of cylinder, normalized by tree
                    % extreme points.
                    case 'relative_height'
                        
                        CylMeanHeight = ...
                                  mean([ob.cylinder_start_point(JCyl,3),...
                                        ob.cylinder_end_point(JCyl,3)],2);
                        %-
            
                        Props.(Label) = ...
                               (CylMeanHeight - ob.tree_limits(1,3)) / ...
                               (ob.tree_limits(2,3) - ob.tree_limits(1,3));
                        %-

                    % Relative position of cylinder along the branch.
                    % Computed from the end point of the cylinder. Thus,
                    % the relative position of the last cylinder in a
                    % branch is always one.
                    case 'relative_position'

                        BranchIndex   = ob.cylinder_branch_index(JCyl);
                        IndexInBranch = ob.cylinder_index_in_branch(JCyl);
                        Length        = ob.cylinder_length(JCyl);


                        RelativePosition = zeros(length(JCyl),1);

                        % Compute cylinder position on branch.
                        % Iterate over all branch numbers.
                        for i = 1:max(BranchIndex)
                            
                            % Indices of cylinder in given branch number.
                            IBranch = (BranchIndex == i);

                            if not(any(IBranch))
                                continue;
                            end
                            
                            % Sort from start to finnish.
                            [~,JCylinder] = sort(IndexInBranch(IBranch));
                            
                            JBranch = find(IBranch);
                            
                            % Cumulative position of end points.
                            CumulativeLength = ...
                                        cumsum(Length(JBranch(JCylinder)));

                            % Total length.
                            TotalBranchLength = CumulativeLength(end);

                            % Relative position of cylinders.
                            RelativePosition(JBranch) = CumulativeLength...
                                                      / TotalBranchLength;
                        end

                        Props.(Label) = RelativePosition;

                end

            end

        end
        
        
        function [TwigStart, TwigEnd] = generate_twigs(ob, LeafParent)
        % Generate twigs connecting leaves to the cylinders. The input
        % <LeafParent> is a sorted vector of cylinder indices of parents 
        % of the leaves. A single cylinder can occur multiple times. The 
        % function returns the start and end points of the twigs.
            
            % Number of twigs to genereate.
            NTwig = length(LeafParent);
            
            % Parameters of cylinders with twigs. Same parent can be
            % repeated multiple times.
            sp      = ob.cylinder_start_point(LeafParent,:);
            ax      = ob.cylinder_axis(LeafParent,:);
            h       = ob.cylinder_length(LeafParent);
            r       = ob.cylinder_radius(LeafParent);
            is_last = ob.cylinder_is_last(LeafParent);
            
            % Parameters of the twigs. Row = single twig.
            % Columns:
            %   1: position on axis on the inverval [0,1].
            %   2: position on radial axis on the inverval [0,1]. 
            %      Equals 1 if on side, less than 1 when on head.
            %   3: rotation around axis. 0 upwards, Pi downwards.
            %   4: twig elevation. -pi backwards, Pi towards tip.
            %   5: twig azimuth. -Pi/2 left, Pi/2 right.
            %   6: twig length.
            TwigParam = ob.fun_twig_distribution(sp, ax, h, r, is_last,...
                                                 ob.twig_length_limits);
            %-
            
            % Index of the last parent cylinder, used to check if cylinder
            % coordinate system needs to be updated.
            LastParent = [];
            
            % Initialize return values.
            TwigStart = zeros(NTwig,3);
            TwigEnd   = zeros(NTwig,3);
            
            % Iterate over twigs.
            for iTwig = 1:NTwig
                
                % Always compute on first run, and when parent index
                % changes.
                if iTwig == 1 || LeafParent(iTwig) ~= LastParent
                    
                    % Update parent index.
                    LastParent = LeafParent(iTwig);
                
                    % Axis pointing to the side, defining radial direction.
                    if all(ax(iTwig,:) == [0 0 1])
                        TwigSideBase = [0 1 0];
                        TwigUpBase = [1 0 0];
                    else
                        TwigSideBase = cross([0 0 1],ax(iTwig,:));
                        TwigSideBase = TwigSideBase/norm(TwigSideBase);

                        TwigUpBase = cross(ax(iTwig,:),TwigSideBase);
                    end
                    
                end
                
                % Debug: plot cylinder coordinate system axis vectors.
                if ob.debug
                    figure(2);
                    conf = [1 4];
                    TwigFront = ax(iTwig,:);
                    TwigUp = TwigUpBase;
                    TwigSide = TwigSideBase;
                    plot_axes(conf,1,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Rotation matrix to rotate around cylinder axis.
                Rax = rotation_matrix(ax(iTwig,:),TwigParam(iTwig,3));

                % Rotate required coordinate axes around cylinder axis.
                TwigUp   = (Rax*TwigUpBase')';
                TwigSide = (Rax*TwigSideBase')';

                % Compute twig start point.
                TwigStart(iTwig,:) = sp(iTwig,:) ...
                                   + ax(iTwig,:)*h(iTwig)...
                                                *TwigParam(iTwig,1) ...
                                   + TwigUp     *r(iTwig)...
                                                *TwigParam(iTwig,2);
                %-
                
                if ob.debug
                    plot_axes(conf,2,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Elevation parameter has to be inverted if the twig is
                % connected to the tip of the cylinder.
                if TwigParam(iTwig,2) == 1
                    ElParam = -TwigParam(iTwig,4);
                else
                    ElParam = TwigParam(iTwig,4);
                end
                
                % Rotation matrix to rotate twig either towards
                % cylinder axis.
                Rel = rotation_matrix(TwigSide,ElParam);

                % Rotate required coordinate axes around side axis.
                TwigUp    = (Rel*TwigUp')';
                TwigFront = (Rel*ax(iTwig,:)')';
                
                % Debug: plot updated coordinate system.
                if ob.debug
                    plot_axes(conf,3,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Debug: plot twig direction unit vector.
                if ob.debug
                    
                    colors = eye(3);
                    E = [TwigUp; TwigSide; TwigFront];
                    
                    figure(2);
                    for i = 1:3
                        plott([TwigStart(iTwig,:); ...
                               TwigStart(iTwig,:)+E(i,:)],...
                              '-','Color',colors(i,:));
                        %-
                        hold on;
                    end
                    hold off;
                    axis equal;
                    
                end

                % Twig is attached to the envelope.
                if TwigParam(iTwig,2) == 1
                                        
                    % Rotation matrix to rotate to wanted azimuth, i.e.,
                    % turn the twig axis off cylinder center.
                    Raz = rotation_matrix(TwigFront,TwigParam(iTwig,5));
                    
                    % Rotate required coordinate axes around 
                    TwigUp = (Raz*TwigUp')';
                    
                    % Debug: update for later plot.
                    if ob.debug
                        TwigSide = (Raz*TwigSide')';
                    end
                    
                    % Compute twig end point.
                    TwigEnd(iTwig,:) = TwigStart(iTwig,:) ...
                                     + TwigUp*TwigParam(iTwig,6);
                    %-
                    
                    
                % Twig is attached to the tip of the cylinder.
                else
                    
                    % Rotation matrix to rotate to wanted azimuth, i.e.,
                    % turn the twig axis off cylinder center.
                    Raz = rotation_matrix(TwigUp,-TwigParam(iTwig,5));
                    
                    % Rotate required coordinate axes around.
                    TwigFront = (Raz*TwigFront')';
                    
                    % Debug: update for later plot.
                    if ob.debug
                        TwigSide = (Raz*TwigSide')';
                    end
                    
                    % Compute twig end point.
                    TwigEnd(iTwig,:) = TwigStart(iTwig,:) ...
                                     + TwigFront*TwigParam(iTwig,6);
                    %-
                    
                end
                
                % Plot updated coordinate system.
                if ob.debug
                    plot_axes(conf,4,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
            end
            
        end


        function fIntersect = block_triangle_intersection(ob, JBlock, Tris)
        % Check if triangles intersect with any of the blocks with given
        % indices. Returns True if at least one of the given triangles
        % intersects any of the given blocks.

            % Number of cylinders.
            NCyl = length(JBlock);
            
            % Number of triangles.
            NTri = size(Tris,1);

            % Flag for intersection.
            fIntersect = false;

            % Iterate over each cylinder-triangle pair.
            for iCyl = 1:NCyl

                for iTri = 1:NTri

                    % Check intersection of each pair.
                    fIntersect = ob.cylinder_triangle_intersect(...
                                               JBlock(iCyl),...
                                               Tris(iTri,:),...
                                               false);
                    %-
                    
                    % If an intersection occurs, skip the rest and return.
                    if fIntersect
                        return;
                    end
                end


            end

        end
        
        % Function to check the intersection of a single cylinder and as
        % single triangle.
        intersect = cylinder_triangle_intersect(ob, jBlock, Tri,...
                                                edgehit, debug)
        %-


        %% toVoxels: Return center coordinates of voxels that contain parts
        % of the blocks.
        function BlockVoxelization = toVoxels(ob, edge, minp, maxp)

            % Initialize voxelized space object.
            BlockVoxelization = CubeVoxelization(edge, minp, maxp);
            
            % Iterate over blocks in the model to populate voxelization.
            for iCyl = 1:ob.block_count

                % Get cube coordinates of voxels occupied by cylinder.
                cc = ob.CylinderToVoxels(iCyl, edge, minp);

                % Add cylinder index to occupied voxels.
                BlockVoxelization.add_object_by_cc(cc,iCyl);

            end

        end
        
        function export_blender(ob, file, d, origin, varargin)
        % Print cylinder model parameters to file, e.g., for exporting 
        % to Blender, using the Blender QSM import addon.

            % Flag if extra parameters are given to print to the file.
            extracols = false;

            if nargin > 4
                extracols = true;
                % Number of extra columns.
                NCol = numel(varargin);
            end

            % Set origin override as zeros, if not given.
            if nargin < 4 || isempty(origin)
                origin = zeros(1,3);
            end

            % Set precision formatter.
            if length(d) > 1
                ft = [' %' num2str(d(2)) '.' num2str(d(1)) 'f'];
            else
                ft = [' %.' num2str(d(1)) 'f'];
            end

            % Number of cylinders.
            NCyl = ob.block_count;

            % Open file stream.
            fid = fopen(file,'w');

            for iCyl = 1:NCyl

                % Print to file:
                % - Branch ID
                % - Starting point
                % - Axis direction
                % - Length
                % - Radius
                fprintf(fid,['%d' repmat(ft,1,8)], ...
                        ob.cylinder_branch_index(iCyl), ...
                        ob.cylinder_start_point(iCyl,:) - origin, ...
                        ob.cylinder_axis(iCyl,:), ...
                        ob.cylinder_length(iCyl), ...
                        ob.cylinder_radius(iCyl));
                %-

                % Print extra columns if present.
                if extracols

                    for iCol = 1:NCol
                        fprintf(fid,ft,varargin{iCol}(iCyl,:));
                    end
                end

                % Print line end.
                fprintf(fid,'\n');
            end

            % Close stream.
            fclose(fid);

        end
        

    end

    methods(Access=protected)

        function cc = CylinderToVoxels(ob, iCyl, edge, minp)
        % Find voxel indices of all voxels a single cylinder occupies. The
        % cylinder is defined by an index, and the voxelization with an
        % edge length and a minimum point.

            % Get cylinder properties to create point samples.
            sp = ob.cylinder_start_point(iCyl,:);
            ax = ob.cylinder_axis(iCyl,:);
            h  = ob.cylinder_length(iCyl,:);
            r  = ob.cylinder_radius(iCyl,:);

            % Generate test points in and on cylinder.
            P = QSMBCylindrical.VoxelSamples(edge, sp, ax, h, r);

            % Compute cube coordinates of each sample point.
            cc = floor(bsxfun(@minus,P,minp)./edge) + 1;

            % Only include each voxel once.
            cc = unique(cc,'rows');

        end
        
    end
    
    methods(Static)

        function P = VoxelSamples(edge, sp, ax, h, r)
        % Generate points inside a cylinder with the given parameters:
        % start point <sp>, axis direction <ax>, length <h>, and radius
        % <r>. The <edge> parameter is used to chose the number of samples.

            % Average distance between points.
            PointDist = edge/2;

            % Number vertex rings inside cylinder.
            NRing = ceil(r/PointDist);

            % Radii of each ring. Note that zero is computed separately.
            RingRadius = linspace(0,1,NRing + 1);
            RingRadius = RingRadius(2:end);

            % Number of vertices on ring loop.
            NVertex = max(ceil( pi./asin(PointDist./(2*r*RingRadius)) ),3);

            % Number of vertical layers.
            NLayer = ceil(h/PointDist) + 1;

            % Heights at which ring layers are distributed at.
            z = linspace(0,h,NLayer);

            % Total number of points.
            NPoint = sum(NVertex)*NLayer;

            % Vertices.
            P = zeros(NPoint,3);

            iPoint = 1;

            for iRing = 1:NRing

                NPointRing = NVertex(iRing)*NLayer;

                ri = RingRadius(iRing);

                Ang = linspace(0,2*pi,NVertex(iRing) + 1);
                Ang = Ang(1:end-1);

                [x,y] = pol2cart(Ang,ri*r);

                X = repmat(x,NLayer,1);
                Y = repmat(y,NLayer,1);
                Z = repmat(z(:),1,NVertex(iRing));

                P(iPoint:iPoint+NPointRing-1,:) = [X(:), Y(:), Z(:)];

                iPoint = iPoint + NPointRing;

            end

            % Add points on axis.
            P = vertcat(P,[zeros(NLayer,2), z(:)]);

            % Coordinate system matrix.

            if ax(3) == 1
                R = eye(3);
            elseif ax(3) == -1
                R = [-1 0 0; 0 1 0; 0 0 -1];
            else
                axx = cross(ax,[0 0 1]);
                axx = axx./norm(axx);

                axy = cross(axx,ax);

                R = [axx; axy; ax];
            end


            % Rotate to cylinder coordinate system and translate.
            P = bsxfun(@plus,P*R,sp);
            
            % Debug: plot resulting points.
            if QSMBCylindrical.debug
                plott(P,'b.');
                hold on;
                axis equal;
            end

        end

    end

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

function TwigParam = default_twig_param_dist(sp, ax, h, r, is_last, len)
% Parameters of the twigs. Row = single twig.
% Columns:
%   1: position on axis on the inverval [0,1].
%      Equals 1 if on head.
%   2: position on radial axis on the inverval [0,1]. 
%      Equals 1 if on side, less than 1 when on head.
%   3: rotation around axis. From -Pi to Pi.
%   4: twig elevation. -Pi/2 backwards, Pi/2 towards tip.
%   5: twig azimuth. -Pi/2 left, Pi/2 right.
%   6: twig length.

    % Limits for even distributions.
    Limits = [0 1; 0 1; -pi pi; -pi/2 pi/2; -pi/2  pi/2; len];
    %Limits = [0 1; 0 1; -pi pi;     0 pi/2; -pi/4 -pi/4; len];

    NTwig = size(ax,1);
    NTwig = max(NTwig,size(sp,1));
    
    TwigParam = zeros(NTwig,6);
    
    for iParam = 1:6
        
        TwigParam(:,iParam) = Limits(iParam,1) ...
                            + (Limits(iParam,2) - Limits(iParam,1)) ...
                            * rand(NTwig,1);
        %-
        
    end
    
    % Ratio of area of circle and complete area.
    ratio = r./(r + 2*h);

    % Leafs that are connected to cylinder envelope and not the end 
    % circle.
    onhead = rand(NTwig,1) < ratio(:) & is_last(:);

    % Set axial translation for head-connected leaves to full.
    TwigParam(onhead,1) = 1;
    % Set radial translation for side-connected leaves to full.
    TwigParam(not(onhead),2) = 1;
    
    

end

function plot_axes(config, index, Origin, TwigUp,TwigSide,TwigFront)
% Function to plot coordinate system axis vectors. Used for debugging.

    colors = eye(3);
    E = [TwigUp; TwigSide; TwigFront];

    subplot(config(1),config(2),index);
    
    for i = 1:3
        plott([Origin; Origin+E(i,:)],'-','Color',colors(i,:));
        hold on;
    end
    
    hold off;
    axis equal;
    legend({'up','side','front'});

end