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

classdef LeafModelTriangle < LeafModel
% Class for storing triangular leaf geometry. Class properties define a
% base geometry that is positioned, scaled and rotated to receive exact
% geometry of included leaves. The base can have any number of triangles.
% LeafModelTriangle is a subclass of LeafModel. 
    
    properties

        % Coordinates of the vertices forming a leaf.
        base_vertices = zeros(0,3);

        % Relative dimensions of the base. Base leaf length always equals
        % one, i.e., base_dimensions(2) == 1.
        base_dimensions = zeros(0,3);

        % Indices of vertices forming the triangular faces.
        base_triangles = zeros(0,3);
        
        % Area of the base.
        base_area = 0;

        % Cell array of indices forming the ngons of the base.
        base_ngons;

        % Number of triangles per leaf.
        triangle_count = 0;

    end
    
    properties(Access=private)
        
        % Computed vertices of included leaves.
        % Used mainly for optimization.
        leaf_triangle_vertices;

        % Computed normals of included leaf triangles.
        leaf_triangle_normals;
        
        % Projected values of the included leaf triangles of their
        % respective normals.
        leaf_triangle_zvalue;

    end


    methods

        function ob = LeafModelTriangle(vertices, tris, varargin)
        % Constructor of class object.

            % Store base geometry.
            ob.base_vertices = vertices;
            ob.base_triangles = tris;
            
            % Get triangle count.
            ob.triangle_count = size(tris,1);
            
            % Compute base dimensions from extreme points.
            mi = min(vertices,[],1);
            ma = max(vertices,[],1);
            ob.base_dimensions = ma - mi;
            
            % Compute base area as the sum of all triangles.
            v1 = bsxfun(@minus, vertices(tris(:,2),:),...
                                vertices(tris(:,1),:));
            %-
            v2 = bsxfun(@minus, vertices(tris(:,3),:),...
                                vertices(tris(:,1),:));
            %-
            
            A = cross(v1,v2,2);
            ob.base_area = sum(sqrt(sum(A.^2,2)))/2;

            % Leaf includes ngons.
            % Cell-list of vertex indices defining them.
            if nargin > 2
                if iscell(varargin{1}) && numel(varargin{1})
                    ob.base_ngons = varargin{1};
                end
            end
            
            % Initialize matrices for optimization.
            if nargin > 3
                assert(isnumeric(varargin{2}) && isscalar(varargin{2}),...
                 'Initialization parameter should be a positive integer.');
                %-
                N = varargin{2};
                
                ob.leaf_start_point = zeros(N,3);
                ob.leaf_direction   = zeros(N,3);
                ob.leaf_normal      = zeros(N,3);
                ob.leaf_scale       = zeros(N,3);
                
                ob.leaf_parent      = zeros(N,1);
                ob.twig_start_point = zeros(N,3);
                
                ob.leaf_triangle_vertices = nan(N,9,ob.triangle_count);
                ob.leaf_triangle_normals = nan(N,3,ob.triangle_count);
                ob.leaf_triangle_zvalue = nan(N,ob.triangle_count);
            end

        end

        function bounding_box(ob)

            MaxLeafSize = max(ob.leaf_scale(:));

            % Extreme points.
            Sp_pr = bsxfun(@plus,ob.leaf_start_point, MaxLeafSize);
            Sp_mr = bsxfun(@plus,ob.leaf_start_point,-MaxLeafSize);
            
            % Upper and lower limits for the tree in the z-axis.
            ob.tree_limits = [min([Sp_pr; Sp_mr],[],1); ...
                              max([Sp_pr; Sp_mr],[],1)];
            %-

        end
        
        function trim_slack(ob)
        % As the main matrices can be initialized with empty rows, this
        % function can be used to trim the excess rows after all leaves
        % have been added.
        
            ob.leaf_start_point = ob.leaf_start_point(1:ob.leaf_count,:);
            ob.leaf_direction   = ob.leaf_direction(1:ob.leaf_count,:);
            ob.leaf_normal      = ob.leaf_normal(1:ob.leaf_count,:);
            ob.leaf_scale       = ob.leaf_scale(1:ob.leaf_count,:);
            
            ob.leaf_parent      = ob.leaf_parent(1:ob.leaf_count,:);
            ob.twig_start_point = ob.twig_start_point(1:ob.leaf_count,:);

            ob.leaf_triangle_vertices = ...
                            ob.leaf_triangle_vertices(1:ob.leaf_count,:,:);
            %-
            ob.leaf_triangle_normals = ...
                            ob.leaf_triangle_normals(1:ob.leaf_count,:,:);
            %-
            ob.leaf_triangle_zvalue = ...
                            ob.leaf_triangle_zvalue(1:ob.leaf_count,:);
            %-
        end

        function [index,area] = add_leaf(ob, origin, dir, normal, scale,...
                                         parent, twig, tris)
        % Add a leaf with input parameters to the leaf collection.

            % Object index for new leaf.
            index = ob.leaf_count + 1;
            
            % Increase included leaf count.
            ob.leaf_count = index;
            
            % If triangles have not been given, compute them.
            if nargin < 8
                tris = ob.triangles(origin,dir,normal,scale);
            end

            % Store leaf triangle vertices for possible later computations.
            ob.leaf_triangle_vertices(index,:,:) = tris';

            % Two vectors on triangle: v2 - v1 and v3 - v1.
            Vec1 = tris(:,4:6) - tris(:,1:3);
            Vec2 = tris(:,7:9) - tris(:,1:3);

            % Compute triangle normals.
            TriNormals = cross(Vec1,Vec2,2);

            % Normalize normal the normal way.
            TriNormals = bsxfun(@times,...
                                TriNormals,...
                                1./sqrt(sum(TriNormals.^2,2)));
            %-
            
            % Store normal(s).
            ob.leaf_triangle_normals(index,:,:) = TriNormals';

            % Compute and store z-value.
            ob.leaf_triangle_zvalue(index,:) = dot(TriNormals,...
                                                   tris(:,1:3),...
                                                   2);
            %-

            ob.leaf_start_point(index,:) = origin;
            ob.leaf_direction(index,:) = dir;
            ob.leaf_normal(index,:)  = normal;
            ob.leaf_scale(index,:) = scale;
            
            ob.leaf_parent(index,:) = parent;
            ob.twig_start_point(index,:) = twig;
            
            % Leaf area of new leaf is received by
            % scaling the area of the base by the 
            % Y-scale of the leaf.
            area = ob.base_area*(scale(2)^2);
            
            % Increase leaf area. 
            ob.leaf_area = ob.leaf_area + area;

        end

        
        function hit = leaf_intersect(ob, Candidates, NewTris)
        % Check if leaf with given triangles intersects any of the included
        % leaves.
        
            % Number of candidate leaves.
            NCandidates = length(Candidates);

            % Number of triangles in a leaf. 
            NTri = ob.triangle_count;

            % Hit detected.
            hit = false;

            CandatesTris    = ob.leaf_triangle_vertices(Candidates,:,:);
            CandatesNormals = ob.leaf_triangle_normals(Candidates,:,:);
            CandatesZvalues = ob.leaf_triangle_zvalue(Candidates,:);

            % Normals of triangles of new leaves.
            NewNormal  = zeros(NTri,3);
            NewZvalues = nan(NTri,1);

            % Iterate over candidates.
            for iCandidate = 1:NCandidates
                
                % Values for the triangle(s) of current candidate.
                CandNormals = CandatesNormals(iCandidate,:,:);
                CandTris    = CandatesTris(iCandidate,:,:);
                CandLimits  = CandatesZvalues(iCandidate,:);
                
                if ob.triangle_count > 1
                    CandNormals = squeeze(CandNormals)';
                    CandTris    = squeeze(CandTris)';
                end

                % Iterate over triangles in candidate.
                for iCandTri = 1:NTri
                    
                    CandTriNormal = CandNormals(iCandTri,:);
                   
                    % Iterate over input triangles.
                    for iNewTri = 1:NTri
                       
                        TriZvalues = CandTriNormal ...
                                   * reshape(NewTris(iNewTri,:),3,3)...
                                   - CandLimits(iCandTri);
                        %-

                        if all(TriZvalues < 0) || all(TriZvalues >= 0)
                            continue;
                        end

                        if isnan( NewZvalues(iNewTri) )
                            TriNormal = cross(NewTris(iNewTri,7:9) ...
                                              - NewTris(iNewTri,1:3),...
                                              NewTris(iNewTri,4:6) ...
                                              -NewTris(iNewTri,1:3));
                            %-
                            
                            TriNormal = TriNormal./norm(TriNormal);

                            % Project one vertex to new triangle normal.
                            ZLimit = TriNormal*NewTris(iNewTri,1:3)';

                            % Store normal and z-limit.
                            NewNormal(iNewTri,:) = TriNormal;
                            NewZvalues(iNewTri) = ZLimit;
                        else
                            % Get triangle normal and z-limit.
                            TriNormal = NewNormal(iNewTri,:);
                            ZLimit = NewZvalues(iNewTri);
                        end

                        CandZvalues = TriNormal...
                                    * reshape(CandTris(iCandTri,:),3,3)...
                                    - ZLimit;
                        %-

                        if all(CandZvalues < 0) || all(CandZvalues >= 0)
                            continue;
                        end

                        % Iterate over edges in new triangle.
                        for iEdge = 1:3

                            % Form origin and direction.
                            if iEdge == 1
                                RayOrigin = NewTris(iNewTri,1:3);
                                RayEnd    = NewTris(iNewTri,4:6);
                            elseif iEdge == 2
                                RayOrigin = NewTris(iNewTri,4:6);
                                RayEnd    = NewTris(iNewTri,7:9);
                            elseif iEdge == 3
                                RayOrigin = NewTris(iNewTri,7:9);
                                RayEnd    = NewTris(iNewTri,1:3);
                            end

                            RayDir = RayEnd - RayOrigin;
                            
                            % Normalize direction vector.
                            RayDir = RayDir./norm(RayDir);

                            % Compute limiting distance.
                            DistanceLimit = RayDir*(RayEnd-RayOrigin)';

                            % Check line segment - triangle intersection.
                            hit = LeafModelTriangle.LineSegTriIntersect(...
                                                RayOrigin,...
                                                RayDir,...
                                                DistanceLimit,...
                                                CandTris(iCandTri,1:3),...
                                                CandTris(iCandTri,4:6),...
                                                CandTris(iCandTri,7:9));
                            %-
                            
                            % Return if intersection happens.
                            if hit
                                return;
                            end

                        end

                        % Iterate over edges in candidate triangle.
                        for iEdge = 1:3

                            % Form origin and direction.
                            if iEdge == 1
                                RayOrigin = CandTris(iCandTri,1:3);
                                RayEnd    = CandTris(iCandTri,4:6);
                            elseif iEdge == 2
                                RayOrigin = CandTris(iCandTri,4:6);
                                RayEnd    = CandTris(iCandTri,7:9);
                            elseif iEdge == 3
                                RayOrigin = CandTris(iCandTri,7:9);
                                RayEnd    = CandTris(iCandTri,1:3);
                            end
                            
                            RayDir = RayEnd - RayOrigin;
                            
                            % Normalize direction vector.
                            RayDir = RayDir./norm(RayDir);

                            % Compute limiting distance.
                            DistanceLimit = RayDir*(RayEnd-RayOrigin)';

                            % Check line segment - triangle intersection.
                            hit = LeafModelTriangle.LineSegTriIntersect(...
                                                RayOrigin,...
                                                RayDir,...
                                                DistanceLimit,...
                                                NewTris(iNewTri,1:3),...
                                                NewTris(iNewTri,4:6),...
                                                NewTris(iNewTri,7:9));
                            %-
                            
                            % Return if intersection happens.
                            if hit
                                return;
                            end

                        end

                    end
                    
                end

            end

        end

        
        function varargout = triangles(ob, origin, dir, normal, scale)
        % Function for computing the geometry of a leaf based on the 
        % transformation parameter inputs. 

            [Vertices, Faces] = compute_geometry(ob, false, origin, ...
                                                 dir, normal, scale);
            %-

            if nargout == 2
                varargout{1} = Vertices;
                varargout{2} = Faces;
            else

                varargout{1} = cat(2,Vertices(Faces(:,1),:),...
                                     Vertices(Faces(:,2),:),...
                                     Vertices(Faces(:,3),:));
                %-
            end

        end

        
        function h = plot_leaves(ob, varargin)
        % Plot accepted leaves using the PATCH function.
            
            [Vertices, Faces] = compute_geometry(ob, false);

            tris = cat(2,Vertices(Faces(:,1),:),...
                         Vertices(Faces(:,2),:),...
                         Vertices(Faces(:,3),:));
            
            % Reshape the triangle data for plotting.
            X = tris(:,[1 4 7])';
            Y = tris(:,[2 5 8])';
            Z = tris(:,[3 6 9])';
            
            % Plot triangles.
            h = patch(X,Y,Z,1,varargin{:});

        end

        function [Vertices, Faces] = compute_geometry(ob, fNgon, varargin)
        % Function to compute exact geometry of leaves by transforming the
        % basis geometry.
        %
        % Possible formats:
        % [Vert, Faces] = ob.triangles(fNgon)
        % [Vert, Faces] = ob.triangles(fNgon, Filter)
        % [Vert, Faces] = ob.triangles(fNgon, origin, dir, normal, scale)

            % If leaf transformation parameters are given as input,
            % use those single parameters.
            if nargin == 6
                
                origin = varargin{1};
                dir    = varargin{2};
                normal = varargin{3};
                scale  = varargin{4};

                % Base vertices that are modified to get
                % vertices of given leaf.
                vert  = ob.base_vertices;

                % Scaling.
                Vertices = bsxfun(@times,vert,scale);

                % Coordinate system.
                E = [cross(normal,dir); dir; normal];

                % Rotation.
                Vertices = Vertices*E;

                % Translation.
                Vertices = bsxfun(@plus,Vertices,origin);

                % Indices of faces in leaf.
                Faces = ob.base_triangles;

            % If transformation parameters are not given, 
            % use all the leaves in the model.
            else
            
                % Number of leaves in model.
                NLeaf = ob.leaf_count;

                % Logical vector for filtering only some of the leaves for
                % export.
                if nargin < 3
                    Filter = true(NLeaf,1);
                else
                    Filter = varargin{1};
                end

                if isempty(ob.base_ngons)
                    fNgon = false;
                end

                % Number of vertices in base.
                NBaseVert = size(ob.base_vertices,1);

                % Matrix of vertices, row == vertex.
                Vertices = zeros(nnz(Filter)*NBaseVert,3);

                % Number of exported leaves.
                jLeaf = 0;

                % Iterate over leaves in model.
                for iLeaf = 1:NLeaf

                    % Check if filtered in.
                    if Filter(iLeaf)
                        % Increase number of exported leaves.
                        jLeaf = jLeaf + 1;
                    else
                        % Ignore filtered leaf.
                        continue;
                    end

                    % Get single leaf parameters and convert to vertices.
                    origin = ob.leaf_start_point(iLeaf,:);
                    scale  = ob.leaf_scale(iLeaf,:);
                    dir    = ob.leaf_direction(iLeaf,:);
                    normal = ob.leaf_normal(iLeaf,:);

                    % Compute vertices by transforming leaf base.
                    vert  = ob.base_vertices;

                    % Scaling.
                    vert = bsxfun(@times,vert,scale);

                    % Coordinate system.
                    E = [cross(normal,dir); dir; normal];

                    % Rotation.
                    vert = vert*E;

                    % Translation.
                    vert = bsxfun(@plus,vert,origin);

                    % Store resulting vertices.
                    Vertices((jLeaf-1)*NBaseVert+1:jLeaf*NBaseVert,:) ...
                        = vert;
                    %-
                end

                % Number of included leaves.
                NLeaf = jLeaf;

                % Using ngons
                if fNgon

                    % Flag: equal-sized faces.
                    fEqualNgons = true;

                    % Vertical catenation only works if faces have equal number
                    % of vertices.
                    try
                        BaseNgons = vertcat(ob.base_ngons{:});
                    catch
                        fEqualNgons = false;
                    end

                    % All ngons have equal number of vertices.
                    if fEqualNgons

                        % Number of triangles in base.
                        NNgon = size(BaseNgons,1);

                        % Indices of base triangle face vertices.
                        ngons = repmat(BaseNgons,NLeaf,1);

                        add = repmat(0:1:NLeaf-1,NNgon,1);
                        add = add(:)*NBaseVert;

                        Faces = bsxfun(@plus,ngons,add);

                    % Ngon vertex count varies.
                    else

                        % Base ngons as cell array.
                        BaseNgons = ob.base_ngons;

                        % Number of ngons in base.
                        NNgon = numel(BaseNgons);

                        % Cell to store faces.
                        Faces = cell(NLeaf*NNgon,1);

                        for iLeaf = 1:NLeaf

                            for iNgon = 1:NNgon

                                % Offset index by previous leaf count,
                                % and store to cell array.
                                Faces{(iLeaf-1)*NNgon+iNgon} = BaseNgons{iNgon} ...
                                                             + (iLeaf-1)*NNgon;
                                %-

                            end

                        end

                    end

                else

                    % Number of triangles in base.
                    NTri = size(ob.base_triangles,1);

                    % Indices of base triangle face vertices.
                    tris = repmat(ob.base_triangles,NLeaf,1);

                    % Offset vector to account for previous leaves.
                    add = repmat(0:1:NLeaf-1,NTri,1);
                    add = add(:)*NBaseVert;

                    % Store offset face indices to matrix.
                    Faces = bsxfun(@plus,tris,add);

                end
            end

            

        end

        function export_geometry(ob, format, fNgon, file, d, OriginOffset, Filter, varargin)
        % Compute accepted leaf geometry and export to a file. Currently supports
        % Wavefront OBJ-format.

            % Flag if extra parameters are given to print to the file.
            fExtraCols = false;

            if nargin > 7
                fExtraCols = true;

                % Number of extra columns.
                NCol = numel(varargin);
            end

            % Set default precision.        
            if nargin < 5
                d = 3;
            end
            
            if nargin < 6
                fOriginOverride = false;
                OriginOffset = [];
            end
            
            if nargin < 7 || isempty(Filter)
                Filter = true(ob.leaf_count,1);
            end

            switch lower(format)

                % Export leaves in Wavefront OBJ-format.
                case 'obj'

                    % Convert leaf transformation parameters into vertices and
                    % faces.
                    [Vertices, Faces] = ob.compute_geometry(fNgon, Filter);

                    % Override origin if given.
                    if fOriginOverride && not(isempty(OriginOffset))
                        Vertices = bsxfun(@minus,Vertices,OriginOffset);
                    end
                    
                    % Write resulting vertices and faces to file.
                    LeafModelTriangle.export_vert_face_obj(Vertices,Faces,d,file);

                case {'ext_obj', 'extobj'}

                    % Base vertices.
                    BaseVertices = ob.base_vertices;

                    % Base vertices either as polygons or triangles.
                    if fNgon
                        BaseFaces = ob.base_ngons;
                    else
                        BaseFaces = ob.base_triangles;
                    end

                    % Open file stream.
                    fid = fopen(file,'w');

                    % Write base geometry to file in OBJ syntax.
                    LeafModelTriangle.export_vert_face_obj(BaseVertices,BaseFaces,d,fid);

                    % Initialize format strings:
                    % Single element.
                    ft = ['%' num2str(d+2) '.' num2str(d) 'g '];
                    % Three element vector.
                    fmt = strtrim(repmat(ft,1,3));
                    ft = strtrim(ft);

                    % Iterate over leaves.
                    for iLeaf = 1:ob.leaf_count

                        % Skip filtered leaves.
                        if ~Filter(iLeaf)
                            continue;
                        end

                        % Print line type.
                        fprintf(fid,'L ');
                        % Print leaf parameters.
                        fprintf(fid,[fmt ' '],ob.twig_start_point(iLeaf,:));
                        fprintf(fid,[fmt ' '],ob.leaf_start_point(iLeaf,:));
                        fprintf(fid,[fmt ' '],ob.leaf_direction(iLeaf,:));
                        fprintf(fid,[fmt ' '],ob.leaf_normal(iLeaf,:));
                        fprintf(fid,fmt,ob.leaf_scale(iLeaf,:));

                        % Print extra columns if present.
                        if fExtraCols
                            % Add extra space to separate extra columns.
                            fprintf(fid, ' ');

                            % Print extra columns.
                            for iCol = 1:NCol
                                fprintf(fid,ft,varargin{iCol}(iLeaf,:));
                            end
                        end

                        % End line.
                        fprintf(fid, '\n');

                    end

                    % Close file stream.
                    fclose(fid);

                    % Print number of exported leaf configurations.
                    disp(['Exported '      num2str(nnz(Filter)) ...
                          ' leaf transformation parameters.']);

                % Unknown export format identifier.
                otherwise
                    fprintf('Unknown export format: "%s"\n', format);
            end
        end

        function export_obj(ob, file, d, OriginOffset, Filter)
        % Legacy alias for exporting geometry in Wavefront OBJ format.
        % Instead use EXPORT_GEOMETRY method when possible.

            ob.export_geometry(ob, 'OBJ', true, file, d, OriginOffset, Filter);
        end

    end

    methods(Access=protected)

    end
    
    methods(Static)
        
        function fIntersect = LineSegTriIntersect(Origin, Dir, ...
                                                          DistLim, ...
                                                          p1, p2, p3)
        %-
        % Function that checks whether a line segment intersects with a
        % triangle defined by three vertices.

            % Intersection flag: true if hit between line segments points.
            fIntersect = false;

            % Vectors on triangle plane.
            vec1 = p2 - p1;
            vec2 = p3 - p1;
            
            % Cross product of line segment direction and one vector on
            % plane.
            cross1 = cross(Dir,vec2);
            
            % Check if vectors were parallel.
            proj1  = dot(vec1,cross1);

            % Return false.
            if abs(proj1) < eps
                return;
            end
            
            % Vector from triangle corner to segment origin.
            op1 = Origin - p1;
            cross2 = cross(op1,vec1);
            
            % Barycentric coordinates of a line parallel to the line
            % segment on the plane defined by the two triangle vectors.
            u = dot(op1,cross1) / proj1;
            v = dot(Dir,cross2) / proj1;
            
            % Check if hit point is inside triangle, return false
            % otherwise.
            if u < 0 || v < 0 || u + v > 1
                return;
            end
            
            % Compute hit distance from line segment origin.
            Dist = dot(vec2,cross2) / proj1;
            
            % Hit occurred if hit distance between limits. Edge hits are
            % discarded as equal operators are not included.
            fIntersect = Dist > 0 && Dist < DistLim;
            
        end
        
        function export_vert_face_obj(Vertices,Faces,d,file)
        % Export vertices and faces in Wavefront OBJ-format.
            
            % Set precision formatter.
            if length(d) > 1
                ft = ['%' num2str(d(2)) '.' num2str(d(1)) 'g'];
            else
                ft = ['%.' num2str(d(1)) 'g'];
            end

            % Flag to close file and the end.
            closefile = false;

            % Check if file name and not file stream.
            if ischar(file)

                % Open file stream with file name.
                fid = fopen(file,'w');
                % Set file to close at the end.
                closefile = true;

            else
                % Otherwise a file stream is given as input.
                fid = file;
            end

            % Number of vertices.
            NVertex = size(Vertices,1);

            % Print format.
            ft = ['v ' ft ' ' ft ' ' ft '\n'];

            % Print vertices.
            for iVertex = 1:NVertex

                % Print vertex to file.
                fprintf(fid,ft,Vertices(iVertex,:));

            end

            % Number of faces.
            NFace = size(Faces,1);

            % Faces can have constant number of vertices or the number can
            % vary. Constant count implies a matrix, varying count needs a
            % cell array.
            fCellFaces = false;

            if iscell(Faces)
                fCellFaces = true;
            end

            % Number of vertices in face.
            NFaceVertex = size(Faces,2);

            % Face vertex print format.
            ft = repmat('%d ', 1, NFaceVertex);
            % Remove trailing space.
            ft = ft(1:end);

            % Print faces.
            for iFace = 1:NFace

                if fCellFaces

                    % Vertices of face.
                    Face = Faces{iFace};

                    % Vertices in face
                    NFaceVertex = length(Face);

                    % Number of vertices affects print format.
                    ft = repmat('%d ', 1, NFaceVertex);
                    ft = ft(1:end);

                else
                    % Vertices of face.
                    Face = Faces(iFace,:);
                end

                % Print face vertex indices to file.
                fprintf(fid,['f ' ft '\n'],Face);

            end

            % If the file stream was opened in this function, 
            % close the stream.
            if closefile
                fclose(fid);
            end

            % Display the number of vertices and faces printed to the file.
            disp(['Exported '      num2str(NVertex) ...
                  ' vertices and ' num2str(NFace) ' faces.']);
            %-

        end
        
    end

end