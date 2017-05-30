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

classdef CubeVoxelization < handle
% CubeVoxelization objects contains the voxelization parameters,
% i.e., minimum corner, edge length, and object distribution 
% information.

    properties

        % Edge length of cube.
        edge_length;
        % Extreme points covered by voxelized space.
        minimum_point;
        maximum_point;

        % Number of voxels in each dimension.
        size;

        % Cell array to contain indices of objects when
        % voxel is occupied.
        occupied_table;
        
        % Number of occupied cells.
        occupied_count = 0;

        % Number of unique objects in voxelized space.
        object_count = 0;

        % Indices of unique objects.
        unique_objects = [];

        % Number of dimensions in space.
        dimensionality = 3;

    end
    
    properties(Access=private)
       
        % Index of last object checked. Used for optimization. 
        last_index_checked = inf;
        
    end


    methods

        function ob = CubeVoxelization(edge, minp, maxp)
        % Constructor of class object. 
        %
        % Inputs:
        % edge      Edge length of the voxels.
        % minp      Minimum point of voxelized space.
        % maxp      Maximum point of voxelized space.
        %
        % Number of cubes is computed from the inputs.

            % Set object properties.
            ob.edge_length = edge;
            ob.minimum_point = minp;

            % Compute number of required cubes and store.
            ob.size = ceil((maxp - minp)/edge);

            % Initialize table of contained elements.
            ob.occupied_table = cell(ob.size);

            % Recompute maximum point.
            ob.maximum_point = minp + ob.size*edge;

        end

        function cc = get_coordinates(ob, cen)
        % Convert Cartesian point <cen> into 3D cube coordinates.

            % Cube coordinates.
            cc = floor(bsxfun(@minus,cen,ob.minimum_point)...
                       ./ob.edge_length) + 1;
            %-

            % Map coordinates outside of space to nearest edge voxels.
            Upper = ob.size;
            Lower = ones(size(Upper)); %#ok<CPROPLC>
            
            % If coordinate is larger than existing limits, return a voxel
            % on the edge.
            cc = bsxfun(@max,cc,Lower);
            cc = bsxfun(@min,cc,Upper);

        end

        function objects = get_neighbor_objects(ob, cc, index)
        % Get object indices of objects in a given cube or any of its
        % 27 neighbour cubes.

            % Optimization when no or few objects in space.
            if ob.object_count == 0 || ...
                    (ob.object_count == 1 && nargin > 2 ...
                     && ob.unique_objects == index)
                %-
                
                objects = [];
                return;
            end

            % Get limits of sub cube.
            DimLimits = ob.neighbor_limits(cc);

            % Get object list, may contain dublicates.
            objects = horzcat(ob.occupied_table{...
                                    DimLimits(1,1):DimLimits(2,1),...
                                    DimLimits(1,2):DimLimits(2,2),...
                                    DimLimits(1,3):DimLimits(2,3)});
            %-

            % Remove dublicates and sort.
            objects = unique(objects);

            % If index defined, remove from object list,
            % otherwise return list as is.
            if nargin < 3
                return;
            end

            % Number of objects.
            NObj = length(objects);

            % If no objects, all is good.
            if NObj == 0
                return;

            % If only self, return empty.
            elseif NObj == 1 && objects == index
                objects = [];
                return;
            end

            % Make list with removed index.
            if objects(1) == index
                objects = objects(2:end);
            elseif objects(end) == index
                objects = objects(1:end-1);
            else

                % Remove own index.
                for iObj = 2:NObj - 1

                    if objects(iObj) == index
                        objects = [objects(1:iObj-1); objects(iObj+1:end)];
                    end

                end
            end

        end
        
        function add_object_by_cc(ob, cc, index)
        % Add an object/objects with a given cube coordinates and object
        % index to the voxelized space.
        
            % If index is given as a scalar when multiple objects,
            % replicate.
            if isscalar(index) && size(cc,1) > 1 %#ok<CPROPLC>
                index = repmat(index,size(cc,1),1); %#ok<CPROPLC>
            end

            % Add objects one at a time.
            NObj = size(cc,1); %#ok<CPROPLC>

            for iObj = 1:NObj

                % Cube coordinate of current object.
                cci = cc(iObj,:);

                % If the coordinate voxel is empty.
                if isempty(ob.occupied_table{cci(1),cci(2),cci(3)})
                    
                    % Initialize the list to have only the new element.
                    ob.occupied_table{cci(1),cci(2),cci(3)} = index(iObj);

                    % Increase count by one as new cell gets occupied.
                    ob.occupied_count = ob.occupied_count + 1;
                else
                    % Existing indices in cell.
                    ExtInd = ob.occupied_table{cci(1),cci(2),cci(3)};

                    % Add new object to list and keep it sorted.
                    [ExtInd, isNew] = CubeVoxelization.add2sortedlist(...
                                                   ExtInd, ...
                                                   index(iObj));
                    %-
                    
                    if isNew
                        % Store update object list to cell array.
                        ob.occupied_table{cci(1),cci(2),cci(3)} = ExtInd;
                    end
                end

                % If object was just checked for uniqueness, skip check.
                if index(iObj) == ob.last_index_checked
                    continue;

                % If object is the first, no need for check.
                elseif ob.object_count == 0
                    
                    % Set as only unique object.
                    ob.unique_objects = index(iObj);
                    
                    % Increase unique object count.
                    ob.object_count = 1;
                    
                    % Set last checked index.
                    ob.last_index_checked = index(iObj);
                    
                    continue;
                end

                % Update total unique object list as well.
                ExtInd = ob.unique_objects;

                % Add and keep sorted.
                [ExtInd, isNew] = CubeVoxelization.add2sortedlist(...
                                                        ExtInd,...
                                                        index(iObj));
                %-

                if isNew
                    % Update unique object list.
                    ob.unique_objects = ExtInd;
                    
                    % Increase unique object count.
                    ob.object_count = ob.object_count + 1;
                    
                    % Set last checked index.
                    ob.last_index_checked = index(iObj);
                end
            end
            
        end

        function cc = add_object_by_cen(ob, cen, index)
        % Add an object/objects with a given center and object index to the
        % voxelized space.

            % If point is out of space bounds, remove empty.
            if any(any(bsxfun(@gt,cen,ob.maximum_point))) || ...
               any(any(bsxfun(@lt,cen,ob.minimum_point)))
           
                % Return empty.
                cc = [];
                warning(['Point(s) out of voxel space bounds. '...
                         'Omitting all points.']);
                %-
                return;
            else

                % Get cube coordinates of given point.
                cc = ob.get_coordinates(cen);
                
                % Add object with a known cube coordinate.
                ob.add_object_by_cc(cc,index);
            end

        end

    end

    methods(Access=protected)

        function DimLimits = neighbor_limits(ob,cc)
            % Get index limits of neighboring voxels. If cc is not on the
            % edge, the resulting limits give the 27 neighboring voxels.
            % Less on the edge.

            % By default on all dimensions the index goes 
            % from cc-1 to cc+1.
            DimLimits = repmat([-1; 1],1,ob.dimensionality);

            % Check which dimensions reach the edge.
            IFirst = cc == 1;
            ILast  = cc == ob.size;
            
            % Set start/end index to self.
            DimLimits(1,IFirst) = 0;
            DimLimits(2,ILast)  = 0;

            % Add center to get final index limits.
            DimLimits = bsxfun(@plus,DimLimits,cc);

        end
    end
    
    methods(Static)
        
        function [List, isNew] = add2sortedlist(List, Element)
        % Add Element to List while keeping sorted order.

            isNew = true;
        
            % Add to beginning.
            if Element < List(1)
                List = [Element, List];
                
            % Add to end.
            elseif Element > List(end)
                List = [List, Element];
                
            % Add to in between.
            else
                % Find first larger element.
                iNew = find(List > Element,1,'first');
                
                % Element is new only if it is not in the list.
                if List(iNew-1) == Element
                    isNew = false;
                else
                    % Augment the list with new element.
                    List = [
                            List(1:iNew-1), ...
                            Element, ...
                            List(iNew:end)
                           ];
                    %-
                end
            end
                
        end
        
    end

end