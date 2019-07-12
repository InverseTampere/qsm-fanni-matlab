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

classdef LeafModel < handle
% LeafModel is an abstract class for holding leaf information. Leaf
% properties describing location, scale and orientation of included leaves
% are included as properties. LeafModel is a subclass of <handle>.

	properties

		% Number of included leaves.
		leaf_count = 0;
        
        % Total leaf area of included leaves.
        leaf_area = 0;

		% Parameters of included leaves.
        
        % Origin of leaf.
        leaf_start_point = zeros(0,3);
        % Scale (x,y,z).
        leaf_scale = zeros(0,3);
        % Direction from origin to leaf tip.
        leaf_direction = zeros(0,3);
        % Direction of leaf normal.
        leaf_normal = zeros(0,3);
        
        % Index of parent block.
        leaf_parent = zeros(0,1);
        
        % Start point of petiole (on block surface).
        twig_start_point = zeros(0,3);

	end

	methods

		function ob = LeafModel()
            % Constructor.

        end

	end

	methods (Abstract)

        % Add single leaf with given parameters.
		add_leaf(ob, origin, dir, normal, scale, parent, twig, tris)
        
        % Convert leaves to triangles for intersection computations with
        % blocks.
		triangles(ob, origin, dir, normal, scale)
        
        % Plot accepted leaves.
		plot_leaves(ob, varargin)

		% Check if leaf with given parameters intersects leafs with 
		% given indices <Candidates>.
		leaf_intersect(ob, origin, dir, normal, scale, Candidates)

        % Compute the minimum and maximum corner of the model.
        bounding_box(ob)

    end


end