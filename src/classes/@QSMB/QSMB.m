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

classdef QSMB
% QSMB - Quantitative Structure Model Blocks. Abstract class for holding
% quantitative structure information.

    properties

        % Number of geometric primitives (blocks).
        block_count = 0;

        % Minimum and maximum point of axis-aligned bounding box.
        tree_limits = [];
        
        % Function handle to a function that defines a twig parameter
        % distribution. The inputs are the parameters of the blocks.
        fun_twig_distribution;
        
        % Lower and upper limit for twig length.
        twig_length_limits = [];
    end


    methods
        
        % Constructor.
        function ob = QSMB()
            
        end
        
    end
    
    methods (Abstract)
        
        % Get selected properties of the blocks of the model.
        Props = get_block_properties(ob,PropNames,vargin)
        
        % Generate twigs connecting leaves to the blocks.
        [TwigStart, TwigEnd] = generate_twigs(ob, LeafParent)

        % Detect intersection between a block and a triangle.
        IntersectFlag = block_triangle_intersection(ob, JBlock, Tris)

        % Convert QSM into simple voxelization that shows which voxels are
        % occupied by any part of the blocks.
        BlockVoxelization = toVoxels(ob, edge, minp, maxp)        

    end

end