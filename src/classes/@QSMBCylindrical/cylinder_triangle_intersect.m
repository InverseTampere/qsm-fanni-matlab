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

function intersect = cylinder_triangle_intersect(ob,jCyl,tri,edgehit,debug)
% Detect the intersection between a cylinder in the QSM and a triangle. The
% cylinder is identified by an index. The triangle is defined by three
% vertices (1 x 9) double. The <edgehit> boolean input can be used to
% include (true) or exclude (false) edgehits. A debug flag is also included
% for development purposes.
%
% This MATLAB implementation is loosely based on the algorithm presented by
% Dave Eberly in Geometric Tools (www.geometrictools.com/).

    % Default debug status is false.
    if nargin < 5
        debug = false;
    end
    
    % Initialize return value as false, i.e., no hit occurs.
    intersect = false;

    % Length and radius of given cylinder.
    h  = ob.cylinder_length(jCyl);
    r  = ob.cylinder_radius(jCyl);
    
    % Half of the cylinder height.
    h2 = h/2;

    % Optimization step, where the coordinate system can be pre-computed 
    % and stored in the object.
    if isnan(ob.cylinder_coordinate_system(1,1,jCyl))
        
        % Cylinder axis direction.
        ax = ob.cylinder_axis(jCyl,:);

        % Center point of cylinder.
        C = ob.cylinder_mid_point(jCyl,:);

        % z-axis of the coordinate system.
        D = ax;
        
        e3 = [0 0 1];

        % Form x- and y-axes.
        if all(D == e3)
            U = [1 0 0];
            V = [0 1 0];
        else
            U = cross(D,e3);
            U = U./norm(U);

            V = cross(U,D);
        end

        % Coordinate sysyem matrix.
        E = [U; V; D];
        
        % Store for later use.
        ob.cylinder_coordinate_system(:,:,jCyl) = E;
    else
        % Recall stored coordinate system matrix.
        E = ob.cylinder_coordinate_system(:,:,jCyl);
    end
    
    % Points of triangle.
    P = reshape(tri,3,3)';
    
    % Triangle vertices are given with respect to the cylinder midpoint.
    PC = bsxfun(@minus,P,C);
    
    % Convert points to cylinder coordinate system.
    Q = PC*E';
    
    % Sort triangle vertices by height in the cylinder CS.
    [~,Iz] = sort(Q(:,3));
    Q = Q(Iz,:);
    
    % Shorthand for the three vertices in the cylinder CS.
    [x,y,z] = deal(Q(:,1),Q(:,2),Q(:,3));
    
    % All points below or above cylinder.
    if z(3) < -h2 || z(1) > h2
        if debug
            disp('Case 0a, 0b');
        end

        return;
    end
    
    % All points vertically between cylinder slabs.
    if z(1) >= -h2 && z(3) <= h2
        if debug
            disp('Case 3a');
        end

        intersect = polygon_check(Q(:,1:2),r);
        return;
    end
    
    % xy-vectors of coordinates.
    Q0 = Q(1,1:2);
    Q1 = Q(2,1:2);
    Q2 = Q(3,1:2);
    
    if z(1) < -h2
        if z(3) > h2
            
            % One or two points above cylinder top and one under cylinder
            % bottom.
            if z(2) >= h2
                if debug
                    disp('Case 4a/4b');
                end

                
                n0 = -h2 - z(1);
                n1 =  h2 - z(1);
                
                denom0 = 1/(z(2) - z(1));
                denom1 = 1/(z(3) - z(1));
                
                p = zeros(4,2);
                
                t = n0*denom1;
                p(1,:) = Q0 + t*(Q2 - Q0);
                t = n0*denom0;
                p(2,:) = Q0 + t*(Q1 - Q0);
                t = n1*denom0;
                p(3,:) = Q0 + t*(Q1 - Q0);
                t = n1*denom1;
                p(4,:) = Q0 + t*(Q2 - Q0);
                
                intersect = polygon_check(p,r);
                
            % One or two points below cylinder bottom and one above the
            % top.
            elseif z(2) <= -h2
                if debug
                    disp('Case 4c/4d');
                end

                
                n0 = -h2 - z(3);
                n1 =  h2 - z(3);
                
                denom0 = 1/(z(2) - z(3));
                denom1 = 1/(z(1) - z(3));
                
                p = zeros(4,2);
                
                t = n0*denom1;
                p(1,:) = Q2 + t*(Q0 - Q2);
                t = n0*denom0;
                p(2,:) = Q2 + t*(Q1 - Q2);
                t = n1*denom0;
                p(3,:) = Q2 + t*(Q1 - Q2);
                t = n1*denom1;
                p(4,:) = Q2 + t*(Q0 - Q2);
                
                intersect = polygon_check(p,r);
                
            % One point above, one below, one between.
            else
                if debug
                    disp('Case 5');
                end

                
                n0 = -h2 - z(1);
                n1 =  h2 - z(1);
                
                denom0 = 1/(z(2) - z(1));
                denom1 = 1/(z(3) - z(1));
                
                p = zeros(5,2);
                
                t = n0*denom1;
                p(1,:) = Q0 + t*(Q2 - Q0);
                t = n0*denom0;
                p(2,:) = Q0 + t*(Q1 - Q0);
                p(3,:) = Q1;
                t = n1*denom0;
                p(4,:) = Q0 + t*(Q1 - Q0);
                t = n1*denom1;
                p(5,:) = Q0 + t*(Q2 - Q0);
                
                intersect = polygon_check(p,r);
            end
            
        elseif z(3) > -h2
            
            % One point inside, one or two below.
            if z(2) <= -h2
                if debug
                    disp('Case 3b/3c');
                end

                
                n0 = -h2 - z(3);
                
                denom0 = 1/(z(2) - z(3));
                denom1 = 1/(z(1) - z(3));
                
                p = zeros(3,2);
                
                t = n0*denom0;
                p(1,:) = Q2 + t*(Q1 - Q2);
                t = n0*denom1;
                p(2,:) = Q2 + t*(Q0 - Q2);
                p(3,:) = Q2;
                
                intersect = polygon_check(p,r);
                
            % Two points inside, one below.
            else
                if debug
                    disp('Case 4e');
                end

                
                n0 = -h2 - z(1);
                
                denom0 = 1/(z(3) - z(1));
                denom1 = 1/(z(2) - z(1));
                
                p = zeros(4,2);
                
                t = n0*denom0;
                p(1,:) = Q0 + t*(Q2 - Q0);
                t = n0*denom1;
                p(2,:) = Q0 + t*(Q1 - Q0);
                p(3,:) = Q1;
                p(4,:) = Q2;
                
                intersect = polygon_check(p,r);
            end
        else 
            
            % Two points below, one on lower bound.
            if z(2) < -h2 
                if debug
                    disp('Case 1a');
                end

                if edgehit
                    intersect = x(3)^2 + y(3)^2 <= r;
                end
                return;
                
            % One point below, two on bottom bound.
            else
                if debug
                    disp('Case 2a');
                end

                if edgehit
                    Q1 = Q(3,1:2);
                    Q2 = Q(2,1:2);
                    
                    intersect = line_segment_check(Q1,Q2,r);
                end
                return;
            end
        end
            
    elseif z(1) < h2
        
        % One point inside, one or two point above.
        if z(2) >= h2
            if debug
                disp('Case 3d/3e');
            end

            
            n0 = -h2 - z(1);
                
            denom0 = 1/(z(3) - z(1));
            denom1 = 1/(z(2) - z(1));

            p = zeros(3,2);

            t = n0*denom0;
            p(1,:) = Q0 + t*(Q2 - Q0);
            t = n0*denom1;
            p(2,:) = Q0 + t*(Q1 - Q0);
            p(3,:) = Q0;

            intersect = polygon_check(p,r);
            
        % Two points inside, one above.
        else
            if debug
                disp('Case 4f');
            end

            
            n0 = -h2 - z(3);
                
            denom0 = 1/(z(2) - z(3));
            denom1 = 1/(z(1) - z(3));

            p = zeros(4,2);

            t = n0*denom0;
            p(1,:) = Q2 + t*(Q1 - Q2);
            t = n0*denom1;
            p(2,:) = Q2 + t*(Q0 - Q2);
            p(3,:) = Q0;
            p(4,:) = Q1;

            intersect = polygon_check(p,r);
        end
    else
        
        % One point on top bound.
        if z(2) > h2
            if debug
                disp('Case 1b');
            end

            if edgehit
                intersect = x(1)^2 + y(1)^2 <= r;
            end
            return;
            
        % Two point on top bound.
        else
            if debug
                disp('Case 2b');
            end

            if edgehit
                Q1 = Q(2,1:2);
                Q2 = Q(1,1:2);
                
                intersect = line_segment_check(Q1,Q2,r);
            end
            return;
        end
    end
    
end

function intersect = polygon_check(p,r)
% Check if polygon with vertices <p> contains the origin or if any polygon
% edge is closer than <r> to the origin.
        
    intersect = false;

    % If polygon contains origin => intersect.
    if polyarea(p(:,1),p(:,2)) > 0 && inpolygon(0,0,p(:,1),p(:,2))
        intersect = true;
    else
        NEdge = size(p,1);

        % Check edge distance to origin.
        for iEdge = 1:NEdge
            Q1 = p(iEdge,:);
            if iEdge < NEdge
                Q2 = p(iEdge+1,:);
            else
                Q2 = p(1,:);
            end
            
            % Check distance.
            intersect = line_segment_check(Q1,Q2,r);
            
            if intersect
                return;
            end
        end
    end

end

function intersect = line_segment_check(a,b,r)
% Check if line segment with end points <a> and <b> is closer than <r> to
% the origin.

    if dot(a-b,-b)*dot(b-a,-a) >= 0
        intersect = abs(det([a 1; b 1; 0 0 1])) / norm(a-b) < r;
    else
        intersect = min(norm(a),norm(b)) < r;
    end

end