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

function h = plot_cylinders(ob,varargin)
% Draw cylinders of the QSM. The following extra parameters are
% available for customization, passed in label-value pairs:
%
% 'Points'      Number of vertices in cylinder rings.
%
% 'Curve'       Taper curve passed to the CYLINDER command.
%
% 'Color'       Color for the cylinders.
%

cdata = [];
curve = 1;
points = 20;

% Process additional arguments.
i = 1;
while i <= size(varargin,2)

    if ischar(varargin{i})

        switch lower(varargin{i})
            case 'points'
                points = varargin{i+1};
                i = i + 1;
            case 'curve'
                curve = varargin{i+1};
                i = i + 1;
            case 'color'
                cdata = varargin{i+1};
                i = i + 1;
            otherwise
                disp(['Ignoring unknown parameter: ' varargin{i}]);
        end
    end
    i = i + 1;
end

% Generate vertices and faces for a unit cylinder with a given number of
% vertices in each ring and a given taper curve.
[x, y, z] = cylinder(curve,points);

% Number of vertices in cylinder.
s = size(x);

% Number of cylinders in QSM.
n = ob.block_count;

% Reshape vertices.
x = repmat(x(:)',n,1);
y = repmat(y(:)',n,1);
z = repmat(z(:)',n,1);

% Scale to adjust correct radius and length.
x = bsxfun(@times,x,ob.cylinder_radius);
y = bsxfun(@times,y,ob.cylinder_radius);
z = bsxfun(@times,z,ob.cylinder_length);

u = [0; 0; 1];

% Origins of cylinders.
start_points = ob.cylinder_start_point;

% Directions of cylinders.
ax = ob.cylinder_axis;

% Find rotation angle to rotate the cylinders.
raxis = (bsxfun(@cross,u,ax'))';
angle = acos(ax(:,3)./sqrt(sum(ax.^2,2)));

% Store figure hold status.
set_hold_off = ~ishold;

% Variable to hold handles to graphical objects.
h = zeros(n,1);

% Iterate over cylinders.
for i = 1:n
    
    % Set hold on after plotting first cylinder.
    if i == 2
        hold on
    end
   
    % If the cylinder is not parallel to z-direction. Rotate vertices by
    % previously computed angle.
    if any(raxis(i,:))
        X = [x(i,:)' y(i,:)' z(i,:)'] ...
          * rotation_matrix(raxis(i,:)...
          / norm(raxis(i,:)),angle(i))';
    %- 
    
    % Downward pointing cylinder.
    elseif ax(i,3) == -1
        X = [x(i,:)' y(i,:)' -z(i,:)'];
    % Upward pointing cylinder.
    else
        X = [x(i,:)' y(i,:)' z(i,:)'];
    end

    % Reshape to size which is required by the surf-command.
    xi = reshape(X(:,1),s);
    yi = reshape(X(:,2),s);
    zi = reshape(X(:,3),s);

    % Transfer computed cylinder to the starting point of the branch.
    xi = xi + start_points(i,1);
    yi = yi + start_points(i,2);
    zi = zi + start_points(i,3);

    % Draw the surface of the branch.
    if ~isempty(cdata)
        h(i) = surf(xi,yi,zi,repmat(cdata(i,:),size(zi,1),size(zi,2)));
    else
        h(i) = surf(xi,yi,zi);
    end
    
end

% Set axis scaling equal.
axis equal;

% Restore the hold status if necessary.
if set_hold_off
    hold off;
end
