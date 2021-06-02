classdef obstacle < handle                       % Class simulating an experiment with N nanoparticles linearly interacting with a single optical cavity field
  properties
    boundary_x                                   % Closed curve (x coords.) delimiting the boundaries of the obstacle
    boundary_y                                   % Closed curve (y coords.) delimiting the boundaries of the obstacle
    
    valid_inside                                 % If the region where the particle can not enter is the inside or outside of the obstacle
  end
  
  % Se o obstáculo for um poligono: https://www.mathworks.com/help/matlab/ref/inpolygon.html
  % Aonde intersecta o poligono:    https://www.mathworks.com/help/map/ref/polyxpoly.html
  % Se o obstáulo for um circulo:   https://www.mathworks.com/help/map/ref/linecirc.html
  
  methods
    function obj = obstacle(x_bound, y_bound, valid_in)
      % Class constructor for an obstacle delimited by the polygon given by the points (x_bound, y_bound)
      % The region if the obstacle is inside the polygon, then valid_in = true.
      
      % Make sure the boundary is a closed polygon
      if x_bound(end) ~= x_bound(1)
        x_bound(end) = x_bound(1);
      end
      if y_bound(end) ~= y_bound(1)
        y_bound(end) = y_bound(1);
      end
      
      % Make sure these arrays are row vectors (I will need this for InterX!)
      if iscolumn(x_bound)
        x_bound = x_bound.';
      end
      if iscolumn(y_bound)
        y_bound = y_bound.';
      end
      
      obj.boundary_x = x_bound;
      obj.boundary_y = y_bound;
      
      obj.valid_inside = valid_in;
    end
    
    function inside_boundary = was_penetrated(obj, xq, yq)
      
      inside_polygon = inpolygon(xq, yq, obj.boundary_x, obj.boundary_y);
      
      inside_boundary = (inside_polygon == obj.valid_inside);
      
    end
    
    function [x_correct, y_correct] = reflect(obj, x_old, y_old, x_wrong, y_wrong)
      
      line = [[x_old, x_wrong];
              [y_old, y_wrong]];
      
      boundary = [obj.boundary_x; obj.boundary_y];
      
      intersection = InterX(line, boundary);
      
    end
    
  end  % End methods
end  % End classdef

% pgon = polyshape(obj.boundary_x, obj.boundary_y);
% plot(pgon)
% hold on
% plot([x_old, x_wrong],[y_old, y_wrong])
% plot(intersection(1,1), intersection(2,1), 'k', 'Marker', '*')

% x_line = [x_old, x_wrong];
      % y_line = [y_old, y_wrong];
      % Following line requires the Mapping Toolbox. I am moving to the InterX from FileExchange !
      % [x_intersect, y_intersect, idx] = polyxpoly(x_line, y_line, obj.boundary_x, obj.boundary_y);



