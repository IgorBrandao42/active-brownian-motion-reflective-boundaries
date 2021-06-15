classdef obstacle < handle                       % Class simulating an experiment with N nanoparticles linearly interacting with a single optical cavity field
  properties
    boundary_x                                   % Closed curve (x coords.) delimiting the boundaries of the obstacle
    boundary_y                                   % Closed curve (y coords.) delimiting the boundaries of the obstacle
    normal_vec                                   % Unitary normal vector for each line segment of the boundary
    N_sides                                      % Number of sides of the boundary
    interior_is_inside                           % If the region where the particle can not enter is the inside or outside of the obstacle
    
    color                                        % Color for plotting  
    alpha                                        % How intense is the color
  end
  
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
      
      dx = diff(obj.boundary_x);
      dy = diff(obj.boundary_y);
      
      norm_vec = sqrt(dx.^2 + dy.^2);
      
      normal_x =  dy./norm_vec;
      normal_y = -dx./norm_vec;
      
      obj.normal_vec = [normal_x; normal_y];
      
      obj.N_sides = size(obj.normal_vec, 2);
      
      obj.interior_is_inside = valid_in;
      
      if obj.interior_is_inside
        obj.normal_vec = - obj.normal_vec;
        obj.color = "red";
        obj.alpha = 1;
      else
        obj.color = "blue";
        obj.alpha = 0.05;
      end
    end
    
    function inside_boundary = was_penetrated(obj, xq, yq)
      
      particle_is_inside_polygon = inpolygon(xq, yq, obj.boundary_x, obj.boundary_y);
      
      inside_boundary = (particle_is_inside_polygon == obj.interior_is_inside);
      
    end
    
    function [x_correct, y_correct] = reflect(obj, particle, k)
      x_old   = particle.x(k);
      y_old   = particle.y(k);
      x_wrong = particle.x(k+1);
      y_wrong = particle.y(k+1);
      
      line_x = [x_old; x_wrong];
      line_y = [y_old; y_wrong];
      
      [p_x, p_y, ~, idx] = intersections(line_x, line_y, obj.boundary_x, obj.boundary_y);
      
      if isempty(idx)
        disp("Aconteceu aquele erro, tô fingindo que não vi!")
        x_correct = x_old;
        y_correct = y_old;
        return
      end
      
      idx = floor(idx(1,1));
      
      r_wrong = [x_wrong; y_wrong];
      p = [p_x(1); p_y(1)];
      n = obj.normal_vec(:, idx);
      
      r_correct = r_wrong - 2*dot(r_wrong - p, n)*n;
      
      x_correct = r_correct(1);
      y_correct = r_correct(2);
      
      if obj.was_penetrated(x_correct, y_correct)                       % After the reflection, the particle can leak through another segment of the boundary !
        line_x = [x_old; x_correct];                                    % Instead of infinite reflections, just stop it at the intersection with this other segment
        line_y = [y_old; y_correct];
        
        [p_x, p_y, ~, ~] = intersections(line_x, line_y, obj.boundary_x, obj.boundary_y);
        
        if isempty(p_x)
          disp("Aconteceu aquele OUTRO erro, tô fingindo que não vi tbm!")
          x_correct = x_old + (x_old - x_wrong);
          y_correct = y_old + (y_old - y_wrong);
          return
        end
        
        x_correct = p_x(1);
        y_correct = p_y(1);
      end
      
      
%       if isempty(idx)
%         x_correct = 2;
%       end
      
%       if obj.was_penetrated(x_correct, y_correct)
%         obj.show()
%         hold on
%         plot([x_old, x_wrong],[y_old, y_wrong])
%         plot(x_old,y_old, 'k*')
%         plot(x_correct, y_correct, 'b*')
%         disp("Boundary is too skew and resolution too low!")
%       end
    end
    
    function h = show(obj)
      pgon = polyshape(obj.boundary_x, obj.boundary_y);
      h = plot(pgon, 'FaceColor', obj.color, 'FaceAlpha', obj.alpha);
      
%       hold on
%       
%       mean_side_x = (obj.boundary_x(2:end) + obj.boundary_x(1:end-1) )/2;
%       mean_side_y = (obj.boundary_y(2:end) + obj.boundary_y(1:end-1) )/2;
%       quiver(mean_side_x, mean_side_y, obj.normal_vec(1,:), obj.normal_vec(2,:))
      
%       hold off
    end
    
  end  % End methods
end  % End classdef



% function [x_correct, y_correct] = reflect_OLD(obj, x_old, y_old, x_wrong, y_wrong)
%       
%       line_x = [x_old; x_wrong];
%       line_y = [y_old; y_wrong];
%       
%       [p_x, p_y, ~, idx] = intersections(line_x, line_y, obj.boundary_x, obj.boundary_y);
%       
%       if isempty(idx)
%         disp("Aconteceu aquele erro, tô fingindo que não vi!")
%         x_correct = x_old;
%         y_correct = y_old;
%         return
%       end
%       
%       idx = floor(idx(1,1));
%       
%       r_wrong = [x_wrong; y_wrong];
%       p = [p_x(1); p_y(1)];
%       n = obj.normal_vec(:, idx);
%       
%       r_correct = r_wrong - 2*dot(r_wrong - p, n)*n;
%       
%       x_correct = r_correct(1);
%       y_correct = r_correct(2);
%       
%       if obj.was_penetrated(x_correct, y_correct)                       % After the reflection, the particle can leak through another segment of the boundary !
%         line_x = [x_old; x_correct];                                    % Instead of infinite reflections, just stop it at the intersection with this other segment
%         line_y = [y_old; y_correct];
%         
%         [p_x, p_y, ~, ~] = intersections(line_x, line_y, obj.boundary_x, obj.boundary_y);
%         
%         if isempty(p_x)
%           disp("Aconteceu aquele OUTRO erro, tô fingindo que não vi tbm!")
%           x_correct = x_old + (x_old - x_wrong);
%           y_correct = y_old + (y_old - y_wrong);
%           return
%         end
%         
%         x_correct = p_x(1);
%         y_correct = p_y(1);
%       end
%       
%       
% %       if isempty(idx)
% %         x_correct = 2;
% %       end
%       
% %       if obj.was_penetrated(x_correct, y_correct)
% %         obj.show()
% %         hold on
% %         plot([x_old, x_wrong],[y_old, y_wrong])
% %         plot(x_old,y_old, 'k*')
% %         plot(x_correct, y_correct, 'b*')
% %         disp("Boundary is too skew and resolution too low!")
% %       end
%     end
    

 
% line = [[x_old, x_wrong];
%   [y_old, y_wrong]];
% 
% boundary = [obj.boundary_x; obj.boundary_y];
% 
% p_old = InterX(line, boundary); % Intersection between last trajectory and boundary

% Use InterX to find intersection of boundary and last trajectory step
% line = [[x_old, x_wrong];
%   [y_old, y_wrong]];
% 
% boundary = [obj.boundary_x; obj.boundary_y];
% 
% p = InterX(line, boundary); % Intersection between last trajectory and boundary



% idx = obj.which_boundary_side_is_in(p_x, p_y);
% function idx = which_boundary_side_is_in(p_x, p_y)
% for idx=1:obj.N_sides                      % For each side
%   m = ;                                    % Find its line equation y = m*x + b
%   b = ;
%   
%   if m*p_x + b == p_y
%     return
%   end
%   
%   error("Could not find intersection on boundary!")
%   
% end
% end



% pgon = polyshape(obj.boundary_x, obj.boundary_y);
% plot(pgon)
% hold on
% plot([x_old, x_wrong],[y_old, y_wrong])
% plot(intersection(1,1), intersection(2,1), 'k', 'Marker', '*')

% x_line = [x_old, x_wrong];
      % y_line = [y_old, y_wrong];
      % Following line requires the Mapping Toolbox. I am moving to the InterX from FileExchange !
      % [x_intersect, y_intersect, idx] = polyxpoly(x_line, y_line, obj.boundary_x, obj.boundary_y);



