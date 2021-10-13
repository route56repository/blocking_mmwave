function [x_coord, y_coord] = rotate(x_coord, y_coord, alpha_0)
% The modulus from the origin to the first point of the segment
mod_1 = sqrt(x_coord(:,1).^2+y_coord(:,1).^2);
mod_2 = sqrt(x_coord(:,2).^2+y_coord(:,2).^2);

% Now, the angles of all the points with respect of the origin are
% computed
alpha_1 = angle(exp(1i*(angle(x_coord(:,1)+1i*y_coord(:,1)))));
alpha_2 = angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)))));

% Then, the scenario is rotated
x_coord = [mod_1.*cos(alpha_1-alpha_0), mod_2.*cos(alpha_2-alpha_0)];
y_coord = [mod_1.*sin(alpha_1-alpha_0), mod_2.*sin(alpha_2-alpha_0)];

end

