function [E_L, E_W, NumberBuildings] = obtain_averages_densities(A)

load('ChicagoDowntownBuildingStats.mat')

x_coord = [];
y_coord = [];

% We loop over all the buildings of the database
for n_building = 1:NumberBuildings
    % We obtain the number of edges that form the building n_building
    N_edges = sum(BuildingXYCoord(n_building,:,1)~=0)-1;
    % Then, we loop over all the edges that form a given building
    for n_edge = 1:N_edges
        % We obtain the 2 values of x_coord that define an edge
        x_coord = [x_coord; BuildingXYCoord(n_building,n_edge,1) BuildingXYCoord(n_building,n_edge+1,1)];
        % We obtain the 2 values of y_coord that define an edge
        y_coord = [y_coord; BuildingXYCoord(n_building,n_edge,2) BuildingXYCoord(n_building,n_edge+1,2)];
    end
end

% We obtain the angles of each segment with the current coordinates system
alpha = mod(angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)-(x_coord(:,1)+1i*y_coord(:,1)))))),pi);

alpha_aux = alpha;

% We obtain the angles of each segment with the current coordinates system
alpha_aux(alpha_aux>deg2rad(10),:) = [];

% We obtain the reference angle to which all segments are turned
alpha_ref = mean(alpha_aux);

% The modulus from the origin to the first point of the segment
mod_1 = sqrt(x_coord(:,1).^2+y_coord(:,1).^2);
mod_2 = sqrt(x_coord(:,2).^2+y_coord(:,2).^2);

% Now, the angles of all the points with respect of the origin are
% computed
alpha_1 = angle(exp(1i*(angle(x_coord(:,1)+1i*y_coord(:,1)))));
alpha_2 = angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)))));

% Then, the scenario is rotated
x_coord = [mod_1.*cos(alpha_1-alpha_ref), mod_2.*cos(alpha_2-alpha_ref)];
y_coord = [mod_1.*sin(alpha_1-alpha_ref), mod_2.*sin(alpha_2-alpha_ref)];

% Once all the segments have been rotated accordingly, the mean average and
% widths of the blockages are obtained
% We obtain the angles of each segment with the current coordinates system
alpha = angle(exp(1i*(angle(x_coord(:,2)+1i*y_coord(:,2)-(x_coord(:,1)+1i*y_coord(:,1))))));

% Finally, we obtain the averages of both the lengths and widths of the
% lines (width is understood as the length of the vertical lines) and their
% densities
L = sqrt((x_coord((alpha<deg2rad(10) | alpha>deg2rad(170)),1)-x_coord((alpha<deg2rad(10) | alpha>deg2rad(170)),2)).^2+(y_coord((alpha<deg2rad(10) | alpha>deg2rad(170)),1)-y_coord((alpha<deg2rad(10) | alpha>deg2rad(170)),2)).^2);
E_L = mean(L);

W = sqrt((x_coord((alpha<deg2rad(100) & alpha>deg2rad(80)),1)-x_coord((alpha<deg2rad(100) & alpha>deg2rad(80)),2)).^2+(y_coord((alpha<deg2rad(100) & alpha>deg2rad(80)),1)-y_coord((alpha<deg2rad(100) & alpha>deg2rad(80)),2)).^2);
E_W = mean(W);

end

