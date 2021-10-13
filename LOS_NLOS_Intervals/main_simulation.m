close all
clear all
clc

% We load the parameters
[r, scaling_factor, gamma, d, ~] = load_parameters();

% We define the distance of the grid of positions where the BS is placed
d_B = 1;

S_emp = cell(length(r),length(scaling_factor));
Z_emp = cell(length(r),length(scaling_factor));
N_z = zeros(length(r),length(scaling_factor));
N_s = zeros(length(r),length(scaling_factor));
d_emp = zeros(length(r),length(scaling_factor));
E_Z_emp = zeros(length(r),length(scaling_factor));
E_S_emp = zeros(length(r),length(scaling_factor));

N_scenarios = 1;

% We loop over N_scenarios generated scenarios
for n = 1:length(scaling_factor)
    for n_scenario = 1:N_scenarios
        % Each time that a simulation runs, it means that some buildings of the
        % scenario are imported with the decimation specified but in a random basis
        [N_segments, x_coord, y_coord] = import_decimate_align(scaling_factor(n));
        
        % We consider the scenario oriented with no rotation anb reversed with
        % 180 degrees of rotation
        for alpha = [0 90 180 270]
            % The scenario is rotated in steps of alpha degrees
            [x_coord, y_coord] = rotate(x_coord, y_coord, deg2rad(alpha));
            
            % Once the scenario is generated accordingly, the BS is placed and
            % moved where desired
            for y_B = -d/2:d_B:d-min(r)
                for x_B = -d/2:d_B:d/2
                    % We define the distance to the trajectory r for each step of the
                    % simulation
                    for m = 1:length(r)
                        % We only compute the lights and shadows for a
                        % trajectory if this trajectory is within the limits of
                        % the deployment
                        if y_B+r(m) <= d/2
                            % Once that the trajectory and the BS are placed,
                            % and the trajectory is within of the deployment,
                            % the limits of the intervals of lights and shadows
                            % are computed
                            [x_min_light, x_max_light, x_min_shadow, x_max_shadow] = obtain_lights_shadows(x_coord, y_coord, x_B, y_B, r(m), d);
                            
                            %We add the lights and the shadows to their corresponding vectors
                            Z_emp{m,n} = [Z_emp{m,n} x_max_light-x_min_light];
                            S_emp{m,n} = [S_emp{m,n} x_max_shadow-x_min_shadow];
                            
                            %We obtain the number of lights and shadows in the
                            %simulation
                            N_z(m,n) = N_z(m,n) + length(x_min_light);
                            N_s(m,n) = N_s(m,n) + length(x_min_shadow);
                            
                            % Then, we accummulate the total length of th
                            E_Z_emp(m,n) = E_Z_emp(m,n) + sum(x_max_light-x_min_light);
                            E_S_emp(m,n) = E_S_emp(m,n) + sum(x_max_shadow-x_min_shadow);
                            if ~isempty(x_max_shadow) && ~isempty(x_max_light)
                                d_emp(m,n) = d_emp(m,n) + max([x_max_light(end),x_max_shadow(end)]) - min([x_min_light(1) x_min_shadow(1)]);
                            elseif ~isempty(x_max_shadow)
                                d_emp(m,n) = d_emp(m,n) + x_max_shadow(end) - x_min_shadow(1);
                            elseif ~isempty(x_max_light)
                                d_emp(m,n) = d_emp(m,n) + x_max_light(end) - x_min_light(1);
                            end
                        end
                    end
                end
                alpha
                y_B
            end
        end
    end
end
% Once all the tests are carried out, the CDFs are computed
% We preset the size of the cell arrays
F_Z_emp = cell(length(r),length(scaling_factor));
z_emp = cell(length(r),length(scaling_factor));
F_S_emp = cell(length(r),length(scaling_factor));
s_emp = cell(length(r),length(scaling_factor));

for m = 1:length(r)
    for n = 1:length(scaling_factor)
        [F_Z_emp{m,n}, z_emp{m,n}] = ecdf(Z_emp{m,n});
        [F_S_emp{m,n}, s_emp{m,n}] = ecdf(S_emp{m,n});
    end
end

% We compute the average length of the shadow and light intervals
E_Z_emp = E_Z_emp./N_z;
E_S_emp = E_S_emp./N_s;

% We compute the average of the density of the number of shadow and light
% intervals
E_N_z_density_emp = N_z./d_emp;
E_N_s_density_emp = N_s./d_emp;

save('shadowing_intervals_empiric.mat','F_Z_emp','z_emp','F_S_emp','s_emp','E_Z_emp','E_S_emp','E_N_z_density_emp','E_N_s_density_emp');