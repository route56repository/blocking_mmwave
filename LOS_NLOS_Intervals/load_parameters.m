function [r, scaling_factor, gamma, d, A] = load_parameters()
% We define the distance from the trajectory to the BS r that are analysed
r = [50; 150; 300];

% We define the differents scaling factors that are used in the simulation
% code so as to coorectly compute the analytic lambda
scaling_factor = [1 2 4];

% We define a corrective factor called gamma to later on compute the
% density lambda. For each radius, a gamma factor is defined
gamma = [2; 1.3; 1];
gamma = repmat(gamma,1,length(scaling_factor));

% We define the dimensions of the simulation deployment so as to have a
% introduce the area in the following function
d = 1000;
A = d^2;

end

