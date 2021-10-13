close all
clear all
clc

% We load the parameters
[r, scaling_factor, gamma, d, A] = load_parameters();

% We obtain the average of the length and width of the buildings
[E_L, E_W, NumberBuildings] = obtain_averages_densities(A);

% Different vsalues of lambda for testing are defined
lambda = NumberBuildings*gamma./repmat(scaling_factor,length(r),1)/A;

%We obtain an upper bound of the CDF of the length of the ligtht
z = 0:0.01:1000;

% We preset the size for F_Z_ub
F_Z_ub = cell(length(r),length(scaling_factor));

% For each value of lambda, the CDFs are computed
for m = 1:length(r)
    for n = 1:length(scaling_factor)
        F_Z_ub{m,n}  = 1-exp(-lambda(m,n)*r(m)/2.*z);
    end
end

save('shadowing_intervals_analytic.mat','lambda','E_L','E_W','NumberBuildings','z','F_Z_ub');