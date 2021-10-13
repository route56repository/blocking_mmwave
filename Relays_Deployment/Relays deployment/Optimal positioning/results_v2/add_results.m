load('P_allKO_no_sensitivity_simulation_a.mat');
load('P_allKO_no_sensitivity_simulation_b.mat');
load('P_allKO_no_sensitivity_simulation_c.mat');
load('P_allKO_no_sensitivity_simulation_d.mat');
load('P_allKO_no_sensitivity_simulation_e.mat');

P_allKO_no_sensitivity_simulation = (P_allKO_no_sensitivity_simulation_a + P_allKO_no_sensitivity_simulation_b + P_allKO_no_sensitivity_simulation_c + P_allKO_no_sensitivity_simulation_d + P_allKO_no_sensitivity_simulation_e)/5;

load('P_allKO_sensitivity_simulation_a.mat');
load('P_allKO_sensitivity_simulation_b.mat');
load('P_allKO_sensitivity_simulation_c.mat');
load('P_allKO_sensitivity_simulation_d.mat');
load('P_allKO_sensitivity_simulation_e.mat');

P_allKO_sensitivity_simulation = (P_allKO_sensitivity_simulation_a + P_allKO_sensitivity_simulation_b + P_allKO_sensitivity_simulation_c + P_allKO_sensitivity_simulation_d + P_allKO_sensitivity_simulation_e)/5;

save('P_allKO_no_sensitivity_simulation.mat','P_allKO_no_sensitivity_simulation');
save('P_allKO_sensitivity_simulation.mat','P_allKO_sensitivity_simulation');



