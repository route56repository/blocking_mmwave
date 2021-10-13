clear;

%In this script, we try to simulate the case in which the buildings are
%modeled as rectangles with height. We only evaluate the link between the
%BS and the UE.
%THE VECTORS THAT DEFINE THE COORDINATES ARE COLUMN VECTORS OF SIZE 3x1.

%The radius of the cell is 100 m.
R = 300;

%We define some differentials that will be used to evaluate the
%granularities in which the UE positions are analised.
d_d = 0.05*R;

%Base station at the origin. Its height is 40 m.
x_B_1 = 0;
x_B_2 = 0;
H_B = 40;
x_B = [x_B_1; x_B_2; H_B];

%We define the vector where the users are going to be located.
d = 0:d_d:R;
phi = 0;
pos = d*cos(phi)+1i*d*sin(phi);
num_pos = length(pos);
H_U = 1.5;

x_U = [real(pos); imag(pos); H_U*ones(1,num_pos)];

%Length l is a r.v. such that: l~U[0,L_max].
L_max = 30;
%Width w is a r.v. such that: w~U[0,W_max].
W_max = 30;
%Hieght h is a r.v. such that: h~U[0,H_max].
H_max = 30;
%Orientation \theta is a r.v. such that: \theta~U[0,\Theta_max].
Theta_max = pi;

%We make simulations for different blockage densities.
lambda = 1e-4;

N_runs = 4000;

P_allKO_simulation = zeros(1,num_pos);

  
for n_run = 1:N_runs
    %We calculate both the number and the characteristics of the blocking
    %elements.
    B = poissrnd(lambda*pi*(1.2*R)^2); %It calculates the number of buildings of
    %every simulation, according to the fact
    %that B is a Poisson variable, and its
    %mean value can be obtained by multplying
    %the density for the area.
    D = 1.2*R*sqrt(rand(1,B));  %Distance from the center of the different
    %buildings to the BS.
    gamma = 2*pi*rand(1,B); %Angle of the different buildings/obstacles
    %respective to the BS.
    
    x_b_c = [D.*cos(gamma); D.*sin(gamma); zeros(1,B)]; %Coordinates of the
    %center of the base
    %of the blockages.
    l = L_max*rand(1,B); %We obtain the length of every building.
    w = W_max*rand(1,B); %We obtain the width of every building.
    h = H_max*rand(1,B); %We obtain the height of every building.
    theta = Theta_max*rand(1,B); %We obtain the orientation of every
    %building.
    
    %We calculate some auxiliary angles so as to obtain the different
    %vertices.
    a = sqrt(w.^2+l.^2)/2;
    rho = atan(w./l);
    gamma_1 = theta-rho;
    gamma_2 = theta+rho;
    gamma_3 = pi+theta-rho;
    gamma_4 = pi+theta+rho;
    
    %With all this data, we are ready to define which are the 8 points
    %that define the building (its vertices).
    x_b_1 = [x_b_c(1,:)+a.*cos(gamma_1); x_b_c(2,:)+a.*sin(gamma_1); zeros(1,B)];
    x_b_2 = [x_b_c(1,:)+a.*cos(gamma_2); x_b_c(2,:)+a.*sin(gamma_2); zeros(1,B)];
    x_b_3 = [x_b_c(1,:)+a.*cos(gamma_3); x_b_c(2,:)+a.*sin(gamma_3); zeros(1,B)];
    x_b_4 = [x_b_c(1,:)+a.*cos(gamma_4); x_b_c(2,:)+a.*sin(gamma_4); zeros(1,B)];
    x_b_5 = [x_b_c(1,:)+a.*cos(gamma_1); x_b_c(2,:)+a.*sin(gamma_1); h];
    x_b_6 = [x_b_c(1,:)+a.*cos(gamma_2); x_b_c(2,:)+a.*sin(gamma_2); h];
    x_b_7 = [x_b_c(1,:)+a.*cos(gamma_3); x_b_c(2,:)+a.*sin(gamma_3); h];
    x_b_8 = [x_b_c(1,:)+a.*cos(gamma_4); x_b_c(2,:)+a.*sin(gamma_4); h];
    
    %For every simulation we check whether the BU link is blocked for every
    %position of the vector pos.
    BS_UE_blocked = zeros(1,num_pos);
    for s = 1:num_pos
        BS_UE_blocked(s) = link_blocked(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,x_B,x_U(:,s));
    end
    
    P_allKO_simulation = P_allKO_simulation+BS_UE_blocked;
    save('n_run.mat','n_run');
    n_run
    save('P_allKO_simulation.mat','P_allKO_simulation');
end
P_allKO_simulation = P_allKO_simulation/N_runs;
save('P_allKO_simulation.mat','P_allKO_simulation');


