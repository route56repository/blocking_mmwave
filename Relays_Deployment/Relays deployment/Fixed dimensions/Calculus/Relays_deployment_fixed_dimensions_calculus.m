clear;

warning('off','all');

% In this script, we try to simulate the case in which the buildings are
% modeled as rectangles with height. We also take into account the cases in
% which there is not a LOS link between the BS and the RS. Then we evaluate
% for every position of the relay possible (for every distance to the BS, r
% and for every height of the RS, h_R) which is the specific position in
% which the probability of no successful transmission is minimized.

% The radius of the cell is 300 m.
R = 300;

%We define some differentials that will be used to evaluate the
%granularities in which both the UE positions and the RS positions are
%analised.
d_d = 0.1*R;

% We have to obtain the azimut of every relay. There are going to be 3 RSs.
N = 3;
psi_n = ((1:N)-1)*2*pi/N;
r = 180;
h_R = 20;

%We get the position for the relays
x_n_1 = r*cos(psi_n); x_n_2 = r*sin(psi_n);

% Base station at the origin. Its height is 40 m.
x_B_1 = 0;
x_B_2 = 0;
H_B = 40;

% We define the height of the UE
H_U = 1.5;

%We define parameters related to the blocking elements.
%The density of blocking elements is \lambda.
lambda = 2.2e-4;
% Length l is fixed to 15 m.
l_0 = 15;
% Width w is fixed to 15 m.
w_0 = 15;
% Height h is a r.v. such that: h~U[0,H_max]
H_max = 30;
% Orientation \theta is a r.v. such that: l~U[0,\Theta_max].
Theta_max = pi;

% We compute the beta parameter
beta = 2*lambda*(l_0+w_0)/pi;
p = lambda*l_0*w_0;

n = 1; %Number of the section/relay that is being evaluated
psi_n_s = psi_n(n)-pi/N; %Angle where the sector (or area of integration) begins
psi_n_e = psi_n(n)+pi/N; %Angle where the sector (or area of integration) ends

%We define the azimuth and distances that we want to analyze
phi = pi/12;
d = 0:d_d:R;

P_allKO_calculus = zeros(1,length(d));

%For each value of d, the we obtain a value of the function
for s = 1:length(d)
    P_allKO_calculus(s) = P_allKO_d_phi_function(d(s),phi,x_n_1(n),x_n_2(n),l_0,w_0,H_max,Theta_max,H_B,H_U,h_R,x_B_1,x_B_2,beta,p,lambda);
    save('P_allKO_calculus.mat','P_allKO_calculus');
    distancia = d(s)
    save('distancia.mat','distancia')
end       
save('P_allKO_calculus.mat','P_allKO_calculus');


