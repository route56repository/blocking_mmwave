%clear;

%In this script, we try to compute which is the probability, in a
%sectorized case of N=3 sectors, that a user is blocked from the BS in a
%scenario where there are no RSs.

%The radius of the cell is 100 m.
R = 300;
d_d = 0.05*R;

%We create a vector to obtain the probability of blockage at given
%distances
d = 0:d_d:R;

%Base station's height is 40 m.
H_B = 40;

%UE's height of 1.5 m.
H_U = 1.5;
    
%Propagation model coefficients.
alpha = 4;
k = 1.35;

%We define parameters related to the blocking elements.
%The density of blocking elements is \lambda.
lambda = 2.2e-4;
% Length l is a r.v. such that: l~U[0,L_max]
L_max = 30;
E_L = L_max/2;
% Width w is a r.v. such that: w~U[0,W_max]
W_max = 30;
E_W = W_max/2;
% Height h is a r.v. such that: h~U[0,H_max]
H_max = 30;

% We compute the beta, p, eta and mu parameters
beta = 2*lambda*(E_L+E_W)/pi;
p = lambda*E_L*E_W;
eta = 1-1/(H_B-H_U)*((H_max^2-H_U^2)/(2*H_max)+H_B-H_max);
mu = 1-H_U/H_max;      
        
P_allKO_calculus_2_2 = 1-exp(-(eta*beta*d+mu*p));