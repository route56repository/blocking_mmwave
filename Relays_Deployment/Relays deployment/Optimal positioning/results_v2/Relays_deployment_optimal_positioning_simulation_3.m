clear;

%In this script, we try to simulate the case in which the buildings are
%modeled as rectangles with height. We take into account the cases in which
%there is not always a LOS link between the BS and the RS. Then we evaluate
%for every position of the relay possible (for every distance to the BS, r
%and for every height of the RS, h_R) which is the specific position in
%which the probability of no successful transmission is minimized.
%THE VECTORS THAT DEFINE THE COORDINATES ARE COLUMN VECTORS OF SIZE 3x1.

%The radius of the cell is 100 m.
R = 300;

%We define some differentials that will be used to evaluate the
%granularities in which both the UE positions and the RS positions are
%analised.
d_d = 0.02*R;
d_r = 50;
d_h_R = 10;

%We obtain the positions where the RSs are going to be evaluated.
r_min = 0;
r_max = R;

h_R_min = 10;
h_R_max = 20;

[r,h_R] = meshgrid(r_min:d_r:r_max,h_R_min:d_h_R:h_R_max);
[U,V] = size(r);

%We have to obtain the azimuth of every relay. There are going to be 3 RSs.
N = 3;
psi_n = ((1:N)-1)*2*pi/N;

%Base station at the origin. Its height is 40 m.
x_B_1 = 0;
x_B_2 = 0;
H_B = 40;
x_B = [x_B_1; x_B_2; H_B];

%We define the grid where the users are going to be located.
%[x_U_1, x_U_2] = meshgrid(-R-d_d:d_d:R+d_d, -R-d_d:d_d:R+d_d);
[x_U_1, x_U_2] = meshgrid(-R:d_d:R, -R:d_d:R);
pos = reshape(x_U_1+1i*x_U_2,1,[]);
%pos(abs(pos(:))>R+d_d/2) = [];
pos(abs(pos(:))>R) = [];
num_pos = length(pos);
H_U = 1.5;

x_U = [real(pos); imag(pos); H_U*ones(1,num_pos)];

%For each position possible of the user, we can calculate its distance d to
%the BS and its azimuth phi.
d = abs(pos);
phi = wrapTo2Pi(angle(pos));

%It has to be seen whether the user can be connected to the nearest relay,
%that is the relay of the sector in which it is found.
sector = zeros(1,num_pos);
for s = 1:num_pos
    [~,sector(s)] = min(abs(phi(s)-psi_n));
end

%So as to add sensitivity as a constraint, we have to take into account
%some parameters and convert them to linear.
P_T_B = 25; %Tx power of the BS in dBm.
    P_T_B = 10^(P_T_B/10)/1000;
    
P_T_R = 20; %Tx power of the RS in dBm.
    P_T_R = 10^(P_T_R/10)/1000;
    
S_R = -80 + 15; %Sensitivity of the RS in dBm + security margin for fading and atmospheric effects (dB).
    S_R = 10^(S_R/10)/1000;
    
S_U = -80 + 15; %Sensitivity of the UE in dBm + security margin for fading and atmospheric effects (dB).
    S_U = 10^(S_U/10)/1000;
    
G_B = 23; %BS's gain in dB.
    G_B = 10^(G_B/10);
    
G_R = 23; %RS's gain in dB.
    G_R = 10^(G_R/10);
    
G_U = 0; %UE's gain in dB.
    G_U = 10^(G_U/10);
    
%Propagation model coefficients.
alpha = 2.3;
k = 70.59;   % path-loss (in dB) at unitary distance
    k = 10^(k/10);

    
%After these parameters, we are going to define some other parameters that
%are key so as to implement the indicator functions that are set for the
%power restrictions.

%BS-RS restrictions.
r_min = 0;
r_max = min([R, ((P_T_B*G_B*G_R)/(k*S_R))^(1/alpha)]);

h_R_min = zeros(U,V)+((H_B-sqrt(((P_T_B*G_B*G_R)/(k*S_R))^(2/alpha)-r.^2))>=0).*(H_B-sqrt(((P_T_B*G_B*G_R)/(k*S_R))^(2/alpha)-r.^2));
h_R_max = H_B+sqrt(((P_T_B*G_B*G_R)/(k*S_R))^(2/alpha)-r.^2);

%BS-UE restriction.
d_min = 0;
d_max = min([R, sqrt(((P_T_B*G_B*G_U)/(k*S_U))^(2/alpha)-(H_B-H_U)^2)]);

%RS-UE restriction.
RU_min = 0;
RU_max = ((P_T_R*G_R*G_U)/(k*S_U))^(2/alpha)-(h_R-H_U).^2;

%Sensitivity indicator functions.
%If we want to add the sensitivity part, we must check the constraint for
%the BR link.
indicator_BR = r>=r_min & r<=r_max & h_R>=h_R_min & h_R<=h_R_max;

%We compute another indicative function.
indicator_BU = d>=d_min & d<=d_max & (H_B-H_U)^2<=((P_T_B*G_B*G_U)/(k*S_U))^(2/alpha);

%We define parameters related to the blocking elements.
%The density of blocking elements is \lambda.
lambda = 2.2e-4;
%Length l is a r.v. such that: l~U[0,L_max].
L_max = 30;
%Width w is a r.v. such that: w~U[0,W_max].
W_max = 30;
%Height h is a r.v. such that: h~U[0,H_max].
H_max = 30;
%Orientation \theta is a r.v. such that: \theta~U[0,\Theta_max].
Theta_max = pi;

N_runs = 100;

P_allKO_no_sensitivity_simulation = zeros(U,V);
P_allKO_sensitivity_simulation = zeros(U,V);

for n_run = 1:N_runs
    P_allKO_no_sensitivity_temporary = zeros(U,V);
    P_allKO_sensitivity_temporary = zeros(U,V);
    
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
                                 
    for u = 1:U
        for v = 1:V
            % We obtain the coordinates of all the RSs.
            x_n = [r(u,v)*cos(psi_n); r(u,v)*sin(psi_n); h_R(u,v)*ones(1,N)];
            
            %At the given position of the RS, we calculate the remaining
            %indicative function.
            indicator_RU = d.^2+r(u,v)^2-2*d*r(u,v).*cos(phi-psi_n(sector))>=RU_min & d.^2+r(u,v)^2-2*d*r(u,v).*cos(phi-psi_n(sector))<=RU_max(u,v);
            
            %Now, we can calculate how many RS are blocked, that is, in
            %which links from RS to BS have blockages.
            BS_RS_n_blocked = zeros(1,N);
            for n = 1:N
                BS_RS_n_blocked(n) = link_blocked(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,x_B,x_n(:,n));
            end
            
            %Once the relay is fixed in a specific position, we check if 
            %the R_n_U link is blocked for every position in the vector
            %pos.
            RS_n_UE_blocked = zeros(1,num_pos);
            for s = 1:num_pos
                RS_n_UE_blocked(s) = link_blocked(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,x_n(:,sector(s)),x_U(:,s));
            end
            
            UE_blocked_no_sensitivity = BS_UE_blocked & (BS_RS_n_blocked(sector) | RS_n_UE_blocked);
            UE_blocked_sensitivity = (~((~BS_UE_blocked) & indicator_BU)) & ((~((~BS_RS_n_blocked(sector)) & indicator_BR(u,v))) | (~((~RS_n_UE_blocked) & indicator_RU)));
            
            P_allKO_no_sensitivity_temporary(u,v) = sum(UE_blocked_no_sensitivity==1)/num_pos;
            P_allKO_sensitivity_temporary(u,v) = sum(UE_blocked_sensitivity==1)/num_pos;
            
            distance_relay_simulation = r(u,v);
            height_relay_simulation = h_R(u,v);
            save('distance_relay_simulation.mat','distance_relay_simulation');
            save('height_relay_simulation.mat','height_relay_simulation');
        end
    end
    
    P_allKO_no_sensitivity_simulation = P_allKO_no_sensitivity_simulation+P_allKO_no_sensitivity_temporary;
    P_allKO_sensitivity_simulation = P_allKO_sensitivity_simulation+P_allKO_sensitivity_temporary;
    save('n_run.mat','n_run');
    n_run
    save('P_allKO_no_sensitivity_simulation.mat','P_allKO_no_sensitivity_simulation');
    save('P_allKO_sensitivity_simulation.mat','P_allKO_sensitivity_simulation');
end

P_allKO_no_sensitivity_simulation = P_allKO_no_sensitivity_simulation/N_runs;
P_allKO_sensitivity_simulation = P_allKO_sensitivity_simulation/N_runs;
save('P_allKO_no_sensitivity_simulation.mat','P_allKO_no_sensitivity_simulation');
save('P_allKO_sensitivity_simulation.mat','P_allKO_sensitivity_simulation');
save('r.mat','r');