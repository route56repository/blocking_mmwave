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
d_d = 0.1*R;

%We have to obtain the azimuth of every relay. There are going to be 3 RSs.
N = 3;
psi_n = ((1:N)-1)*2*pi/N;
r = 180;
h_R = 20;

% We obtain the coordinates of all the RSs.
x_n = [r*cos(psi_n); r*sin(psi_n); h_R*ones(1,N)];

%Base station at the origin. Its height is 40 m.
x_B_1 = 0;
x_B_2 = 0;
H_B = 40;
x_B = [x_B_1; x_B_2; H_B];

%We define the grid where the users are going to be located.
d = 0:d_d:R;
phi = pi/4;
pos = d*cos(phi)+1i*d*sin(phi);
num_pos = length(pos);
H_U = 1.5;

x_U = [real(pos); imag(pos); H_U*ones(1,num_pos)];

%It has to be seen whether the user can be connected to the nearest relay,
%that is the relay of the sector in which it is found.
sector = zeros(1,num_pos);
for s = 1:num_pos
    [~,sector(s)] = min(abs(phi-psi_n));
end

%We define parameters related to the blocking elements.
%The density of blocking elements is \lambda.
lambda = 2.2e-4;
% Length l is fixed to 15 m.
l_0 = 15;
% Width w is fixed to 15 m.
w_0 = 15;
%Height h is a r.v. such that: h~U[0,H_max].
H_max = 30;
%Orientation \theta is a r.v. such that: \theta~U[0,\Theta_max].
Theta_max = pi;

N_runs = 20000;

P_allKO_simulation = zeros(1,num_pos);

for n_run = 1:N_runs
    P_allKO_temporary = zeros(1,num_pos);
    
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
    l = l_0*ones(1,B); %We obtain the length of every building.
    w = w_0*ones(1,B); %We obtain the width of every building.
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
    
    UE_blocked = BS_UE_blocked & (BS_RS_n_blocked(sector) | RS_n_UE_blocked);
    
    P_allKO_simulation = P_allKO_simulation+UE_blocked;
    n_run
    save('n_run.mat','n_run');
    save('P_allKO_simulation.mat','P_allKO_simulation');
end

P_allKO_simulation = P_allKO_simulation/N_runs;
save('P_allKO_simulation.mat','P_allKO_simulation');