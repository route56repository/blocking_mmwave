clear;

%BS coordinates
x_B_1 = 0;
x_B_2 = 0;
H_B = 40;

%UE coordinates
d = 180;
phi = pi/6;
x_U_1 = d*cos(phi);
x_U_2 = d*sin(phi);
H_U = 1.5;

%Dimensions of the blocking
l = 145;
w = 15;
h = 70;
theta = 0;

%We define some some angles and dimensions that will be used throughout the
%function
a = sqrt(w^2+l^2)/2;
rho = atan(w/l);
gamma_1 = theta-rho;
gamma_2 = theta+rho;
gamma_3 = pi+theta-rho;
gamma_4 = pi+theta+rho;

%The positions of the rectangle related with the UE
x_U_1_1 = 0+(x_U_1+a.*cos(gamma_1))*(h>=H_U);
x_U_2_1 = 0+(x_U_2+a.*sin(gamma_1))*(h>=H_U);

x_U_1_2 = 0+(x_U_1+a.*cos(gamma_2))*(h>=H_U);
x_U_2_2 = 0+(x_U_2+a.*sin(gamma_2))*(h>=H_U);

x_U_1_3 = 0+(x_U_1+a.*cos(gamma_3))*(h>=H_U);
x_U_2_3 = 0+(x_U_2+a.*sin(gamma_3))*(h>=H_U);

x_U_1_4 = 0+(x_U_1+a.*cos(gamma_4))*(h>=H_U);
x_U_2_4 = 0+(x_U_2+a.*sin(gamma_4))*(h>=H_U);

%The positions of the rectangle related with the BS in the BU link
if (x_U_1+1i*x_U_2)==(x_B_1+1i*x_B_2) 
    x_B_BU_1_1 = 0+(x_B_1+a*cos(gamma_1))*(h>=H_B);
    x_B_BU_2_1 = 0+(x_B_2+a*sin(gamma_1))*(h>=H_B);
    
    x_B_BU_1_2 = 0+(x_B_1+a*cos(gamma_2))*(h>=H_B);
    x_B_BU_2_2 = 0+(x_B_2+a*sin(gamma_2))*(h>=H_B);
    
    x_B_BU_1_3 = 0+(x_B_1+a*cos(gamma_3))*(h>=H_B);
    x_B_BU_2_3 = 0+(x_B_2+a*sin(gamma_3))*(h>=H_B);
    
    x_B_BU_1_4 = 0+(x_B_1+a*cos(gamma_4))*(h>=H_B);
    x_B_BU_2_4 = 0+(x_B_2+a*sin(gamma_4))*(h>=H_B);
else
    x_B_BU_1_1 = 0+real(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_1*(h>=H_B)+a*cos(gamma_1)*(h>=H_U);
    x_B_BU_2_1 = 0+imag(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_2*(h>=H_B)+a*sin(gamma_1)*(h>=H_U);
    
    x_B_BU_1_2 = 0+real(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_1*(h>=H_B)+a*cos(gamma_2)*(h>=H_U);
    x_B_BU_2_2 = 0+imag(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_2*(h>=H_B)+a*sin(gamma_2)*(h>=H_U);
    
    x_B_BU_1_3 = 0+real(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_1*(h>=H_B)+a*cos(gamma_3)*(h>=H_U);
    x_B_BU_2_3 = 0+imag(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_2*(h>=H_B)+a*sin(gamma_3)*(h>=H_U);
    
    x_B_BU_1_4 = 0+real(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_1*(h>=H_B)+a*cos(gamma_4)*(h>=H_U);
    x_B_BU_2_4 = 0+imag(x_B_1+1i*x_B_2+((x_U_1+1i*x_U_2-(x_B_1+1i*x_B_2))*(H_B-h)/(H_B-H_U)))*(h>=H_U & h<H_B)+x_B_2*(h>=H_B)+a*sin(gamma_4)*(h>=H_U);
end

%First of all, we obtain the S_BU_l_w_h_theta polygon:
S_BU_l_w_h_theta_1 = polyshape([x_B_BU_1_1 x_B_BU_1_2 x_B_BU_1_3 x_B_BU_1_4],[x_B_BU_2_1 x_B_BU_2_2 x_B_BU_2_3 x_B_BU_2_4]);
S_BU_l_w_h_theta_2 = polyshape([x_U_1_1 x_U_1_2 x_U_1_3 x_U_1_4],[x_U_2_1 x_U_2_2 x_U_2_3 x_U_2_4]);
S_BU_l_w_h_theta_3 = polyshape([x_B_BU_1_1 x_B_BU_1_3 x_U_1_3 x_U_1_1],[x_B_BU_2_1 x_B_BU_2_3 x_U_2_3 x_U_2_1]);
S_BU_l_w_h_theta_4 = polyshape([x_B_BU_1_2 x_B_BU_1_4 x_U_1_4 x_U_1_2],[x_B_BU_2_2 x_B_BU_2_4 x_U_2_4 x_U_2_2]);
        
S_BU_l_w_h_theta = union(union(S_BU_l_w_h_theta_1,S_BU_l_w_h_theta_2),union(S_BU_l_w_h_theta_3,S_BU_l_w_h_theta_4));

plot(S_BU_l_w_h_theta)
hold on;
scatter3(x_B_1, x_B_2, H_B,'o')
scatter3(x_U_1, x_U_2, H_U,'x','black')
hold off;

A_analytic = 0+((h-H_U)/(H_B-H_U)*d*(l*sin(theta)+w*abs(cos(theta)))+w*l)*(h>=H_U & h<H_B)+(d*(l*sin(theta)+w*abs(cos(theta)))+w*l)*(h>=H_B)
A_numeric = area(S_BU_l_w_h_theta)
        