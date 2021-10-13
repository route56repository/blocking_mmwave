function E_K_BR_RU_l_w_h_theta = E_K_BR_RU_l_w_h_theta_function(l,w,h,theta,x_U_1,x_U_2,x_n_1,x_n_2,H_max,Theta_max,H_B,H_U,h_R,x_B_1,x_B_2,lambda)
%This function calculates the mean number of blockages both in the BR
%and/or in the RU link.
[P,Q]=size(h);
E_K_BR_RU_l_w_h_theta=zeros(P,Q);

l = l*ones(P,Q);
w = w*ones(P,Q);

%We define some some angles and dimensions that will be used throughout the
%functiºon
a = sqrt(w.^2+l.^2)/2;
rho = atan(w./l);
gamma_1 = theta-rho;
gamma_2 = theta+rho;
gamma_3 = pi+theta-rho;
gamma_4 = pi+theta+rho;

%The positions of the rectangle related with the UE
x_U_1_1 = zeros(P,Q)+(x_U_1+a.*cos(gamma_1)).*(h>=H_U);
x_U_2_1 = zeros(P,Q)+(x_U_2+a.*sin(gamma_1)).*(h>=H_U);

x_U_1_2 = zeros(P,Q)+(x_U_1+a.*cos(gamma_2)).*(h>=H_U);
x_U_2_2 = zeros(P,Q)+(x_U_2+a.*sin(gamma_2)).*(h>=H_U);

x_U_1_3 = zeros(P,Q)+(x_U_1+a.*cos(gamma_3)).*(h>=H_U);
x_U_2_3 = zeros(P,Q)+(x_U_2+a.*sin(gamma_3)).*(h>=H_U);

x_U_1_4 = zeros(P,Q)+(x_U_1+a.*cos(gamma_4)).*(h>=H_U);
x_U_2_4 = zeros(P,Q)+(x_U_2+a.*sin(gamma_4)).*(h>=H_U);

%The positions of the rectangle related with the RS in the BR link
x_n_BR_1_1 = zeros(P,Q)+(x_n_1+a.*cos(gamma_1)).*(h>=h_R);
x_n_BR_2_1 = zeros(P,Q)+(x_n_2+a.*sin(gamma_1)).*(h>=h_R);

x_n_BR_1_2 = zeros(P,Q)+(x_n_1+a.*cos(gamma_2)).*(h>=h_R);
x_n_BR_2_2 = zeros(P,Q)+(x_n_2+a.*sin(gamma_2)).*(h>=h_R);

x_n_BR_1_3 = zeros(P,Q)+(x_n_1+a.*cos(gamma_3)).*(h>=h_R);
x_n_BR_2_3 = zeros(P,Q)+(x_n_2+a.*sin(gamma_3)).*(h>=h_R);

x_n_BR_1_4 = zeros(P,Q)+(x_n_1+a.*cos(gamma_4)).*(h>=h_R);
x_n_BR_2_4 = zeros(P,Q)+(x_n_2+a.*sin(gamma_4)).*(h>=h_R);

%The positions of the rectangle related with the BS in the BR link
x_B_BR_1_1 = zeros(P,Q)+real(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_1.*(h>=H_B)+a.*cos(gamma_1).*(h>=h_R);
x_B_BR_2_1 = zeros(P,Q)+imag(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_2.*(h>=H_B)+a.*sin(gamma_1).*(h>=h_R);

x_B_BR_1_2 = zeros(P,Q)+real(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_1.*(h>=H_B)+a.*cos(gamma_2).*(h>=h_R);
x_B_BR_2_2 = zeros(P,Q)+imag(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_2.*(h>=H_B)+a.*sin(gamma_2).*(h>=h_R);

x_B_BR_1_3 = zeros(P,Q)+real(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_1.*(h>=H_B)+a.*cos(gamma_3).*(h>=h_R);
x_B_BR_2_3 = zeros(P,Q)+imag(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_2.*(h>=H_B)+a.*sin(gamma_3).*(h>=h_R);

x_B_BR_1_4 = zeros(P,Q)+real(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_1.*(h>=H_B)+a.*cos(gamma_4).*(h>=h_R);
x_B_BR_2_4 = zeros(P,Q)+imag(x_B_1+1i*x_B_2+((x_n_1+1i*x_n_2-(x_B_1+1i*x_B_2)).*(H_B-h)./(H_B-h_R))).*(h>=h_R & h<H_B)+x_B_2.*(h>=H_B)+a.*sin(gamma_4).*(h>=h_R);

%The positions of the rectangle related with the RS in the RU link
x_n_RU_1_1 = zeros(P,Q)+real(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_1.*(h>=h_R)+a.*cos(gamma_1).*(h>=H_U);
x_n_RU_2_1 = zeros(P,Q)+imag(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_2.*(h>=h_R)+a.*sin(gamma_1).*(h>=H_U);

x_n_RU_1_2 = zeros(P,Q)+real(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_1.*(h>=h_R)+a.*cos(gamma_2).*(h>=H_U);
x_n_RU_2_2 = zeros(P,Q)+imag(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_2.*(h>=h_R)+a.*sin(gamma_2).*(h>=H_U);

x_n_RU_1_3 = zeros(P,Q)+real(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_1.*(h>=h_R)+a.*cos(gamma_3).*(h>=H_U);
x_n_RU_2_3 = zeros(P,Q)+imag(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_2.*(h>=h_R)+a.*sin(gamma_3).*(h>=H_U);

x_n_RU_1_4 = zeros(P,Q)+real(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_1.*(h>=h_R)+a.*cos(gamma_4).*(h>=H_U);
x_n_RU_2_4 = zeros(P,Q)+imag(x_n_1+1i*x_n_2+((x_U_1+1i*x_U_2-(x_n_1+1i*x_n_2)).*(h_R-h)./(h_R-H_U))).*(h>=H_U & h<h_R)+x_n_2.*(h>=h_R)+a.*sin(gamma_4).*(h>=H_U);

for p=1:P
    for q=1:Q
        %First of all, we obtain the S_BR_l_w_h_theta polygon:
        S_BR_l_w_h_theta_1 = polyshape([x_B_BR_1_1(p,q) x_B_BR_1_2(p,q) x_B_BR_1_3(p,q) x_B_BR_1_4(p,q)],[x_B_BR_2_1(p,q) x_B_BR_2_2(p,q) x_B_BR_2_3(p,q) x_B_BR_2_4(p,q)]);
        S_BR_l_w_h_theta_2 = polyshape([x_n_BR_1_1(p,q) x_n_BR_1_2(p,q) x_n_BR_1_3(p,q) x_n_BR_1_4(p,q)],[x_n_BR_2_1(p,q) x_n_BR_2_2(p,q) x_n_BR_2_3(p,q) x_n_BR_2_4(p,q)]);
        S_BR_l_w_h_theta_3 = polyshape([x_B_BR_1_1(p,q) x_B_BR_1_3(p,q) x_n_BR_1_3(p,q) x_n_BR_1_1(p,q)],[x_B_BR_2_1(p,q) x_B_BR_2_3(p,q) x_n_BR_2_3(p,q) x_n_BR_2_1(p,q)]);
        S_BR_l_w_h_theta_4 = polyshape([x_B_BR_1_2(p,q) x_B_BR_1_4(p,q) x_n_BR_1_4(p,q) x_n_BR_1_2(p,q)],[x_B_BR_2_2(p,q) x_B_BR_2_4(p,q) x_n_BR_2_4(p,q) x_n_BR_2_2(p,q)]);
        
        S_BR_l_w_h_theta = union(union(S_BR_l_w_h_theta_1,S_BR_l_w_h_theta_2),union(S_BR_l_w_h_theta_3,S_BR_l_w_h_theta_4));
        
        %Now, we do the same for the S_RU_l_w_h_theta polygon
        S_RU_l_w_h_theta_1 = polyshape([x_n_RU_1_1(p,q) x_n_RU_1_2(p,q) x_n_RU_1_3(p,q) x_n_RU_1_4(p,q)],[x_n_RU_2_1(p,q) x_n_RU_2_2(p,q) x_n_RU_2_3(p,q) x_n_RU_2_4(p,q)]);
        S_RU_l_w_h_theta_2 = polyshape([x_U_1_1(p,q) x_U_1_2(p,q) x_U_1_3(p,q) x_U_1_4(p,q)],[x_U_2_1(p,q) x_U_2_2(p,q) x_U_2_3(p,q) x_U_2_4(p,q)]);
        S_RU_l_w_h_theta_3 = polyshape([x_n_RU_1_1(p,q) x_n_RU_1_3(p,q) x_U_1_3(p,q) x_U_1_1(p,q)],[x_n_RU_2_1(p,q) x_n_RU_2_3(p,q) x_U_2_3(p,q) x_U_2_1(p,q)]);
        S_RU_l_w_h_theta_4 = polyshape([x_n_RU_1_2(p,q) x_n_RU_1_4(p,q) x_U_1_4(p,q) x_U_1_2(p,q)],[x_n_RU_2_2(p,q) x_n_RU_2_4(p,q) x_U_2_4(p,q) x_U_2_2(p,q)]);
        
        S_RU_l_w_h_theta = union(union(S_RU_l_w_h_theta_1,S_RU_l_w_h_theta_2),union(S_RU_l_w_h_theta_3,S_RU_l_w_h_theta_4));
        
        E_K_BR_RU_l_w_h_theta(p,q) = lambda*area(union(S_BR_l_w_h_theta,S_RU_l_w_h_theta))*1/H_max*1/Theta_max;   
    end
end

end

