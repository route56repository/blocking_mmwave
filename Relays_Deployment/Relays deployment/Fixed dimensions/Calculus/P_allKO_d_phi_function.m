function P_allKO_d_phi = P_allKO_d_phi_function(d,phi,x_n_1,x_n_2,l,w,H_max,Theta_max,H_B,H_U,h_R,x_B_1,x_B_2,beta,p,lambda)
%UNTITLED3 Summary of this function goes here

eta_BU = 1-1/(H_B-H_U)*((H_max^2-H_U^2)/(2*H_max)+H_B-H_max);
mu_BU = 1-H_U/H_max;

x_U_1 = d*cos(phi);     x_U_2 = d*sin(phi);

%For each position of the user, we calculate
E_K_BU = eta_BU*beta*d+mu_BU*p;
E_K_BR_RU = integral2(@(h,theta)E_K_BR_RU_l_w_h_theta_function(l,w,h,theta,x_U_1,x_U_2,x_n_1,x_n_2,H_max,Theta_max,H_B,H_U,h_R,x_B_1,x_B_2,lambda),0,H_max,0,Theta_max,'Method','tiled','AbsTol', 5e-3,'RelTol',5e-3);
E_K_BR_RU_BU = integral2(@(h,theta)E_K_BR_RU_BU_l_w_h_theta_function(l,w,h,theta,x_U_1,x_U_2,x_n_1,x_n_2,H_max,Theta_max,H_B,H_U,h_R,x_B_1,x_B_2,lambda),0,H_max,0,Theta_max,'Method','tiled','AbsTol', 5e-3,'RelTol',5e-3);

P_allKO_d_phi = 1-exp(-E_K_BU)-exp(-E_K_BR_RU)+exp(-E_K_BR_RU_BU);

end

