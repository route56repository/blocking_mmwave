function [x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,B] = discard_buildings(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,P1,P2)
%This function discards all the buildings defined by the matrices x_b_c
%(a 3xB matrix in which the coordinates of the centre of the base of each 
%building are stored) and x_b_i with i=1,...,8 (3xB matrices that contain 
%the 8 vertices that define each building) that will surely not block the
%link given by P2 and P1. It returns the matrices x_b_i with i=1,...,8 that
%correspond to the buildings that may block the link and B, which is the
%number of these buildings.

%First, we take P1 as the origin and calculate the angles that the centers
%of the base of each building taking P1 as the origin form with respect to
%the horizontal axis.
alpha_c = angle(x_b_c(1,:)+1i*x_b_c(2,:)-(P1(1,:)+1i*P1(2,:)));

%Then, we calculate the angles that the vertices of the base of each
%building taking P1 as the origin form with respect the anles alpha_c.
alpha_1 = angle(exp(1i*(angle(x_b_1(1,:)+1i*x_b_1(2,:)-(P1(1,:)+1i*P1(2,:)))-alpha_c)));
alpha_2 = angle(exp(1i*(angle(x_b_2(1,:)+1i*x_b_2(2,:)-(P1(1,:)+1i*P1(2,:)))-alpha_c)));
alpha_3 = angle(exp(1i*(angle(x_b_3(1,:)+1i*x_b_3(2,:)-(P1(1,:)+1i*P1(2,:)))-alpha_c)));
alpha_4 = angle(exp(1i*(angle(x_b_4(1,:)+1i*x_b_4(2,:)-(P1(1,:)+1i*P1(2,:)))-alpha_c)));

%Then, we calculate the angle that P2 taking P1 as the origin forms with
%respect the anles alpha_c.
alpha_P2 = angle(exp(1i*(angle(P2(1,:)+1i*P2(2,:)-(P1(1,:)+1i*P1(2,:)))-alpha_c)));

%Now, we obtain the maximum and the minimum values of the angles alpha_i
%with i=1,...,4, that is, we obtain the delimiting angles of the possible
%occlusion sector that each building is associated with.
alpha_max = max([alpha_1; alpha_2; alpha_3; alpha_4],[],1);
alpha_min = min([alpha_1; alpha_2; alpha_3; alpha_4],[],1);

%We erase from the matrices that contain both the vertices of the base and
%the center of the buildings that will surely not block the link. 
x_b_1(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_2(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_3(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_4(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_5(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_6(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_7(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_8(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];
x_b_c(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P2 > alpha_max | alpha_P2 < alpha_min ))=[];


%Now, we do the same thing but taking P2 as the reference
%First, we calculate the angles that the centers of the base of each
%building taking P2 as the origin form with respect to the horizontal axis.
alpha_c = angle(x_b_c(1,:)+1i*x_b_c(2,:)-(P2(1,:)+1i*P2(2,:)));

%Then, we calculate the angles that the vertixes of the base of each
%building taking P2 as the origin form with respect the anles alpha_c.
alpha_1 = angle(exp(1i*(angle(x_b_1(1,:)+1i*x_b_1(2,:)-(P2(1,:)+1i*P2(2,:)))-alpha_c)));
alpha_2 = angle(exp(1i*(angle(x_b_2(1,:)+1i*x_b_2(2,:)-(P2(1,:)+1i*P2(2,:)))-alpha_c)));
alpha_3 = angle(exp(1i*(angle(x_b_3(1,:)+1i*x_b_3(2,:)-(P2(1,:)+1i*P2(2,:)))-alpha_c)));
alpha_4 = angle(exp(1i*(angle(x_b_4(1,:)+1i*x_b_4(2,:)-(P2(1,:)+1i*P2(2,:)))-alpha_c)));

%Then, we calculate the angle that P1 taking P2 as the origin forms with
%respect the anles alpha_c.
alpha_P1 = angle(exp(1i*(angle(P1(1,:)+1i*P1(2,:)-(P2(1,:)+1i*P2(2,:)))-alpha_c)));

%Now, we obtain the maximum and the minimum values of the angles alpha_i
%with i=1,...,4, that is, we obtain the delimiting angles of the possible
%occlusion sector that each building is associated with.
alpha_max = max([alpha_1; alpha_2; alpha_3; alpha_4],[],1);
alpha_min = min([alpha_1; alpha_2; alpha_3; alpha_4],[],1);

%We erase from the matrices that contain both the vertices of the base and
%the center of the buildings that will surely not block the link.
x_b_1(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_2(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_3(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_4(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_5(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_6(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_7(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_8(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];
x_b_c(:,(alpha_max <= pi/2 & alpha_min >= -pi/2) & (alpha_P1 > alpha_max | alpha_P1 < alpha_min ))=[];

%Finally, we obtain the final number of buildings that may block the link.
[~,B]=size(x_b_c);
end