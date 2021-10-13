function [blocked] = link_blocked(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,P1,P2)
%This function computes if the link represented by the segment given by P1
%and P2 is blocked by any of the B buildings defined by the matrices x_b_c
%(a 3xB matrix in which the coordinates of the centre of the base of each 
%building are stored) and x_b_i with i=1,...,8 (3xB matrices that contain 
%the 8 vertices that define each building).

%First of all, most buildings are discarded, since they will surely not
%block the link. After this comparison, we obtain the matrices x_b_i with
%i=1,...,8 that correspond to the buildings that may block the link. We 
%also obtain B, which is the number of these buildings.
[x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,B] = discard_buildings(x_b_c,x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,P1,P2);

blocked = 0;
b = 1;
while blocked == 0 && b <= B
    blocked = segment_intersects_convex_hull(x_b_1(:,b), x_b_2(:,b), x_b_3(:,b), x_b_4(:,b), x_b_5(:,b), x_b_6(:,b), x_b_7(:,b), x_b_8(:,b), P1, P2);
    b = b + 1;
end
end

