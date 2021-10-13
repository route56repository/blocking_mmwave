function [intersect] = segment_intersects_convex_hull(x_b_1,x_b_2,x_b_3,x_b_4,x_b_5,x_b_6,x_b_7,x_b_8,P1,P2)
%This function creates the appropriate matrices to feeed the linprog
%function and obtains the result from such function.
%Basically, the variable intersect equals 1 when the segment from P1 to P2
%intersects the convex hull formed by the points x_b_1, x_b_2, x_b_3,
%x_b_4, x_b_5, x_b_6, x_b_7 and x_b_8. Otherwise, it is 0.
beq = [zeros(3,1)
       1
       -P1];
   
Aeq = [x_b_1 x_b_2 x_b_3 x_b_4 x_b_5 x_b_6 x_b_7 x_b_8 -eye(3) zeros(3,1)
       ones(1,8) zeros(1,4)
       zeros(3,8) -eye(3) P2-P1];

lb = [zeros(8,1)
      -Inf*ones(3,1)
      0];
  
ub = [Inf*ones(11,1)
      1];
  
f=ones(12,1);

options = optimoptions('linprog','Display','off');

[~,~,exitflag,~] = linprog(f,[],[],Aeq,beq,lb,ub,options);
    
intersect = exitflag == 1;

end

