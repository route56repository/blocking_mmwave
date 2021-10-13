function [x_min_light, x_max_light, x_min_shadow, x_max_shadow] = obtain_lights_shadows(x_coord, y_coord, x_B, y_B, r, d)
x_coord_aux = x_coord;
y_coord_aux = y_coord;

% A segment that is totally out of the BS-trajectory area is discarded
x_coord_aux(sum(y_coord>(y_B+r),2)==2 | sum(y_coord<y_B,2)==2,:)=[];
y_coord_aux(sum(y_coord>(y_B+r),2)==2 | sum(y_coord<y_B,2)==2,:)=[];

% If one of its vertices is out of the BS-trajectory area, this vertix is
% moved to the limit
m = repmat((y_coord_aux(:,2)-y_coord_aux(:,1))./(x_coord_aux(:,2)-x_coord_aux(:,1)),1,2);
n = repmat(y_coord_aux(:,1)-m.*x_coord_aux(:,1),1,2);

x_coord_aux(y_coord_aux>(y_B+r)) = (r+y_B-n(y_coord_aux>(y_B+r)))./m(y_coord_aux>(y_B+r));
y_coord_aux(y_coord_aux>(y_B+r)) = r + y_B;

x_coord_aux(y_coord_aux<y_B) = (y_B-n(y_coord_aux<y_B))./m(y_coord_aux<y_B);
y_coord_aux(y_coord_aux<y_B) = y_B;

% We obtain the angles of the segments with respect the BS
alpha_1 = angle(exp(1i*(angle(x_coord_aux(:,1)+1i*y_coord_aux(:,1)-(x_B+1i*y_B)))));
alpha_2 = angle(exp(1i*(angle(x_coord_aux(:,2)+1i*y_coord_aux(:,2)-(x_B+1i*y_B)))));

%We obtain the values for the limits of the shadow for each blockage
x_lim(:,1) = r./tan(alpha_1)+x_B;
x_lim(:,2) = r./tan(alpha_2)+x_B;

% We sort the limits of the shadows, being the value on the left the
% minimum and the maximum on the right
x_lim = sort(x_lim,2);

% Given that the shadows of the elements can be found in any order, they
% are sorted upon the value of the left limit
[x_lim(:,1), I] = sort(x_lim(:,1),'ascend');
x_lim(:,2) = x_lim(I,2);

% % ONLY FOR DEBUGGING PURPOSES
% % We also sort the auxiliary coordinates according to the ordering of the
% % phases
% I = repmat(I,1,2);
% x_coord_aux = x_coord_aux(I);
% y_coord_aux = y_coord_aux(I);

%We take into account when alpha = pi/2
x_lim(isnan(x_lim(:,1)),1) = 0;
x_lim(isnan(x_lim(:,2)),1) = 0;

[N_segments_aux,~] = size(x_lim);

x_min_shadow = [];
x_max_shadow = [];
x_min_light = [];
x_max_light = [];

%Now, the overlapping intervals need to be found
for n_segment=1:N_segments_aux
    %We take the first interval
    if n_segment == 1
        n = 1;
        x_min_shadow(n) = x_lim(n_segment,1);
        x_max_shadow(n) = x_lim(n_segment,2);
    else
        if x_lim(n_segment,1) <= x_max_shadow(n) && x_lim(n_segment,2) >= x_max_shadow(n)
            x_max_shadow(n) = x_lim(n_segment,2);
        elseif x_lim(n_segment,1) > x_max_shadow(n)
            n = n+1;
            x_min_shadow(n) = x_lim(n_segment,1);
            x_max_shadow(n) = x_lim(n_segment,2);
        end
    end
end

%From the shadows, we obtain the lights
x_min_light = x_max_shadow;
x_max_light = x_min_shadow;

%We eliminate the first and last values of the x_max_light and
%x_min_light
if ~isempty(x_max_light)
    x_max_light(1) = [];
    x_min_light(end) = [];
end

%We remove the shadows which are out of the cell (on the left side)
x_max_shadow(x_min_shadow < -d/2) = [];
x_min_shadow(x_min_shadow < -d/2) = [];

%We remove the shadows which are out of the cell (on the right side)
x_min_shadow(x_max_shadow > d/2) = [];
x_max_shadow(x_max_shadow > d/2) = [];

%We remove the lights which are out of the cell (on the left side)
x_max_light(x_min_light < -d/2) = [];
x_min_light(x_min_light < -d/2) = [];

%We remove the lights which are out of the cell (on the right side)
x_min_light(x_max_light > d/2) = [];
x_max_light(x_max_light > d/2) = [];

end

