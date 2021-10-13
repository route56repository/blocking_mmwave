close all
clear all
clc

% We load the files with the results
load('shadowing_intervals_analytic.mat');
load('shadowing_intervals_empiric.mat');

% We load the parameters
[r, scaling_factor, gamma, d, A] = load_parameters();

line_emp{1}='b--';marker_emp{1}='bs';legend_emp{1}='bs--';
line_emp{2}='r--';marker_emp{2}='r*';legend_emp{2}='r*--';
line_emp{3}='k--';marker_emp{3}='ko';legend_emp{3}='ko--';
line_emp{4}='m--';marker_emp{4}='mx';legend_emp{4}='mx--';
line_emp{5}='g--';marker_emp{5}='g^';legend_emp{5}='g^--';
line_emp{6}='c--';marker_emp{6}='c<';legend_emp{6}='c<--';

line_ana{1}='b';marker_ana{1}='b+';legend_ana{1}='b+-';
line_ana{2}='r';marker_ana{2}='rv';legend_ana{2}='rv-';
line_ana{3}='k';marker_ana{3}='kd';legend_ana{3}='kd-';
line_ana{4}='m';marker_ana{4}='m>';legend_ana{4}='m>-';
line_ana{5}='g';marker_ana{5}='gp';legend_ana{5}='gp-';
line_ana{6}='c';marker_ana{6}='ch';legend_ana{6}='ch-';

lw=2;  %Fat lines (such as line_width=2) are nice in the pdf. To see how good are bounds, set lw to 1
fs=12; %Big font size (such as font size =12) are nice in the pdf. It can be changed to 10, to see figures clearly when working
ms=8; %%Big marker size (such as marker size =8) are nice in the pdf.
md=3; %Distance between markers

% To represent the axis, first we can create an array with the values where
% we want to draw the markers for both analytical and empirical curves
d_axis = 100;

% We create one figure for each r in the calculus. Then, for a given r, the
% empiric result is represented together with the different analytic CDFs
% for different values of lambda
p = 1;
q = 1;
figure(1)
for n = 1:length(scaling_factor)
    % We predefine the size of the legend container
    Legend = cell(2*length(r),1);
    subplot(length(scaling_factor),1,n)
    hold on;
    for m = 1:length(r)
        % We define at what positions the markers should be placed
        pos_markers_z_emp = (m-1)*d_axis/length(r):d_axis:max(z);
        pos_markers_z_ana = d_axis/2+(m-1)*d_axis/length(r):d_axis:max(z);
        
        [~,ind_z_emp] = min(abs(z_emp{m,n}-pos_markers_z_emp));
        [~,ind_z_ana] = min(abs(z'-pos_markers_z_ana));
        
        plot(z_emp{m,n},F_Z_emp{m,n},line_emp{q},'LineWidth',lw,'HandleVisibility','off'); plot(z_emp{m,n}(ind_z_emp),F_Z_emp{m,n}(ind_z_emp),marker_emp{q},'LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,legend_emp{q},'LineWidth',lw,'MarkerSize',ms); 
        Legend{p} = ['Empiric, r=',num2str(r(m)),' m'];
        p = p + 1;
        
        plot(z,F_Z_ub{m,n},line_ana{q},'LineWidth',lw,'HandleVisibility','off'); plot(z(ind_z_ana),F_Z_ub{m,n}(ind_z_ana),marker_ana{q},'LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,legend_ana{q},'LineWidth',lw,'MarkerSize',ms); 
        Legend{p} = ['Analytic, \lambda_{lines}=',num2str(gamma(m)),'\lambda_{buildings}, r=', num2str(r(m)),' m'];
        p = p + 1;
        q = q + 1;
    end
    legend(Legend,'Location','southeast');
    title(['CDFs for different distances r, for \lambda_{buildings}=',num2str(NumberBuildings/scaling_factor(n)/A,'%.2e'),' m^{-2}']);
    
    set(gca,'fontsize',10)
    xlim([0 300]);
    ylim([0 1]);
    hold off
    grid on
    grid minor 
    p = 1;
    q = 1;
    Legend = [];
end
xlabel('Length of the LOS intervals, z (m)')
