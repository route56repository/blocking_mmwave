%Load the results
load('P_allKO_no_sensitivity_simulation.mat','P_allKO_no_sensitivity_simulation');
load('P_allKO_sensitivity_simulation.mat','P_allKO_sensitivity_simulation');
load('r.mat','r');

[~,M] = size(r);

lin{1}='b';marca{1}='bs';leg{1}='bs-';lin{2}='r';marca{2}='r*';leg{2}='r*-';lin{3}='k';marca{3}='ko';leg{3}='ko-';lin{4}='m';marca{4}='mx';leg{4}='mx-';lin{5}='g';marca{5}='g^';leg{5}='g^-';
linea{1}='b--';linea{2}='r--';linea{3}='k--';linea{4}='m--';linea{5}='g--';linea{6}='c-.';
lw=2;  %Fat lines (such as line_width=2) are nice in the pdf. To see how good are bounds, set lw to 1
fs=12; %Big font size (such as font size =12) are nice in the pdf. It can be changed to 10, to see figures clearly when working
ms=8; %%Big marker size (such as marker size =8) are nice in the pdf.
md=3; %Distance between markers

figure(1)
hold on
plot(r(1,:),P_allKO_no_sensitivity_simulation(1,1)*100*ones(M,1),linea{5},'LineWidth',lw)
plot(r(1,:),P_allKO_no_sensitivity_simulation(1,:)*100,leg{1},'LineWidth',lw,'MarkerSize',ms)
plot(r(1,:),P_allKO_no_sensitivity_simulation(2,:)*100,leg{2},'LineWidth',lw,'MarkerSize',ms)
plot(r(1,:),P_allKO_sensitivity_simulation(1,1)*100*ones(M,1),linea{6},'LineWidth',lw)
plot(r(1,:),P_allKO_sensitivity_simulation(1,:)*100,leg{3},'LineWidth',lw,'MarkerSize',ms)
plot(r(1,:),P_allKO_sensitivity_simulation(2,:)*100,leg{4},'LineWidth',lw,'MarkerSize',ms)

title('Mean blockage probability for different positions of the RSs', 'FontSize',fs)
xlabel('Distance from RS to BS, r (m)','FontSize',fs)
ylabel('Mean blockage probability, P(allKO) (%)','FontSize',fs)
axis([0 max(r(1,:)) 10 max(max(P_allKO_sensitivity_simulation))*100]);
legend('Results for no RSs, no sensitivity','Results for h_R=10 m, no sensitivity','Results for h_R=20 m, no sensitivity','Results for no RSs, sensitivity','Results for h_R=10 m, sensitivity','Results for h_R=20 m, sensitivity')
grid on
grid minor
hold off