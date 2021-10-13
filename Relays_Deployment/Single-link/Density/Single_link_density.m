%WE HAVE TO LOAD MANUALLY THE RESULTS

lin{1}='b';marca{1}='bs';leg{1}='bs-';lin{2}='r';marca{2}='r*';leg{2}='r*-';lin{3}='k';marca{3}='ko';leg{3}='ko-';lin{4}='m';marca{4}='mx';leg{4}='mx-';lin{5}='g';marca{5}='g^';leg{5}='g^-';
linea{1}='b--';linea{2}='r--';linea{3}='k--';linea{4}='m--';linea{5}='g--';
lw=2;  %Fat lines (such as line_width=2) are nice in the pdf. To see how good are bounds, set lw to 1
fs=12; %Big font size (such as font size =12) are nice in the pdf. It can be changed to 10, to see figures clearly when working
ms=8; %%Big marker size (such as marker size =8) are nice in the pdf.
md=3; %Distance between markers

figure(2)
hold on
plot(lambda,P_KO_simulation_30m*100,leg{1},'LineWidth',lw,'MarkerSize',ms)
plot(lambda,P_KO_calculus_30m*100,linea{1},'LineWidth',lw);
plot(lambda,P_KO_simulation_40m*100,leg{2},'LineWidth',lw,'MarkerSize',ms)
plot(lambda,P_KO_calculus_40m*100,linea{2},'LineWidth',lw);
title('Mean blockage probability for different H_{max}', 'FontSize',fs)
xlabel('Blockage density, \lambda (m^{-2})','FontSize',fs)
ylabel('Mean blockage probability, P(KO) (%)','FontSize',fs,'interpreter','latex')
axis([0 max(lambda) 0 max([max(P_KO_simulation_40m) max(P_KO_calculus_40m)])*100]);
legend('Simulated results for H_{max}=30 m','Analytical results for H_{max}=30 m','Simulated results for H_{max}=40 m','Analytical results for H_{max}=40 m')
grid on
grid minor
hold off