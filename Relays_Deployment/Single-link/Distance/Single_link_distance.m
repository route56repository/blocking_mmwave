%WE HAVE TO LOAD MANUALLY THE RESULTS

lin{1}='b';marca{1}='bs';leg{1}='bs-';lin{2}='r';marca{2}='r*';leg{2}='r*-';lin{3}='k';marca{3}='ko';leg{3}='ko-';lin{4}='m';marca{4}='mx';leg{4}='mx-';lin{5}='g';marca{5}='g^';leg{5}='g^-';
linea{1}='b--';linea{2}='r--';linea{3}='k--';linea{4}='m--';linea{5}='g--';
lw=2;  %Fat lines (such as line_width=2) are nice in the pdf. To see how good are bounds, set lw to 1
fs=12; %Big font size (such as font size =12) are nice in the pdf. It can be changed to 10, to see figures clearly when working
ms=8; %%Big marker size (such as marker size =8) are nice in the pdf.
md=3; %Distance between markers

figure(1)
hold on
plot(d,P_KO_simulation_1*100,lin{1},'LineWidth',lw,'HandleVisibility','off'); plot(d(1:md:end),P_KO_simulation_1(1:md:end)*100,marca{1},'LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,leg{1},'LineWidth',2,'MarkerSize',ms);
plot(d,P_KO_calculus_1*100,linea{1},'LineWidth',lw);
plot(d,P_KO_simulation_2_2*100,lin{2},'LineWidth',lw,'HandleVisibility','off'); plot(d(1:md:end),P_KO_simulation_2_2(1:md:end)*100,marca{2},'LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,leg{2},'LineWidth',2,'MarkerSize',ms);
plot(d,P_KO_calculus_2_2*100,linea{2},'LineWidth',lw);
title('Blockage probability for different blockage densities', 'FontSize',fs)
xlabel('Distance from UE to BS, d (m)','FontSize',fs)
ylabel('Blockage probability, P(KO) (%)','FontSize',fs)
axis([0 max(d) 0 max([max(P_KO_simulation_2_2) max(P_KO_calculus_2_2)])*100]);
legend('Simulated results for \lambda=1e-4 m^{-2}','Analytical results for \lambda=1e-4 m^{-2}','Simulated results for \lambda=2.2e-4 m^{-2}','Analytical results for \lambda=2.2e-4 m^{-2}')
grid on
grid minor
hold off