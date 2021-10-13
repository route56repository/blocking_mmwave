%WE HAVE TO LOAD MANUALLY THE RESULTS

lin{1}='b';marca{1}='bs';leg{1}='bs-';lin{2}='r';marca{2}='r*';leg{2}='r*-';lin{3}='k';marca{3}='ko';leg{3}='ko-';lin{4}='m';marca{4}='mx';leg{4}='mx-';lin{5}='g';marca{5}='g^';leg{5}='g^-';
linea{1}='b--';linea{2}='r--';linea{3}='k--';linea{4}='m--';linea{5}='g--';
linia{1}='bs:';linia{2}='r*:';linia{3}='ko:';linia{4}='m:';linia{5}='g^:';
lw=2;  %Fat lines (such as line_width=2) are nice in the pdf. To see how good are bounds, set lw to 1
fs=12; %Big font size (such as font size =12) are nice in the pdf. It can be changed to 10, to see figures clearly when working
ms=8; %%Big marker size (such as marker size =8) are nice in the pdf.
md=3; %Distance between markers

figure(3)
hold on
plot(d,P_allKO_simulation_0rad*100,leg{1},'LineWidth',lw,'MarkerSize',ms)
plot(d,P_allKO_calculus_0rad*100,linea{1},'LineWidth',lw);
plot(d,P_allKO_calculus_ind_0rad*100,linia{1},'LineWidth',lw,'MarkerSize',ms);
plot(d,P_allKO_simulation_pi_12rad*100,leg{2},'LineWidth',lw,'MarkerSize',ms)
plot(d,P_allKO_calculus_pi_12rad*100,linea{2},'LineWidth',lw);
plot(d,P_allKO_calculus_ind_pi_12rad*100,linia{2},'LineWidth',lw,'MarkerSize',ms);
plot(d,P_allKO_simulation_pi_6rad*100,leg{3},'LineWidth',lw,'MarkerSize',ms)
plot(d,P_allKO_calculus_pi_6rad*100,linea{3},'LineWidth',lw);
plot(d,P_allKO_calculus_ind_pi_6rad*100,linia{3},'LineWidth',lw,'MarkerSize',ms);
title('Blockage probability for different observation angles', 'FontSize',fs)
xlabel('Distance from UE to BS, d (m)','FontSize',fs)
ylabel('Blockage probability, P(allKO) (%)','FontSize',fs)
axis([0 max(d) 0 max([max(P_allKO_simulation_0rad) max(P_allKO_calculus_0rad)])*100]);
legend('Simulated results for \phi=0º','Analytical results for \phi=0º','Analytical results for \phi=0º, assuming independence','Simulated results for \phi=15º','Analytical results for \phi=15º','Analytical results for \phi=15º, assuming independence','Simulated results for \phi=30º','Analytical results for \phi=30º','Analytical results for \phi=30º, assuming independence')
grid on
grid minor
hold off