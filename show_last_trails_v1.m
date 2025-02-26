
function show_last_trails_v1(exper)
ff=33;
figure(ff); clf;
subplot(2,3,1);
draw_exe_data_v2(exper.CATCHexe, exper.body, 0);
ylabel(sprintf('CATCH [trno:%d]',exper.catchtrial_no));
subplot(2,3,2);
draw_exe_data_v2(exper.lastunPERTexe, exper.body, 0);
ylabel(sprintf('last unPERT [trno:%d]',exper.unPERTix(end)));
subplot(2,3,3);
draw_exe_data_v2(exper.lastPERTexe, exper.body, 0);
ylabel(sprintf('last PERT [trno:%d]',exper.PERTix(end)));

subplot(2,3,4);
draw_exe_data_v2(exper.firstDEADexe, exper.body, 0);
ylabel(sprintf('first DEADAPT [trno:%d]',exper.DEADix(1)));
subplot(2,3,5);
draw_exe_data_v2(exper.lastDEADexe, exper.body, 0);
ylabel(sprintf('last DEADAPT [trno:%d]',exper.DEADix(end)));

subplot(2,3,6);
draw_exe_data_v2(exper.lastPERTxexe, exper.body, 0);
ylabel(sprintf('last PERTx [trno:%d]',exper.PERTxix(end)));

figure(ff+1); clf;
exe_data = exper.lastunPERTexe;
%plot(exe_data.des_cartpos(:,1), exe_data.des_cartpos(:,2),'.','Linewidth', 2, 'Color',[0 0 0]); hold on;
plot(exe_data.act_cartpos(:,1), exe_data.act_cartpos(:,2),'--','Linewidth', 2, 'Color',[0.5 0.5 0.5]); hold on;
exe_data = exper.lastPERTexe;
%plot(exe_data.des_cartpos(:,1), exe_data.des_cartpos(:,2),'k-','Linewidth', 1.5); hold on;
plot(exe_data.act_cartpos(:,1), exe_data.act_cartpos(:,2),'-','Linewidth', 2, 'Color',[0.5 0.5 0.5]);

exe_data = exper.CATCHexe;
%plot(exe_data.des_cartpos(:,1), exe_data.des_cartpos(:,2),'y-','Linewidth', 1.5); hold on;
plot(exe_data.act_cartpos(:,1), exe_data.act_cartpos(:,2),'-','Linewidth', 2, 'Color',[0.0 0.3 1]);

exe_data = exper.firstDEADexe;
%plot(exe_data.des_cartpos(:,1), exe_data.des_cartpos(:,2),'-','Linewidth', 2, 'Color',[0.9 0.5 0.0]);
plot(exe_data.act_cartpos(:,1), exe_data.act_cartpos(:,2),'-','Linewidth', 2, 'Color',[0.6 0.8 1]);


axis equal;
axis([-0.5,+0.5, 0.45, 1.3]);
legend('unpert','lastpert','catch','1^{st}dead');
%legend('unpert_d','unpert','lastpert','catch','1^{st}dead_d','1^{st}dead');

