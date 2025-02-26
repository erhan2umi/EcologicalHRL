function analyzeSFS(exper, figno, EXE)
    pert_ix = exper.PERTix;   %
    succ_pert_ix = pert_ix(exper.h_failed(pert_ix)~=1);
    fail_pert_ix = pert_ix(exper.h_failed(pert_ix)==1);
    exper.succ_pert_ix = succ_pert_ix;
    exper.fail_pert_ix = fail_pert_ix;
    exper.succ_pert_cost = exper.h_cost(succ_pert_ix);
    mysix= exper.succ_pert_ix;
    myfix= exper.fail_pert_ix;
    myix = exper.PERTix;
    allix = 1:exper.lasttrialno;
    allfix = allix(exper.h_failed==1);
    allsix = allix(exper.h_failed==1);

    [delTA, delITA, delcost, delix] = check_SFS(exper);
    if ~isempty(delix), delix = delix(:,2)-1; end  % delix is the list of S,S pairs around an F. So -1 gets you to the last F within S's.   
    ITA_incr_cnt = sum(delITA>0) - sum(delITA<0);
    TA_incr_cnt = sum(delTA>0) - sum(delTA<0)  ;
    cost_incr_cnt = sum(delcost>0) - sum(delcost<0)  ;
    ITA_incr_amt = sum(delITA(delITA>0)) + sum(delITA(delITA<0));
    TA_incr_amt = sum(delTA(delTA>0)) + sum(delTA(delTA<0));
    cost_incr_amt = sum(delcost(delcost>0)) + sum(delcost(delcost<0));
    fprintf('\nEXP SFS\n');
    fprintf('**    ITA_incr_cnt:%d  ITA_incr_amt:%3.3f   [cost_incr_cnt:%d]\n',ITA_incr_cnt,ITA_incr_amt, cost_incr_cnt);
    fprintf('**     TA_incr_cnt:%d   TA_incr_amt:%3.3f   [cost_incr_amt:%f]\n',TA_incr_cnt, TA_incr_amt,cost_incr_amt);
    fff=5;
    [ss_delTA, ss_delITA, ss_delcost, ss_delix] = check_SS(exper);
    ss_ITA_incr_cnt = sum(ss_delITA>0) - sum(ss_delITA<0);
    ss_TA_incr_cnt = sum(ss_delTA>0) - sum(ss_delTA<0)  ;
    ss_cost_incr_cnt = sum(ss_delcost>0) - sum(ss_delcost<0)  ;
    ss_ITA_incr_amt = sum(ss_delITA(ss_delITA>0)) + sum(ss_delITA(ss_delITA<0));
    ss_TA_incr_amt = sum(ss_delTA(ss_delTA>0)) + sum(ss_delTA(ss_delTA<0))  ;
    ss_cost_incr_amt = sum(ss_delcost(ss_delcost>0)) + sum(ss_delcost(ss_delcost<0));
    fprintf('\nEXP SS:\n');
    fprintf('**  ITA_incr_cnt:%d  ITA_incr_amt:%3.3f [cost_incr_cnt:%d]\n',ss_ITA_incr_cnt, ss_ITA_incr_amt,ss_cost_incr_cnt);
    fprintf('**   TA_incr_cnt:%d   TA_incr_amt:%3.3f [cost_incr_amt:%f]\n',ss_TA_incr_cnt, ss_TA_incr_amt,ss_cost_incr_amt);


    figure(figno); clf;
    p0=[0 1];
    sctr = (exper.PERTix'-1)/exper.PERTix(end);
    y  = exper.h_ITA(exper.PERTix);
    ITAmin = min(y);
    ITAmax = max(y-ITAmin);
    itay = (y-ITAmin)/ITAmax + 0;
    yy = exper.h_TA(exper.PERTix);
    TAmin = min(yy);
    TAmax = max(yy-TAmin);
    tay = (yy-TAmin)/TAmax + 0;

    [ITAexp_par,ITAexp_err] = fit(sctr,itay','exp1', 'StartPoint', p0 );
    [ TAexp_par, TAexp_err] = fit(sctr,tay,'exp1', 'StartPoint', p0 );

    fprintf('ITA, TA exponent:[%3.3f, %3.3f]\n',ITAexp_par.b,   TAexp_par.b);

    subplot(3,1,1);
    plot(myix, exper.h_TA(myix), 'Color',[0 0 0],'Linewidth',1.2); hold on;
    plot(mysix, exper.h_TA(mysix), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
    plot(myfix, exper.h_TA(myfix),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;

    unpertfellix = find(exper.h_failed(exper.unPERTix)==1);
    unpertsuccix = find(exper.h_failed(exper.unPERTix)==0);
    plot(exper.unPERTix, exper.h_TA(exper.unPERTix),'k-');
    plot(unpertfellix, exper.h_TA(unpertfellix),'.','Color',[1 0.5 0.3],'MarkerSize',15);
    plot(unpertsuccix, exper.h_TA(unpertsuccix),'.','Color',[0.3 1 0.5],'MarkerSize',15);
    plot(myix, (TAexp_par.a*exp(TAexp_par.b*sctr)-0)*TAmax+TAmin,'c-');
    
    frfix = allfix(exper.h_costs(allfix,4) < exper.h_costs(allfix,5));  % front falls 
    plot(frfix, exper.h_TA(frfix)-0.2,'^', 'Color',[0.3 0.5 0.8],'MarkerSize',10); hold on;
    legend('All','Success','Fail','Location','northeast');
    title(sprintf('TA'));
    grid minor;
    subplot(3,1,2);
    plot(myix, exper.h_ITA(myix), 'Color',[0 0 0],'Linewidth',1.2); hold on;
    plot(mysix, exper.h_ITA(mysix), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
    plot(myfix, exper.h_ITA(myfix),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;
    plot(myix, (ITAexp_par.a*exp(ITAexp_par.b*sctr)-0)*ITAmax+ITAmin,'c-');

    plot(delix, exper.h_ITA(delix),'o','MarkerSize',12,'Color',[0 0 1]);
    plot(delix(delITA>=0), exper.h_ITA(delix(delITA>=0)),'+','MarkerSize',12,'Color',[0 0 1]);
    plot(exper.unPERTix, exper.h_ITA(exper.unPERTix),'k-');
    plot(unpertfellix, exper.h_ITA(unpertfellix),'.','Color',[1 0.5 0.3],'MarkerSize',15);
    plot(unpertsuccix, exper.h_ITA(unpertsuccix),'.','Color',[0.3 1 0.5],'MarkerSize',15);

    %frfix = allfix(exper.h_costs(allfix,4) < exper.h_costs(allfix,5));  % front falls 
    plot(frfix, exper.h_ITA(frfix)-0.2,'^', 'Color',[0.3 0.5 0.8],'MarkerSize',10); hold on;

    grid minor;
    legend('All', 'Success','Fail','Location','northeast');
    title('ITA');
    subplot(3,1,3);
    plot(myix, exper.h_costs(myix,EXE.COMP), 'Color',[0 0 0],'Linewidth',1.2); hold on;
    plot(mysix, exper.h_costs(mysix,EXE.COMP), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
    plot(myfix, exper.h_costs(myfix,EXE.COMP),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;

    plot(exper.unPERTix, exper.h_costs(exper.unPERTix, EXE.COMP),'k-');
    plot(unpertfellix, exper.h_costs(unpertfellix, EXE.COMP),'.','Color',[1 0.5 0.3],'MarkerSize',15);
    plot(unpertsuccix, exper.h_costs(unpertsuccix, EXE.COMP),'.','Color',[0.3 1 0.5],'MarkerSize',15);
    %frfix = allfix(exper.h_costs(allfix,4) < exper.h_costs(allfix,5));  % front falls 
    plot(frfix, exper.h_costs(frfix, EXE.COMP)-0.2,'^', 'Color',[0.3 0.5 0.8],'MarkerSize',10); hold on;
    legend('All','Success','Fail','Location','northeast');
    %grid on; 
    grid minor;
    title(sprintf('CompCost  (%1.2f*effort + %1.2f*safety)', exper.body.effortW,1-exper.body.effortW));
    drawnow;