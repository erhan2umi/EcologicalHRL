function x = show_stats_disc_succnorm_sem(foldername, repeats)
    N = 84;
    COMP=1; SAFETY=2; FALL=3;
    if ~exist('foldername','var')
        foldername = 'keep1_REPS2021_4_21_15_0';repeats = 60;    % ERH2025 added this
        fprintf('\n** Assuming the folder is %s with repeat = %d **\n', foldername, repeats);

        fprintf('Need to supply a folder name and a repeat count. See repeat_exp_v1.m\n');
        %foldername = './EXP_REPS2020_12_3_12_44';repeats = 25;
        %foldername = './EXP_REPS2020_12_3_13_23'; repeats = 40;
        %%foldername = './EXP_REPS2020_12_3_16_57'; repeats = 30;
        %foldername = './EXP_REPS2020_12_3_17_24'; repeats=30;
        %%return
    end
    ALLEX = [];
    filenamebase = 'expsim';

    TAnorm =zeros(repeats, N);
    ITAnorm =zeros(repeats, N);
    sfs_ITA_incr_cnt=[];
    sfs_ITA_incr_amt=[];
    sfs_TA_incr_cnt=[];
    sfs_TA_incr_amt=[];

    num_pert_FAIL = [];
    num_pert_SUCC = [];
    k = 0;
    TA_M = zeros(repeats, N);
    ITA_M = TA_M;
    cost_M = TA_M;      % 
    FAIL_M = TA_M+1;    % mark all fails
    FBMUL_M = TA_M+nan;     % feedback gains
    SC = ones(repeats,1);
    realN = zeros(repeats,1);
    okix = [];
    succ_pert_fix=[]; succ_pert_lix=[];
    experL=[];
    first_succ_ix_M = [];
    
    NUMCH = 3;  %this is for the first-3 last-3 succ and fail comparsion
     
    for k_it=1:repeats
        fname = sprintf('%s/%s_rep%d',foldername, filenamebase,k_it);
        DD = load(fname);
        exper = DD.exper;
        experL{k_it} = exper;
        num_trials = exper.lasttrialno;
        %showexper(exper,333);
      

        TA = exper.h_TA(exper.PERTix);
        if (max(TA)>4000)
            fprintf('   ### DISCARDING: Sim %d exploded (huge TA:%4.2f).\n',k_it, max(TA));
            continue;
        elseif max(exper.h_costs(:,COMP))>500
            fprintf('   ### DISCARDING: %d exploded (huge cost:%4.2f) .\n',k_it,exper.h_costs(:,COMP));
            continue;
        elseif isempty(exper.PERTix)
            fprintf('   ### DISCARDING: %d: PERTix empty (so, the required # of successes cannot be reached)!.\n',k_it );
            continue;
        else
            k = k + 1;
            okix(k)= k_it;
        end
        
        % -----------------------------
        %First 3 last 3 for TA (succ and fail)
        Vf = exper.h_TA( exper.fail_pert_ix);  
        Vs = exper.h_TA( exper.succ_pert_ix);
        fail_TA_first3 = mean( Vf(1:NUMCH));
        fail_TA_last3 = mean(Vf(end-NUMCH+1:end));
        succ_TA_first3 = mean( Vs(1:NUMCH));
        succ_TA_last3 = mean(Vs(end-NUMCH+1:end));
        % First 3 last 3 for ITA (succ and fail)
        Vf = exper.h_ITA( exper.fail_pert_ix);  
        Vs = exper.h_ITA( exper.succ_pert_ix);
        fail_ITA_first3 = mean( Vf(1:NUMCH));
        fail_ITA_last3 = mean(Vf(end-NUMCH+1:end));
        succ_ITA_first3 = mean( Vs(1:NUMCH));
        succ_ITA_last3 = mean(Vs(end-NUMCH+1:end));
        
        fail_TA_FL(k,:) = [ fail_TA_first3 , fail_TA_last3 ];
        succ_TA_FL(k,:) = [ succ_TA_first3 , succ_TA_last3 ];
        
        fail_ITA_FL(k,:) = [ fail_ITA_first3 , fail_ITA_last3 ];
        succ_ITA_FL(k,:) = [ succ_ITA_first3 , succ_ITA_last3 ];
        %--------------------------------------
        
        numpert = exper.condTrialCounts(2);
        Rix = 1:numpert;
        Nix = 1: N;
        SC = N/(numpert-1);
        realN(k) = num_trials;
        rSC(k) = SC;
        
       % exper.PERTix(1)

        
        rix = exper.PERTix(1) + round(Nix/SC) - 1;
                % ERHAN temp
        %figure(145);clf; plot(exper.h_ITA(rix),'Linewidth',1); hold on; plot(exper.h_TA(rix), 'r.','Linewidth',2)
       % BELOW hack is for fixing the issue that when subjects fall
       % backward their TA calculation should not include the point after
       % they lose balance (i.e. COM projection leaves the support polygon)
       % Need to run the simulation again to fix this properly. 
        TA_M(k,:)   = exper.h_TA(rix); negfallix = TA_M(k,:)<-400;  TA_M(k,negfallix) = -400; 
        ITA_M(k,:)  = exper.h_ITA(rix);
        cost_M(k,:) = exper.h_costs(rix,COMP);
        FAIL_M(k,:) = exper.h_failed(rix);
        FBMUL_M(k,:) = exper.h_fb_MUL(rix); %xxx
        
        f = FAIL_M(k,:);
        c = cost_M(k,:);
        succ_cost_M(k,:) = cost_M(k,:);
        succ_cost_M(k,f==1)=nan;
        fail_cost_M(k,f==0)=nan;
        
        six = find(FAIL_M(k,:) ~= 1);
        fix = find(FAIL_M(k,:) == 1);
        first_succ_ix_M(k) = six(1);
        
        succ_TA_M(k,:) = TA_M(k,:);
        succ_TA_M(k,f==1)=nan;
        succ_ITA_M(k,:) = ITA_M(k,:);
        succ_ITA_M(k,f==1)=nan;
        

        
        pert_ix = exper.PERTix;   
        succ_pert_ix = pert_ix(exper.h_failed(pert_ix)~=1);
        fail_pert_ix = pert_ix(exper.h_failed(pert_ix)==1);
        
        num_pert_FAIL(k) = length(fail_pert_ix);
        num_pert_SUCC(k)  = length(succ_pert_ix);
        
       % first_pertsuccess_ix(k) = succ_pert_ix(1) - pert_ix(1);
        succ_pert_fix = [succ_pert_fix, succ_pert_ix(1)];
        succ_pert_lix = [succ_pert_lix, succ_pert_ix(end)];

        %succ_ITA(k,:) = resampAT(exper.h_ITA(psix), psix,  N);
        %fail_ITA(k,:) = resampAT(exper.h_ITA(pfix), pfix,  N);

        %succ_TA(k,:) = resampAT(exper.h_TA(psix), psix, N);
        %fail_TA(k,:) = resampAT(exper.h_TA(pfix), pfix, N);
        %figure(44); clf;
        %plot(pert_ix,exper.h_cost(pert_ix),'r.-'); hold on;
        %plot(succ_pert_ix, exper.h_cost(succ_pert_ix),'b.-','Linewidth',1.5); hold on;

        %drawnow;   
        %exper.succ_pert_ix = succ_pert_ix;
        %exper.fail_pert_ix = fail_pert_ix;
        %exper.succ_pert_cost = exper.h_cost(succ_pert_ix);

        ALLEX{k}=exper;

        [sfs_delTA,  sfs_delITA, sfs_delcost, sfs_delix] =  check_SFS(exper);
        sfs_ITA_incr_cnt(k) = sum(sfs_delITA>0) - sum(sfs_delITA<0);
        sfs_TA_incr_cnt(k) = sum(sfs_delTA>0) - sum(sfs_delTA<0)  ;
        sfs_cost_incr_cnt(k) = sum(sfs_delcost>0) - sum(sfs_delcost<0)  ;
        sfs_ITA_incr_amt(k) = sum(sfs_delITA); 
        sfs_TA_incr_amt(k) = sum(sfs_delTA); 
        sfs_cost_incr_amt(k) = sum(sfs_delcost);
        
        fprintf('EXP:[%d/%d] > SFS ITA_incr_cnt:%d   TA_incr_cnt:%d cost_incr_cnt:%d\n',k,k_it,sfs_ITA_incr_cnt(k), sfs_TA_incr_cnt(k), sfs_cost_incr_cnt(k));
        %fprintf('EXP:%d > ITA_incr_amt:%d   TA_incr_amt:%d\n',k,ITA_incr_amt, TA_incr_amt);
        fff=5;

        [ss_delTA,  ss_delITA, ss_delcost, ss_delix] =  check_SS(exper);
        ss_ITA_incr_cnt(k) = sum(ss_delITA>0) - sum(ss_delITA<0);
        ss_TA_incr_cnt(k) = sum(ss_delTA>0) - sum(ss_delTA<0)  ;
        ss_cost_incr_cnt(k) = sum(ss_delcost>0) - sum(ss_delcost<0)  ;
        ss_ITA_incr_amt(k) = sum(ss_delITA);
        ss_TA_incr_amt(k) = sum(ss_delTA);
        ss_cost_incr_amt(k) = sum(ss_delcost); 
        fprintf('EXP:[%d/%d] > SS ITA_incr_cnt:%d   TA_incr_cnt:%d cost_incr_cnt:%d\n\n',k,k_it,ss_ITA_incr_cnt(k), ss_TA_incr_cnt(k),ss_cost_incr_cnt(k));
    
        all_ITA_incr_amt(k) = sfs_ITA_incr_amt(k) + ss_ITA_incr_amt(k);  % is equal to ITA@last_succ - ITA@first_succ
        all_TA_incr_amt(k) = sfs_TA_incr_amt(k) + ss_TA_incr_amt(k);% is equal to TA@last_succ - TA@first_succ
        last_succ_pair = max([sfs_delix; ss_delix]);
        first_succ_pair = min([sfs_delix; ss_delix]);
        last_succ_ix = last_succ_pair(2); first_succ_ix = first_succ_pair(1);
        check_all_TA(k) = exper.h_TA(last_succ_ix) - exper.h_TA(first_succ_ix);
        check_all_ITA(k) = exper.h_ITA(last_succ_ix) - exper.h_ITA(first_succ_ix);
         %plot( exper.h_TA(first_succ_ix:last_succ_ix))
         
         %pert_TA = exper.h_TA(exper.PERTix);
         %pert_succ_TA =  exper.h_TA(exper.h_fallcost<1e-5);
         %allTA = [check_all_TA(k);
         %         sum(pert_TA(2:end)-pert_TA(1:end-1));
         %         exper.h_TA(first_succ_ix) ]
    end  % for k_it
    
    OEC = length(okix);  % this is the number of experiments used in stats.

    avgsuccix = floor(mean(first_succ_ix_M(1:OEC)));
    for k = 1:OEC
        firsuccix = first_succ_ix_M(k);
        
        sft_succ_TA_M = succ_TA_M(k,:)*nan;
        succthatcanfit = succ_TA_M(k,firsuccix:N-max(0,avgsuccix-firsuccix));
        sft_succ_TA_M(avgsuccix:min(avgsuccix + length(succthatcanfit)-1,N)) = succthatcanfit;
        succ_TA_M(k,:) = sft_succ_TA_M;
        
        sft_succ_ITA_M = succ_ITA_M(k,:)*nan;
        succthatcanfit = succ_ITA_M(k,firsuccix:N-max(0,avgsuccix-firsuccix));
        sft_succ_ITA_M(avgsuccix:min(avgsuccix + length(succthatcanfit)-1,N)) = succthatcanfit;
        succ_ITA_M(k,:) = sft_succ_ITA_M;
        
        sft_succ_cost_M = succ_cost_M(k,:)*nan;
        succthatcanfit = succ_cost_M(k,firsuccix:N-max(0,avgsuccix-firsuccix));
        sft_succ_cost_M(avgsuccix:min(avgsuccix + length(succthatcanfit)-1,N)) = succthatcanfit;
        succ_cost_M(k,:) = sft_succ_cost_M;

    end
    

% %     % This reconfirms the indexing can be removed after having stable code 
% %     figure(134);  clf;
% %     subplot(2,1,1);
% %     plot(check_all_ITA,'k'); hold on;
% %     for k=1:OEC
% %         %plot( k, ALLEX{k}.h_ITA(succ_pert_fix(k)),'ro');  
% %         %plot( k, ALLEX{k}.h_ITA(succ_pert_lix(k)),'bo'); 
% %         plot( k, ALLEX{k}.h_ITA(succ_pert_lix(k)) - ALLEX{k}.h_ITA(succ_pert_fix(k)),'go'); 
% %     end
% %     title('ITA last-first success');
% %     
% %     subplot(2,1,2);
% %     plot(check_all_TA,'k'); hold on;
% %     for k=1:OEC
% %         %plot(exper.h_TA(succ_pert_fix),'r-'); hold on;  
% %         %plot(exper.h_TA(succ_pert_lix),'b-'); 
% %         plot( k, ALLEX{k}.h_TA(succ_pert_lix(k)) - ALLEX{k}.h_TA(succ_pert_fix(k)),'go'); 
% %     end
% %     title('TA last-first success');
    
   
    for i=1:N
        f = FAIL_M(1:OEC,i);
        c = cost_M(1:OEC,i);
        ita = ITA_M(1:OEC,i);
        ta = TA_M(1:OEC,i);
        fail_cnt(i) = sum(f==1);
        fail_cost = c(f==1);
        fail_cost_avg(i) = mean(fail_cost);
        fail_cost_std(i) = std(fail_cost);
        fail_cost_sem(i) = fail_cost_std(i)/sqrt(length(fail_cost));
        succ_cnt(i) = sum(f==0);
        
        %succ_cost = c(f==0);   % this does not see the avg norm. for succ
        succ_cost = succ_cost_M(1:OEC,i); succ_cost=succ_cost(~isnan(succ_cost));
        succ_cost_avg(i) = mean(succ_cost);
        succ_cost_std(i) = std(succ_cost);
        succ_cost_sem(i) = succ_cost_std(i)/sqrt(length(succ_cost));
        
        %succ_ITA = ita(f==0); % this does not see the avg norm. for succ
        succ_ITA = succ_ITA_M(1:OEC,i); succ_ITA=succ_ITA(~isnan(succ_ITA));
        succ_ITA_avg(i) = mean(succ_ITA);
        succ_ITA_std(i) = std(succ_ITA);
        succ_ITA_sem(i) = succ_ITA_std(i)/sqrt(length(succ_ITA));
         
        fail_ITA = ita(f==1);
        fail_ITA_avg(i) = mean(fail_ITA);
        fail_ITA_std(i) = std(fail_ITA);
        fail_ITA_sem(i) = fail_ITA_std(i)/sqrt(length(fail_ITA));

        
        %succ_TA = ta(f==0); % this does not see the avg norm. for succ
        succ_TA = succ_TA_M(1:OEC,i); succ_TA=succ_TA(~isnan(succ_TA));
        succ_TA_avg(i) = mean(succ_TA);
        succ_TA_std(i) = std(succ_TA);
        succ_TA_sem(i) = succ_TA_std(i)/sqrt(length(succ_TA));
        
        fail_TA = ta(f==1);
        succ_TA = succ_TA_M(1:OEC,i); succ_TA=succ_TA(~isnan(succ_TA));
        fail_TA_avg(i) = mean(fail_TA);
        fail_TA_std(i) = std(fail_TA);
        fail_TA_sem(i) = fail_TA_std(i)/sqrt(length(fail_TA));
    end

    
    all_ = 1;
    ss_  = 2;
    sfs_ = 3;

    %del_ita = ita(2:end)-ita(1:end-1);
    %del_ta  =  ta(2:end)- ta(1:end-1);
    
    lens = [length(ita), length(ss_ITA_incr_amt), length(sfs_ITA_incr_amt)];
        
    % average over the repeats to get mean increment in (ITA and TA) of areas on ss and sfs
    ita_means   = [mean(all_ITA_incr_amt), mean(ss_ITA_incr_amt), mean(sfs_ITA_incr_amt)];
    ita_devs    = [ std(all_ITA_incr_amt),  std(ss_ITA_incr_amt),  std(sfs_ITA_incr_amt)]./lens.^0.5;
    ta_means    = [ mean(all_TA_incr_amt), mean(ss_TA_incr_amt), mean(sfs_TA_incr_amt)];
    ta_devs     = [ mean(all_TA_incr_amt),  std(ss_TA_incr_amt), std(sfs_TA_incr_amt)]./lens.^0.5;
    
    
    fprintf('\nMean ITA SFS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean(sfs_ITA_incr_cnt),std(sfs_ITA_incr_cnt),ita_means(sfs_),ita_devs(sfs_));
    fprintf('Mean  TA SFS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean( sfs_TA_incr_cnt),std(sfs_TA_incr_cnt),mean( sfs_TA_incr_amt),std(sfs_TA_incr_amt));
    fprintf('Mean  Cost SFS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean( sfs_cost_incr_cnt),std(sfs_cost_incr_cnt),mean( sfs_cost_incr_amt),std(sfs_cost_incr_amt));
    
    fprintf('\nMean ITA SS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean(ss_ITA_incr_cnt),std(ss_ITA_incr_cnt),mean(ss_ITA_incr_amt),std(ss_ITA_incr_amt));
    fprintf('Mean  TA SS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean( ss_TA_incr_cnt),std(ss_TA_incr_cnt),mean( ss_TA_incr_amt),std(ss_TA_incr_amt));
    fprintf('Mean  Cost SS increment [cnt:%3.3f+-%3.3f]   [amount:%3.3f+-%3.3f]\n',mean( ss_cost_incr_cnt),std(ss_cost_incr_cnt),mean( ss_cost_incr_amt),std(ss_cost_incr_amt));

    fprintf('Succes in PERT reached on %2.2f +- %2.2f trials (Total trials: OK/Attempted=%d/%d)\n',mean(succ_pert_fix),std(succ_pert_fix),OEC,repeats);

    fprintf('Number of fails in PERT trials %2.2f +- %2.2f trials (Total trials: OK/Attempted=%d/%d)\n',mean(num_pert_FAIL), std(num_pert_FAIL),OEC,repeats);
    allcol  = [0 0 0];
    failcol = [1 0 0];
    succcol = [0 0.9 0.2];
    
    ntr = 1:N; 
    eps = 0.2;

    figure(46); clf;
    subplot(1,3,1);
    %plot_mean_std_shade(ntr, mean(TA_M), sem(TA_M), allcol,whiten(allcol,0.8),3); hold on;
    plot_mean_std_nan(ntr, mean(TA_M), sem(TA_M), allcol,whiten(allcol,0.5),3);hold on;
    plot_mean_std_nan(ntr+2*eps, fail_TA_avg, fail_TA_sem, failcol,whiten(failcol,0.6), 1.5);
    plot_mean_std_nan(ntr+eps, succ_TA_avg, succ_TA_sem, succcol,whiten(succcol,0.1), 2);
        
    title('Pert cond: Average TA'); %legend('','All','', 'Success', '','Fail');
    xlabel('normalized trial no'); ylabel('TA(cm^2)');axis tight;
    yl = ylim; ylim([max(yl(1),-500), min(yl(2),900)]);
    subplot(1,3,2);
    %plot_mean_std_shade(ntr, mean(ITA_M), sem(ITA_M), allcol,whiten(allcol,0.8),3); hold on;
    plot_mean_std_nan(ntr, mean(ITA_M), sem(ITA_M), allcol,whiten(allcol,0.5),3); hold on;
    plot_mean_std_nan(ntr+2*eps, fail_ITA_avg, fail_ITA_sem, failcol,whiten(failcol,0.6), 1.5);
    plot_mean_std_nan(ntr+eps, succ_ITA_avg, succ_ITA_sem, succcol,whiten(succcol,0.1), 2);
    title('Pert cond: Average ITA'); %legend('','All','', 'Success', '','Fail');
    xlabel('normalized trial no'); ylabel('ITA(cm^2)');axis tight;

    subplot(1,3,3);
    %plot_mean_std_shade(ntr, mean(cost_M), sem(cost_M), allcol,whiten(allcol,0.8),3); hold on;
    plot_mean_std_nan(ntr, mean(cost_M), sem(cost_M), allcol,whiten(allcol,0.5),3); hold on;
    plot_mean_std_nan(ntr+2*eps, fail_cost_avg, fail_cost_sem, failcol,whiten(failcol,0.6), 1.5);
    plot_mean_std_nan(ntr+eps, succ_cost_avg, succ_cost_sem, succcol,whiten(succcol,0.1), 2);
    title(sprintf('Pert cond: Avg. CompCost  %1.2f*effort + %1.2f*safety', exper.body.effortW,1-exper.body.effortW));
    xlabel('normalized trial no'); axis tight

    
    figure (45); clf;
    subplot(1,2,1);
    XX = 1:3;
    h=bar(XX, diag(ita_means),'stacked'); hold on; er = errorbar(XX,ita_means,ita_devs/2,ita_devs/2); 
    set(h(1),'facecolor',[0.6 0.6 0.6]);
    set(h(2),'facecolor',[1 0.6 0.0]);
    set(h(3),'facecolor',[0 0.1 .7]);
    set(gca,'xtick',XX,'xticklabel',{'all','S-S','S-F-S'})
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    ylabel('delta ITA');
    title('mean ITA change /w SEM');
    subplot(1,2,2);
    h=bar(XX, diag(ta_means),'stacked'); hold on; er = errorbar(XX,ta_means,ta_devs/2,ta_devs/2); 
    set(h(1),'facecolor',[0.6 0.6 0.6]);
    set(h(2),'facecolor',[1 0.6 0.0]);
    set(h(3),'facecolor',[0 0.1 .7]);
    set(gca,'xtick',XX,'xticklabel',{'all','S-S','S-F-S'})
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    ylabel('delta TA');
    title('mean TA change /w SEM');
    %return
    
    FAIL_M = FAIL_M(1:OEC,:);   % remove the last rows due to discards
    FBMUL_M = FBMUL_M(1:OEC,:);
    SUCC_M = 1-FAIL_M;
    figure(35); clf; 
    subplot(2,1,1);
    plot(sum(FAIL_M),'r.');  hold on; 
    %plot(sum(SUCC_M),'g.');
    xlabel('normazlized Tiral no'); ylabel('number of failures'); %legend('Fails','Succ.');
    title('Failures over the multiple experiments'); 
     subplot(2,1,2);
    plot_mean_std_shade(ntr, mean(FBMUL_M), std(FBMUL_M), allcol,whiten(allcol,0.8),3);
    plot_mean_std_shade(ntr, mean(FBMUL_M.*(1-FAIL_M)), std(FBMUL_M), [0 .8 0.2],whiten([0 .8 0.2],0.8),3);
    plot_mean_std_shade(ntr, mean(FBMUL_M.*FAIL_M), std(FBMUL_M), [1 0.2 0.2],whiten([1 0.2 0.2],0.8),3);
    
    title('FB gain multiplier');
    for selk =18: -18 % repeats, % cmink
        mysix= ALLEX{selk}.succ_pert_ix;
        myfix= ALLEX{selk}.fail_pert_ix;
        myix = ALLEX{selk}.PERTix;
        fpi = myix(1)+1;
        
        [delTA, delITA, delcost, delix] = check_SFS(ALLEX{selk});
        subplot(2,3,4); 
        plot(myix-fpi, ALLEX{selk}.h_TA(myix), 'Color',[0 0 0],'Linewidth',1.2); hold on;
        plot(mysix-fpi, ALLEX{selk}.h_TA(mysix), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
        plot(myfix-fpi, ALLEX{selk}.h_TA(myfix),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;
        legend('All','Fail', 'Success','Location','southeast');
        title(sprintf('Trial %d: TA (sum:%3.3f)',selk,sum(delTA) ));
        axis tight;
        subplot(2,3,5); 
        plot(myix-fpi, ALLEX{selk}.h_ITA(myix), 'Color',[0 0 0],'Linewidth',1.2); hold on;
        plot(mysix-fpi, ALLEX{selk}.h_ITA(mysix), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
        plot(myfix-fpi, ALLEX{selk}.h_ITA(myfix),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;
        plot(delix-fpi, ALLEX{selk}.h_ITA(delix),'o','MarkerSize',12,'Color',[0 0 1]);
        plot(delix(delITA>=0)-fpi, ALLEX{selk}.h_ITA(delix(delITA>=0)),'+','MarkerSize',12,'Color',[0 0 1]);
        legend('All','Fail', 'Success','Location','southeast');
        title(sprintf('Trial %d: ITA (sum:%3.3f)',selk,sum(delITA) ));
        axis tight;
        %title('TA - success');
        subplot(2,3,6);
        plot(myix-fpi, ALLEX{selk}.h_costs(myix,COMP), 'Color',[0 0 0],'Linewidth',1.2); hold on;
        plot(mysix-fpi, ALLEX{selk}.h_costs(mysix,COMP), '.','Color',[0.2 0.8 0.2],'MarkerSize',15); hold on;
        plot(myfix-fpi, ALLEX{selk}.h_costs(myfix,COMP),'.', 'Color',[1 0 0],'MarkerSize',15); hold on;
        legend('All','Fail', 'Success','Location','northeast');
        title(sprintf('Trial %d:  Cost',selk));
        axis tight;
        %input('Press enter for next trial..');  subplot(2,3,4); cla; subplot(2,3,5); cla; subplot(2,3,6); cla; 
    end
    
    
    [H,p] = ttest(fail_ITA_FL(:,1),fail_ITA_FL(:,2));
    fprintf('ITA FAILS: first vs. last trials: [%2.3f +- %2.3f]   [%2.3f +- %2.3f] (H:%d, p:%1.5f)\n', mean(fail_ITA_FL(:,1)),std(fail_ITA_FL(:,1)),mean(fail_ITA_FL(:,2)),std(fail_ITA_FL(:,2)), H, p);
    [H,p] = ttest(succ_ITA_FL(:,1),succ_ITA_FL(:,2));
    fprintf('ITA SUCCS: first vs. last trials: [%2.3f +- %2.3f]   [%2.3f +- %2.3f] (H:%d, p:%1.5f)\n', mean(succ_ITA_FL(:,1)),std(succ_ITA_FL(:,1)),mean(succ_ITA_FL(:,2)),std(succ_ITA_FL(:,2)), H, p);
    
    [H,p] = ttest(fail_TA_FL(:,1),fail_TA_FL(:,2));
    fprintf(' TA FAILS: first vs. last trials: [%2.3f +- %2.3f]   [%2.3f +- %2.3f] (H:%d, p:%1.5f)\n', mean(fail_TA_FL(:,1)),std(fail_TA_FL(:,1)),mean(fail_TA_FL(:,2)),std(fail_TA_FL(:,2)), H, p);
    [H,p] = ttest(succ_TA_FL(:,1),succ_TA_FL(:,2));    
    fprintf(' TA SUCCS: first vs. last trials: [%2.3f +- %2.3f]   [%2.3f +- %2.3f] (H:%d, p:%1.5f)\n', mean(succ_TA_FL(:,1)),std(succ_TA_FL(:,1)),mean(succ_TA_FL(:,2)),std(succ_TA_FL(:,2)), H, p);
    
    
    
end %of function show_stats_v2()


function softcol = whiten(col, w)
    if ~exist('w','var'), w=0.5; end
    softcol = (1-w)*col + w*(col*0+1);
end
 

function B = resamp(A, newsize)
    tr = (1:length(A))/length(A);
    sp = spline(tr, A);
    newtr = (1:newsize)/newsize;
    B = ppval(sp, newtr);
end % of resamp

function B = resampAT(A,T,N)
size(A)
size(T)
N
    sp = spline([0, T]/N, [A(1),A]);
    B = ppval(sp, [0:N-1]/N);
end % of resamp

function B = remap(A, newsize)
    ix = 1:length(A);
    SC = length(A)/newsize;
    B = zeros(newsize,1);
    for k=1:newsize
        acc = max(round(k*SC), 1);
        B(k) = A(acc);
    end
end

function plot_mean_std_nan(x, y, dev, col, shadecol, linewidth)
    gi = ~isnan(y);
    x = x(gi);
    y = y(gi);
    dev = dev(gi);
    plot_mean_std(x, y, dev, col, shadecol, linewidth);
end

function plot_mean_std_shade(x, y, dev, col, shadecol, linewidth)
    curve1 = y + dev;
    curve2 = y - dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, shadecol, 'FaceAlpha',0.5,'LineStyle','none');
    hold on;
    plot(x, y, 'r', 'LineWidth', linewidth,'Color',col);
end

function plot_mean_std(x, y, dev, col, shadecol, linewidth)
    curve1 = y + dev;
    curve2 = y - dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    %fill(x2, inBetween, shadecol, 'FaceAlpha',0.5,'LineStyle','none');
    hold on;
    plot([x;x],[curve1;curve2],'Color',shadecol, 'Linewidth',0.75)
    plot(x, y, '-', 'LineWidth', linewidth,'Color',col);
    %plot(x, curve1, '.', 'LineWidth', 2,'Color',col);
    %plot(x, curve2, '.', 'LineWidth', 2,'Color',col);
end

function plot_mean_std_errorbar_fix(x, y, dev, col, shadecol, linewidth)
    curve1 = y + dev;
    curve2 = y - dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    %fill(x2, inBetween, shadecol, 'FaceAlpha',0.5,'LineStyle','none');
    hold on;
    errorbar([x;x]',[curve1;curve2]','Color',shadecol, 'Linewidth',0.75)
    plot(x, y, '-', 'LineWidth', linewidth,'Color',col);
    %plot(x, curve1, '.', 'LineWidth', 2,'Color',col);
    %plot(x, curve2, '.', 'LineWidth', 2,'Color',col);
end

function sem_val = sem(x)
    sem_val = std(x)/sqrt(length(x));
end

function showexper(exper, figno)
    COMP=1; SAFETY=2; FALL=3;
    it = exper.lasttrialno;
    figure(figno); clf;
    subplot(3,1,1); cla;
    %plot(exper.h_cost, 'Color', [0.6 0.6 0.6],'Linewidth',2);  hold on; 

    %bar((exper.h_fallcost(1:it)~=0)*max([0.1;exper.h_baseline_cost(1:it)+exper.h_adv(1:it)]),'FaceColor',[1 0.7 0.7],'EdgeColor','none'); hold on

    bar((exper.h_failed(1:it))*max([0.1;exper.h_baseline_costs(1:it,1)+exper.h_costadvs(1:it,1)]),'FaceColor',[1 0.7 0.7],'EdgeColor','none'); hold on
    plot(exper.h_baseline_costs(1:it,1), 'Color', [0.2 0.2 0.2],'Linewidth',3);  hold on;
    plot(exper.h_baseline_costs(1:it,1)+exper.h_costadvs(1:it,1),'.-','Color',[0 0 1],'Linewidth',1);
    %plot(exper.h_baseline_cost+exper.h_adv,'.','Color',[0 0 0]);
    plot(exper.h_fb_MUL(1:it)/120,'-','Linewidth',2, 'Color',[0.1 0.7 0.5]);
    ylim([0*min(exper.h_baseline_costs(:,1)+exper.h_costadvs(:,1)),max(exper.h_baseline_costs(:,1)+exper.h_costadvs(:,1))]);
    grid on;
    %xlabel('trial no');
    legend('falls','baseline_ccost','advantage','fpMUL/120')


    subplot(3,1,2); cla;
    plot(exper.unPERTix,exper.h_TA(exper.unPERTix), 'Color', [1 0.8 0.0],'Linewidth',2); hold on;
    plot(exper.PERTix, exper.h_TA(exper.PERTix), 'Color', [1 0.5 0.0],'Linewidth',2); hold on;
    plot(exper.catchtrial_no, exper.h_TA(exper.catchtrial_no),'r.');
    plot(exper.PERTxix, exper.h_TA(exper.PERTxix), 'Color', [1 0.5 0.0],'Linewidth',2); hold on;
    plot(exper.DEADix, exper.h_TA(exper.DEADix), 'Color', [1 0.8 0.0],'Linewidth',2); hold on;
    plot(exper.h_ITA(1:it)*100, 'Color', [0 0.6 0.0],'Linewidth',2); hold on;
    ylabel('cm^2');
    %plot(exper.h_cost(1:it)*10, 'Color', [0 0.2 0.9],'Linewidth',2); hold on;
    %ylabel('fatness'); 
    %xlabel('trial no'); 
    xlim([0,exper.lasttrialno])
    legend('TA_{unp}','TA_{per}','CT','TA_{perX}','TA_{dead}','100*ITA');

    grid on;

    subplot(3,1,3); cla;
    tmp1 = max([exper.h_costs(1:it,COMP);exper.h_costs(1:it,COMP)+exper.h_costadvs(1:it,COMP)]);
    bar(exper.h_failed(1:it)*max([0.02; 1.1*exper.h_costs(1:it,FALL);tmp1]),'FaceColor',[1 0.7 0.7],'EdgeColor','none'); hold on;
    bar(exper.catchtrial_no, max([0.02; 1.1*exper.h_costs(1:it, FALL);tmp1]),0.5,'FaceColor',[1 0 0]);
    bar((1-exper.h_failed(1:it))*max([0.02;1.1*exper.h_costs(1:it, FALL);tmp1]),'FaceColor',[0.8 1 1],'EdgeColor','none'); hold on;
    bar((1-exper.h_failed(exper.unPERTix))*max([0.02;1.1*exper.h_costs(1:it, FALL);tmp1]),'FaceColor',[0.6 0.9 0.6],'EdgeColor','none'); hold on;
    bar(exper.DEADix,(1-exper.h_failed(exper.DEADix))*max([0.02;1.1*exper.h_costs(1:it,FALL);tmp1]),'FaceColor',[0.6 0.9 0.6],'EdgeColor','none'); hold on;

    bar(exper.h_costs(1:it,FALL),0.5 ,'FaceColor',[0.6 0.6 0.6]); hold on;
    bar(exper.DEADix, exper.h_costs(exper.DEADix,FALL),0.5 ,'FaceColor',[0.9 0.9 0.9]); hold on;
    bar(exper.unPERTix, exper.h_costs(exper.unPERTix,FALL),0.5 ,'FaceColor',[0.9 0.9 0.9]); hold on;
    plot(exper.h_baseline_costs(1:it,FALL),'y-','Linewidth',2)
    plot(exper.h_costs(1:it,COMP), 'Color', [0 0.2 0.9],'Linewidth',2); 

    xlim([0,exper.lasttrialno])
    legend('fall','catch','success','unper_0','unper_1','fallcost','fallcost_{unp}','fallcost_{unp}','baselnfallcost','ccost'); xlabel('trial no');
end % of showexpr

