function [C, pval] = crosscorr(D1, D2, figno, name)
% Applies correlation to the columns of D1 and D2 using  corr() function of
% matlab and visualizes the results
subN = size(D2,2)
SHOW_AVG_OVER_HUMANS = false;   % if true it will also make plot for average over human data
colcnt1 = size(D1,2); 
colcnt2 =size(D2,2); 
[C, pval] = corr(D1,D2, 'Type', 'Pearson');  % 'Pearson' %'Kendall''Spearman'
hum_not_correlated_with_any_sim = sum(pval>=0.05);
hh= hum_not_correlated_with_any_sim;
figure(figno); clf;
subplot(2,1,1)

sigix     = pval<0.05;
notsigix = pval>=0.05;

C(notsigix)=nan;

maxC = max(C); minC = min(C); [meanC, stdC] = mymean(C); 
%maxC = max(C); minC = min(C); meanC = mean(C); stdC = std(C);
maxbar = bar(maxC, 'FaceColor',[0.75,0.75,0.75]); hold on;
bar(meanC,'FaceColor',[1,0.8,0.8]); 
bar(minC,'FaceColor',[0.5,0.5,0.5]); 
er = errorbar(meanC,stdC,'LineWidth',1.2,'Color',[0,0,0], 'LineStyle', 'none');
legend('MAX','MEAN','MIN')
ylabel(sprintf('r values of significant correlations [p<0.05]'));
xlabel('Participant no');

%[maxv, maxi] = max(meanC); [minv, mini] = min(meanC);
title(sprintf('Significant [p<0.05] correlations of the %s of the participants with model generated %ss',name,name))
%title(sprintf('Mean correlation with simulation: MIN=%1.3f (part.no:%d); MAX=%1.3f (part.no:%d) [all mean:%1.3f]',minv, mini,maxv,maxi, mean(meanC)));
for k = 1:length(maxbar.XData)
    tt=text(maxbar.XData(k)-0.2, maxbar.YData(k)+0.05,sprintf('(%d)',hh(k)) );
end
xticks(1:subN)
axis([0.25,subN+0.75,-0.1,1.15]);
fprintf('Number of pairs failing to be significant (%s):%d/%d (%2.2f%%)\n',name,sum(hh),colcnt1*colcnt2, 100*sum(hh)/(colcnt1*colcnt2));

subplot(2,1,2)
PLOT_PVAL=1;
if PLOT_PVAL==1
    avg_pval = mean(pval); std_pval = std(pval);
    bar(avg_pval,'LineWidth',2); hold on
    er = errorbar(avg_pval,std_pval,'LineWidth',1.2,'Color',[0,0,0], 'LineStyle', 'none');
    %plot(avg_pval,'LineWidth',2); hold on; % The correlation pvalues averaged over sims
    %%plot(mean(pval,2), 'LineWidth',2); hold on; % The correlation pvalues averaged over subjects
    plot([0,20],[0.05,0.05],'--');
    xlabel('Participant no');
    ylabel(sprintf('p-mean for correlation of human and simulation %s',name))
    %ylabel(sprintf('Min, Mean, Max correlation (%s)', name))
    title(sprintf('Mean p-values averaged over sim.runs of the correlation and their standard deviation'));
else

end


fprintf(sprintf('The best matching simulation for a subject would have these correlation coefficients and p-values for %s:\n', name))
[vv,ii] = max(C);
[vv; pval(ii)];


figure(figno+10); clf;
SHOWONLY_MAX_MATCH = 0;
if SHOWONLY_MAX_MATCH==1
    [kval, kix] = max(vv);
    d2 = D2(:,kix);
    d1 = D1(:,ii(kix));
    plot([d1,d2],'Linewidth',2);
    title(sprintf('Maximum correlating pair. Human:%d, Sim:%d',kix, ii(kix)));
else
   MMIN = min(min(min(D2)), min(min(D1)));
   MMAX = max(max(max(D2)), max(max(D1)));
    for k=1:subN
        subplot(5,4,k)
        d2 = D2(:,k);
        d1 = D1(:,ii(k));
        r = C(ii(k),k); p = pval(ii(k),k);
        plot([d1,d2],'Linewidth',2); 
        title(sprintf('Participant %d : best matching Simulated %s (sim %d)',k,name,ii(k)));
        ylabel(name);xlabel('trial no');
        text(1,MMAX*0.9,sprintf('r=%1.3f p=%1.5f',r,p))
        axis([0,84,MMIN,MMAX])
    end
end


if (SHOW_AVG_OVER_HUMANS)
    %**This time let's average the correlations over the simulation runs ******
    figure(figno+1); clf;
    subplot(2,1,1)
    maxC = max(C,[],2); minC = min(C,[],2); meanC = mean(C,2); stdC = std(C,0,2);
    bar(maxC,'FaceColor',[0.75,0.75,0.75]); hold on;
    bar(meanC,'FaceColor',[1,0.8,0.8]); hold on;
    bar(minC,'FaceColor',[0.5,0.5,0.5]); hold on;
    er = errorbar(meanC,stdC,'LineWidth',1.2,'Color',[0,0,0], 'LineStyle', 'none');
    ylabel(sprintf('Min, Mean, Max correlation averaged over participant (%s)', name))
    xlabel('Model simulation runs');
    [maxv, maxi] = max(meanC); [minv, mini] = min(meanC);
    title(sprintf('Mean correlation with participants: MIN=%1.3f (sim.no:%d); MAX=%1.3f (sim.no:%d)',minv, mini,maxv,maxi));
    subplot(2,1,2)
    avg_pval = mean(pval,2); std_pval = std(pval,0,2);
    bar(avg_pval,'LineWidth',2); hold on
    er = errorbar(avg_pval,std_pval,'LineWidth',1.2,'Color',[0,0,0], 'LineStyle', 'none');
    %plot(avg_pval,'LineWidth',2); hold on; % The correlation pvalues averaged over sims
    %%plot(mean(pval,2), 'LineWidth',2); hold on; % The correlation pvalues averaged over subjects
    plot([0,60],[0.05,0.05],'--');
    xlabel('Simulation run');
    ylabel(sprintf('p-mean for correlation of human and simulation %s',name))
    %ylabel(sprintf('Min, Mean, Max correlation (%s)', name))
    title(sprintf('Mean p-values avereged over participants of the correlation and their standard deviation'));
    
end
    %Except participant 7 all the particapants TA shows significant correlation.
    %Only 26/1180 (2.2%) pairs of (subject,simrun) pairs have low correlation  (<0.21) with p>0.05.
    %16 of these come from participant 7. Removing him/her makes the pairs
    %failing to show significant correlation to less then 1%  (0.8475%)
    
    %In ITA the number of pairs failing to show significant correlations is 
    %340/1180 (28.8%), which is higher then TA. So our model better captured
    %the TA patterns of participants.

display('returning')
return

binCNT=16;
figure(figno); clf; 
subplot(2,1,1); hist(hTA,binCNT); title('subject correlation histogram (TA) ') 
subplot(2,1,2); hist(hITA,binCNT); title('subject correlation histogram (ITA) ') 


I = eye(subN);
avg_corrTA  =  sum(sum(C  - I))/(subN*subN-subN);
avg_corrITA =  sum(sum(corrITA - I))/(subN*subN-subN);

[vals,max_m]=max(C-I);
[maxcor, max_n] =max(vals);
max_m = max_m(max_n);
fprintf('Avg. TA correlation is %1.4f \n',avg_corrTA)
fprintf('Max TA correlation is %1.4f between subject %d and %d\n\n',maxcor,max_m,max_n);

[vals,max_m]=max(corrITA-I);
[maxcor, max_n] =max(vals);
max_m = max_m(max_n);
fprintf('Avg. ITA correlation is %1.4f \n',avg_corrITA)
fprintf('Max ITA correlation is %1.4f between subject %d and %d\n\n',maxcor,max_m,max_n);


%---
%I = eye(84);
%avg_corrTA  =  sum(sum(C  - I))/(subN*subN-subN);
%avg_corrITA =  sum(sum(corrITA - I))/(subN*subN-subN);
%bar(max(C,[],2),'FaceColor',[0.75,0.75,0.75]); hold on;
%bar(mean(C,2),'FaceColor',[1,0.8,0.8]); hold on;
%bar(min(C,[],2),'FaceColor',[0.5,0.5,0.5]); hold on;


function [avg, stdev] = mymean(A)
% Takes the column mean and std of a matrix by ignoring the nans
avg = zeros(1,size(A,2));
stdev = avg;

for k = 1: size(A,2)
    v = A(:,k);
    v(isnan(v))=[];
    avg(k)   = mean(v);
    stdev(k) = std(v);
end
x=1