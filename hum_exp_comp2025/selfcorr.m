function selfcorr(ita, ta, figno)
subN = size(ta,2);
corrTA = corr(ta);
corrITA = corr(ita);

hTA = reshape(triu(corrTA,1),1,subN^2);    hTA(hTA==0)   = []; % take the upper triangle and remove 1s and 0s
hITA = reshape(triu(corrITA,1),1,subN^2); hITA(hITA==0)  = []; 

binCNT=16;
figure(figno); clf; 
subplot(2,1,1); hist(hTA,binCNT); title('subject correlation histogram (TA) ') 
subplot(2,1,2); hist(hITA,binCNT); title('subject correlation histogram (ITA) ') 


I = eye(subN);
avg_corrTA  =  sum(sum(corrTA  - I))/(subN*subN-subN)
avg_corrITA =  sum(sum(corrITA - I))/(subN*subN-subN)

[vals,max_m]=max(corrTA-I);
[maxcor, max_n] =max(vals);
max_m = max_m(max_n);
fprintf('Avg. TA self-correlation is %1.4f \n',avg_corrTA)
fprintf('Max TA self-correlation is %1.4f between subject %d and %d\n\n',maxcor,max_m,max_n);
[vals,max_m]=max(corrITA-I);
[maxcor, max_n] =max(vals);
max_m = max_m(max_n);
fprintf('Avg. ITA self-correlation is %1.4f \n',avg_corrITA)
fprintf('Max ITA self-correlation is %1.4f between subject %d and %d\n\n',maxcor,max_m,max_n);
