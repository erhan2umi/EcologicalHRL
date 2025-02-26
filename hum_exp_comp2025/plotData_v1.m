

load('dataTA')
figure(201); clf;

%subL = 1:20; m = 4; n = 5;
subL = 5:8; m = 4; n= 1;
N = length(subL);

for sub = subL
    subplot(m,n,sub-subL(1)+1);
    plot(1:84, dataAllNormTrials(sub,:), '-'); hold on
    for k = 1:84,
        if ~isnan(dataFailNormTrialsTA(sub,k))
            plot([k k], [dataFailNormTrialsTA(sub,k)-50 dataFailNormTrialsTA(sub,k)+50], '-r'); 
        end
    end
    plot(1:84, dataSuccNormTrialsTA(sub,:)-5, '.g','Markersize',16); hold on
    title(sprintf('sub%d TA',sub));
        
end

load('dataITA')

figure(202); clf

for sub = subL
    subplot(m,n,sub-subL(1)+1);
    plot(1:84, dataAllNormTrialsITA(sub,:), '-'); hold on
    for k = 1:84,
        if ~isnan(dataFailNormTrialsITA(sub,k))
            plot([k k], [dataFailNormTrialsITA(sub,k)-5 dataFailNormTrialsITA(sub,k)+5], '-r'); 
        end
    end
    plot(1:84, dataSuccNormTrialsITA(sub,:), '.g','Markersize',16); hold on
    title(sprintf('sub%d ITA',sub));
        
    
end


    %scatter(1:84, dataAllNormTrialsITA(sub,:), '.k'); hold on
    %scatter(1:84, dataFailNormTrialsITA(sub,:), '.r'); hold on
    %scatter(1:84, dataSuccNormTrialsITA(sub,:), '.g'); hold on