load('dataTA')

figure(10); clf;
for sub = 1:20
    
    subplot(4,5,sub)
    scatter(1:84, dataAllNormTrials(sub,:), '.k'); hold on
    scatter(1:84, dataFailNormTrialsTA(sub,:), '.r'); hold on
    scatter(1:84, dataSuccNormTrialsTA(sub,:), '.g'); hold on
        
end

load('dataITA')

figure(11); clf;

for sub = 1:20
    
    subplot(4,5,sub)
    scatter(1:84, dataAllNormTrialsITA(sub,:), '.k'); hold on
    scatter(1:84, dataFailNormTrialsITA(sub,:), '.r'); hold on
    scatter(1:84, dataSuccNormTrialsITA(sub,:), '.g'); hold on
        
end
