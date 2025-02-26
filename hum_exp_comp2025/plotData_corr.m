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

ta=[];
ita = [];
for s1=1:20
      ta(:,s1)  = dataAllNormTrials(s1,:);

      ita(:,s1) = dataAllNormTrialsITA(s1,:);
      nanmark = isnan(ta(:,s1));
      nanix = find(nanmark);
      if ~isempty(nanix)
          oldta = ta(:,s1);
          fprintf('sub %d data has %d NANs, fixing them by overwriting with neighbor values\n',s1,length(nanix));
      else
          oldta = [];
      end
      while ~isempty(nanix)
        if nanix(1)~=1   % if the first NAN is not at index 1
            if ~isnan(ta(nanix(1)-1))   % if left is ok copy it
                 ta(nanix(1),s1)  =  ta(nanix(1)-1,s1);
                ita(nanix(1),s1)  = ita(nanix(1)-1,s1);
            else
                fprintf('********* Cannot be!\n')
            end
        else             % if first nan is at index 1
            goodix = find(~nanmark);      % replace it with the first nonNan
             ta(1,s1)  =  ta(goodix(1),s1);
            ita(1,s1)  = ita(goodix(1),s1);
        end
        nanmark = isnan(ta(:,s1));  % repeat the process % at each iteration 1 NAN is overwritten by a neighbour
        nanix = find(nanmark);
      end

      %Enable this to see the fixed NANs
      %if ~isempty(oldta)
      %    nansfixed = [oldta, ta(:,s1)]
      %    display("OK?")
      %end
end


selfcorr(ita, ta, 5);