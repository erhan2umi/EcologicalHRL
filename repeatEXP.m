% Erhan Oztop
% Februray 2025: 
% This scripts calls HIREL_model.m to make a virtual experiment over "repeats" many
% times. Each repeat is considered analogous to an experiment by a human
% participant.
% All the virtual participant data is stored in a folder with a name
% constructed from the date and time of the execution.
% The folder is then analysed by show_stats_disc_succnorm_sem.m and some
% plots are generated.


clear;
repeats = 10; % Number of experiments to run (# virtual subjects)
c=clock;
timestamp=sprintf('%4d_%d_%d_%d_%d',c(1),c(2),c(3),c(4),c(5));

foldername = sprintf('%s%s','./REPS',timestamp);
filenamebase = 'expsim';


eval(sprintf('mkdir %s',foldername));
fprintf('Saving the results to folder %s...\n',foldername);
for k=1:repeats
    makeSILENT = true;  % this is to suppress drawing in HIREL_model_v2
    HIREL_model 
    command = sprintf('save %s/%s_rep%d exper',foldername, filenamebase,k);
    fprintf('[%s] is being executed...\n',command);
    eval(command);
end
com = sprintf('show_stats_disc_succnorm_sem(''%s'',%d);', foldername, repeats);  
fprintf('Executing %s...\n',com);
eval(com)
fprintf('** The experiment data is saved to folder %s **\n',foldername);