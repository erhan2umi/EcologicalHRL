function [SFS_delTA,  SFS_delITA, SFS_delcost, SFS_delix] = check_SFS(exper)
    COMP=1; SAFETY=2; FALL=3;
    SFS_delTA  = []; SFS_delix  = [];
    SFS_delITA = []; 
    SFS_delcost = [];
    pert_ix = exper.PERTix;   %
    pert_six= exper.succ_pert_ix;
    pert_fix= exper.fail_pert_ix;


    f = exper.h_failed(pert_ix);
    k0 = 1;
    for k = 1:length(f)
        if k<-100, continue, end
        if f(k)==0, continue, end % a success trial ignore
        pri_s = find(0==f(k0:k-1))+k0-1;
        pos_s = find(0==f(k+1:end))+k+1-1;
        if (isempty(pri_s) || isempty(pos_s)), continue, end
        TAchange  = exper.h_TA(pert_ix(pos_s(1))) - exper.h_TA(pert_ix(pri_s(end)));
        ITAchange = exper.h_ITA(pert_ix(pos_s(1))) - exper.h_ITA(pert_ix(pri_s(end)));
        costchange = exper.h_costs(pert_ix(pos_s(1)),COMP) - exper.h_costs(pert_ix(pri_s(1)),COMP);
        k0 = pos_s(1);  % make sure to count F's as one in such cases SFFFFS  
        SFS_delTA = [SFS_delTA, TAchange];
        SFS_delITA = [SFS_delITA, ITAchange];
        SFS_delcost = [SFS_delcost, costchange];
        %SFS_delix = [SFS_delix; [pert_ix(pri_s(end)) pert_ix(pos_s(1))]];
        %delix = [delix, pert_ix(k)];
        %SFS_delix = [SFS_delix, pert_ix(k)];
        %SFS_delix = [SFS_delix; [pert_ix(pri_s(end)) pert_ix(k)]];
        SFS_delix = [SFS_delix; [pert_ix(pri_s(end)) pert_ix(pos_s(1))]];
    end
end % of check_SFS

 