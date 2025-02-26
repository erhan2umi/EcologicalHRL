function [SS_delTA, SS_delITA, SS_delcost, SS_delix] = check_SS(exper)
    COMP=1; SAFETY=2; FALL=3;
    SS_delTA  = []; SS_delix  = [];
    SS_delITA = []; 
    SS_delcost =[];
    pert_ix = exper.PERTix;   %
    pert_six= exper.succ_pert_ix;
    pert_fix= exper.fail_pert_ix;


    f = exper.h_failed(pert_ix);
    for k = 1:length(f)-1
        if f(k)==1, continue, end % a fail trial ignore
        if f(k+1)==0 %we have now SS case
            TAchange  = exper.h_TA(pert_ix(k+1)) - exper.h_TA(pert_ix(k));
            ITAchange = exper.h_ITA(pert_ix(k+1)) - exper.h_ITA(pert_ix(k));
            costchange = exper.h_costs(pert_ix(k+1),COMP) - exper.h_costs(pert_ix(k),COMP);
            SS_delTA = [SS_delTA, TAchange];
            SS_delITA = [SS_delITA, ITAchange];
            SS_delcost = [SS_delcost, costchange];
            SS_delix = [SS_delix; [pert_ix(k) pert_ix(k+1)]];
        end
    end
end % of check_SS