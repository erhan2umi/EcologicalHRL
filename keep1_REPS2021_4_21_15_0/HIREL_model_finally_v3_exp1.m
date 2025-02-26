% ERHAN memo
% Feb 2,2021: -HIREL_model_finally_V1 seems to settled many things and does an OK job (tested on Perturbed trials) Use it with repeat_exp_v1.m 
% Top level does via point search via simple RL
% (stateless). In parallel idyn learning takes place.
% This model works fine: noperturbation + adaptivepart when the adaptive part is based on pseudoinverse. 
% To achieve this use these settings:
% EXE.learning_type = 0;           % 2: matlab NN,  0: average batch pseuinv, 1:my grad desc. has problems
% EXE.offset_dyn = 1;              % 1:learn dyn as unperturbed + novel dyn   0: learn the usual way
% fb_gainP = fb_gainD = [0,0]


% NOTE: Nov 26 2020: HIREL_model_v2.m started. Will split safety and cost
% optimization as 2 level hirerchical RL

global body SPpar
global W first_roll trialc
global EXE
global perfectIDYN 
global fb_gainP_base fb_gainD_base fb_MUL fb_delay_ms
global RBUF RBUFc RBUFfull  
global SILENT MAKEMOVIE
global myfminpars;  
global cart_grad;
cart_grad=[];
SILENT = false;
MAKEMOVIE = false; %true;
movieFrames = [];
frame_no = 0;
movieRate = 4;
if exist('makeSILENT','var')
    SILENT = true; 
end
clear makeSILENT
exper = [];

MAX_fb_MUL = 2000;
fb_MUL = 1000;
fb_delay_ms = 50;        % delay in milliseconds
fb_gainP_base =  [2, 3];     
fb_gainD_base =  [2, 10]/40;     

perfectIDYN = 0;

createOtherGlobals;

x0 = [body.SQUAT_STATE(1),0,body.SQUAT_STATE(2),0];     % initial state

cx0 = body.SQUAT_CART;
cx1 = body.EXTENDED_CART;
cxm = 0.5*(cx0+cx1);

xlin = linspace(cx0(1),cx1(1),2+SPpar.NumCPm);
ylin = linspace(cx0(2),cx1(2),2+SPpar.NumCPm);
initialTraj_cart_CPm = [xlin(2:end-1);ylin(2:end-1)];  % via points spread linearly between start and end 

%initialTraj_cart_CPm(1,2)=initialTraj_cart_CPm(1,2)+0.02;
if (SPpar.NumCPm==3)
    %initialTraj_cart_CPm = [0.00+cx0(1)  0.00+cxm(1) 0.01+cx1(1) ;     % Cart coordinates for the control point
    %                        0.09+cx0(2)  0.00+cxm(2) -0.0+cx1(2)];
    initialTraj_cart_CPm = [0.0+cx0(1)  0.0+cxm(1) 0.0+cx1(1) ;     % Cart coordinates for the control point
                            0.0+cx0(2)  0.00+cxm(2) -0.0+cx1(2)];
else
    
    initialTraj_cart_CPm = [cx0(1)    cx0(1) ;     % Cart coordinates for the control point
                             0.9    1.2];
end

                    
if (SPpar.NumCPm==3)
    initialTraj_cart_CPm_fb2 = [0.1    0.15  0.1 ;     % We use this as additional data to leanr unperturbed dynamics.
                               0.5    .9    1.2];
end
if (SPpar.NumCPm==2)
initialTraj_cart_CPm_fb2 = [0.1      0.1 ;     % We use this as additional data to leanr unperturbed dynamics.
                           0.5       1.2];
end


trajSP = make_trajSP(initialTraj_cart_CPm, EXE.endtime);
trajSP_fb1 = trajSP;
trajSP_fb2 = make_trajSP(initialTraj_cart_CPm_fb2, EXE.endtime);
N = SPpar.NumCPm*2;
mu_init = [reshape(initialTraj_cart_CPm',N,1)];             % state are the knot configuration (CPm)                                      

  
EXE.falltype = {'B','F'};
EXE.fallstat_txt = {'xF','xB','xS','N/A','S','B','F'}; %-3 -2 -1 0 1 2 3
EXE.fallstat_SUCC = 1;
EXE.fallstat_BFALL = 2;
EXE.fallstat_FFALL = 3;
EXE.fallstat_xSUCC = -1;  %should not happen 
EXE.fallstat_xBFALL = -2; % falls due to xnoiseon
EXE.fallstat_xFFALL = -3;
EXE.compcost_rho = 1; %0.25;                 % avg. decay par for baseline cost calculation
EXE.fallcost_rho = 1;             % avg. decay par for baseline fallcost calculation
EXE.safetycost_rho = 1; %0.5;
EXE.COMP=1; EXE.SAFETY=2; EXE.FALL=3; EXE.FALL_BA=4; EXE.FALL_FR=5;  % constants for indexing the cost vector
EXE.noCONTafterFALL = 1;  % Warning: check if this is manipulated for DEAD and CATCH trials
perturbed = body.par;
unperturbed = [body.par(1), 0];
data=[];

 

% To form the initial internal model. We need free standup movement with no perturbation. 
% However, one roll is not sufficient for this.
% pseudo_inv does not like it. Need to make more runs with noise injection to learn it properly.
% This could be learned and saved as a file and loaded at the start time.
% EXE.offset_dyn=1 enable to learn the novel dynamics as
% unperturbed+novel_dynamics. This works fine with pseudo-inv
% (EXE.learning_type=0) but matlab NN does not work well with it. (probably stopping condition is reached too quick) 
% Anyways:  
%  For NN (EXE.learning_type=2) use EXE.offset_dyn=0
%  For Pseudo-inv (EXE.learning_type=2) use EXE.offset_dyn=1
fprintf('*) Making an UNPERTerbued run with no feedback \n');
[compcost, safetycost, fallcost, ~, exe_data] = execute_movement(x0, trajSP,0, EXE.timesteps, EXE.endtime, unperturbed, 114 );  
exe_data.W = {W};
title('Unperturbed run with true dynamics and NO feedback','FontSize', 8);
set(gcf,'Position',[0,1024,300,350])
first_roll = 1;  % This forces the system to use exact idyn. This is necessary to get the next plot.
fprintf('*) Making an UNPERTerbued run with true dynamics & feedback\n');
[compcost_fb1, safety_cost_fb1, fallcost_fb1, ~, exe_data_fb1] = execute_movement(x0, trajSP_fb1,fb_MUL, EXE.timesteps, EXE.endtime, unperturbed, 115 );  
title('Unperturbed run with true dynamics and WITH feedback (fb1)','FontSize', 8);
set(gcf,'Position',[0,880,300,350])
exe_data_fb1.W = {W}; exe_data_fb1.isfalltrial = 0;  % assume it is success (normally it should be!)
data = merge_data(data, exe_data_fb1);

fprintf('*) Making an UNPERTerbued run with learned dynamics & feedback\n');
[compcost_fb2, safetycost_fb2, fallcost_fb2, ~,  exe_data_fb2] = execute_movement(x0, trajSP_fb2,fb_MUL, EXE.timesteps, EXE.endtime, unperturbed, 1115 );
title('Unperturbed run with trueXLearned? dynamics and WITH feedback (fb2)','FontSize', 8);
set(gcf,'Position',[330,880,300,350])
exe_data_fb2.W = {W}; exe_data_fb2.isfalltrial = 0;  % assume it is success (normally it should be!)
data = merge_data(data, exe_data_fb2);
 
[W,err] = est_idyn(data, EXE.learning_type, unperturbed);   % learn w/ two roll outs only. (fb1 fb2)
maxv = max(abs(exe_data.net_u));
per_err = 100*[err(:,1)/maxv(1), err(:,2)/maxv(2)];
if ~SILENT
    figure(120); 
    clf;
    plot(per_err,'Linewidth',2); 
    legend('force','torque'); title('percent error of the the initial dynamic model');
    %[cost, fallcost, ~, exe_data] = execute_movement(x0, trajSP, EXE.timesteps, EXE.endtime, unperturbed, 1);  % 0*par(2) is the realistic case
    %title('Unperturbed run with learned unperturbed dynamics.');
end
ccpm = initialTraj_cart_CPm;
last_cpert = 0;
errL = ones(EXE.numtrial,1);

%-------------------- HRL ----------------------
mu = mu_init;             % state are the knot configuration (CPm)
a = mu;                   % actions on the knots


success_cnt_reached = 0;
it = 0;  
impc = 0;
MAXIT = 320;    % should be this same as EXE.numtrial?  EXE.numtrial is used for experience buffer size setting. Better check it.
%[baseline_cost, fallcost0, zmp, xy, Th]  = objfun(state2CPm(mu),CPb,CPe,kb,q1s,q2s,t0,t1,par);
trajSP = make_trajSP(initialTraj_cart_CPm, EXE.endtime);
[baseline_compcost_pert, safetycost_pert, fallcost_pert, ~, exe_data] = execute_movement(x0, trajSP, 0, EXE.timesteps, EXE.endtime, perturbed, 116);  
title('perturbed run with learned unpert.dyn and NO feedback','FontSize', 8);
set(gcf,'Position',[0,440,300,350])
[baseline_compcost, baseline_safetycost, baseline_fallcost, ~,  exe_data_fb] = execute_movement(x0, trajSP, fb_MUL, EXE.timesteps, EXE.endtime, perturbed, 117);  
title('perturbed run with learned unpert.dyn and WITH feedback','FontSize', 8);
ss = sprintf('[INITIAL compcost:%2.3f safetycost:%2.3f  fallcost:%2.3f',baseline_compcost, baseline_safetycost,baseline_fallcost);
fprintf('%s\n',ss);
ylabel(ss);
set(gcf,'Position',[0,0,300,350])
[x, ddq1_, ddq2_] = state_from_spline(state2CPm(mu),trajSP.CPb,trajSP.CPe,SPpar.kb,trajSP.q1s,trajSP.q2s,0);

gooda_safety_cost = 1e10;
compcost_adv = -1e-10;
delta = 0;
started = 1;
good_perta = mu;
lastunpertmu = mu;   


cond_unPERT = 1;
cond_PERT = 2;
cond_CATCH = 3;
cond_PERTx = 4;
cond_DEAD = 5;
condtxt = {'unPERT','  PERT',' CATCH',' PERTx','  DEAD'};
EXE.unPERT_cnt = 1  ; %15  ;
EXE.DEAD_cnt = 1 ; %15;
EXE.PERTx_cnt  = 1 % 15;
EXE.reqsucc   = 60;    % subject needs to complete this many succesfull trails in PERT condition



exper.condTrialCounts = [EXE.unPERT_cnt 0 1 EXE.PERTx_cnt  EXE.DEAD_cnt];   % PERT condition is unknown apriori
exper.h_sig = zeros(MAXIT,1);
exper.h_a = zeros(MAXIT,length(mu));
exper.h_fallcost= zeros(MAXIT,1);
exper.h_baseline_fallcost= zeros(MAXIT,1);
exper.h_adv= zeros(MAXIT,1);
exper.h_cost= zeros(MAXIT,1);
exper.h_baseline_cost= zeros(MAXIT,1);
exper.h_TA= zeros(MAXIT,1);
exper.h_ITA= zeros(MAXIT,1);
exper.h_fb_MUL = zeros(MAXIT,1);
exper.h_gradnorms = zeros(MAXIT,3);
exper.body = body;
exper.EXE = EXE;
exper.SPpar = SPpar;
exper.date = datestr(now);
exper.h_failed = [];
exper.h_ITA = [];
exper.catchtrial_no = -1;    
exper.lastPERTix = -1; 
exper.success_cnt_reached = false;
FALLREW = 0; SUCCREW = 1;
exper.lastPERTexe = [];
exper.lastunPERTexe = [];
exper.CATCHexe = [];
exper.firstDEADexe = [];
exper.lastDEADexe = [];

fallcost0 = 0;
failed = 1;
trialc = 0;
succ = [0,0,0,0,0];
good_pertmu = [];
base_sig = 0.02; 
lrate = 1 ;
COST_grad_mul = 5;   % for cost the lrate is multiplied with this
if (SPpar.NumCPm==2)
    COV=[1 0.5 0.5 0.5 ]';
else
    COV=[1 1 1 1 1 1]';
end

GRmax=1.0;

% --------------- START LEARNING TRIALS --------------------------------
for cond= [cond_unPERT cond_PERT cond_CATCH cond_PERTx cond_DEAD] % [ unperturbed(15),perturbed(succ=60),catch(1),unperturbed(15) ];
    cit = 0;         
    done = false;
    sig = base_sig;
    last_failed = 1;  % No success trial from a previous condition should effect the current condition 
    winmul = 900+cond;
    compcost_adv = 0;
    fprintf('STARTing CONDITION: %s...\n',condtxt{cond});
    if cond == cond_PERT | cond == cond_PERTx
        environ = perturbed;
        dangerCond = 1;
    else
        environ = unperturbed;
        dangerCond = 0;
        %sig = base_sig*0;
    end
    
    if cond == cond_PERT 
        unpertW = W;
        good_perta = [];
        good_pertmu = []; 
        body.SAFETY_SWITCH_succ = round(random('norm', body.SAFETY_SWITCH_succ_mu, body.SAFETY_SWITCH_succ_sig));
    end
    if cond == cond_PERTx  %ignore catch trials effect on mu
        mu = good_pertmu;
        if isempty(mu), mu=mu_init; end
        body.SAFETY_SWITCH_succ = round(random('norm', body.SAFETY_SWITCH_succ_mu, body.SAFETY_SWITCH_succ_sig));
    end

    if cond == cond_DEAD
        if (isempty(lastunpertmu))
            fprintf('Strange! No succesful unperturbed trial was achieved, using initial mu.\n');
            mu = [reshape(initialTraj_cart_CPm',N,1)];
        else
            mu = lastunpertmu;
        end
        unpertW_weight = 0.75;
        W = W*(1-unpertW_weight)+ unpertW_weight*unpertW;  % ERH
    end
    
    if cond == cond_CATCH
        dummy=1;
    end
    effort_try = 0;
    while (~done)
        trialc = trialc + 1;
        cit = cit + 1;
      
        if failed && dangerCond && ~isempty(good_pertmu)  % reset to previous success
            mu = good_perta;
        end

        bias = (failed && (cond == cond_PERT))*0.03*[1 1 1 0 0 0]';  %ERH was 0.025*
        aa = random('norm',mu+bias ,sig*COV); 
        aa(3)=aa(2)/1.5;
        aa(1)= aa(2)/1.5;
        % vv Apply constrains to the sampled action
        a = aa; 
        ina = a;
        if SPpar.NumCPm>2, 
            if a(2) < 0.5*(a(1)+a(3)), a(2) = 0.5*(a(1)+a(3))*1.01; end
            if a(6) < a(5), a(6) = a(5); end
            if a(4) > a(5), a(4) = a(5); end
        end
        violix = a < SPpar.slb_cart;  a(violix) = SPpar.slb_cart(violix); violc = sum(violix);
        violix = a > SPpar.sub_cart;  a(violix) = SPpar.sub_cart(violix); violc = violc + sum(violix);   
        if (violc>0)
            %in_outa = [ina a]
            %fprintf('violc:%d\n',violc);
        end
        % ^^ Apply constrains to the sampled action
        % vv get the spline corresponding to action a, an make rollout, get the costs and data for Idyn learning.
        
        delta = a - mu;
        
        
        trajSP = make_trajSP(a(1:N), EXE.endtime);  
        if SILENT 
            showwin = (mod(it,50)==0);
        else
            showwin = (mod(it,1 )==0);   % To see realtime mod(it,1)
        end
        
        %if cit>body.SAFETY_SWITCH_succ, body.xnoiseon = 1; end
        lastEXECW = W;
        %execom = a(1:N)'
        
        %EXE.noCONTafterFALL=1;
        %if (cond==cond_DEAD && cit==1) EXE.noCONTafterFALL=1; end
        %if (cond==cond_CATCH && cit==1) EXE.noCONTafterFALL=1; end
        [compcost, safetycost, fallcost0, fallstat, exe_data]= execute_movement(x0, trajSP, fb_MUL,EXE.timesteps, EXE.endtime, environ, showwin*winmul);
        ylabel(sprintf('Learning...(trial:%d)(%s)', it,condtxt{cond}));
        fallcost = abs(fallcost0);  % treat fall cost as usual.

        if ~SILENT && MAKEMOVIE && cond==cond_PERT
            frame_no = frame_no + 1;
            drawnow;
            movieFrames{frame_no} = getframe(gcf) ;
        end
        last_failed = failed;
        failed = (fallstat ~= EXE.fallstat_SUCC); 

        compcost_adv = (baseline_compcost - compcost);  
        fallcost_adv = (baseline_fallcost -fallcost);
        safetycost_adv = (baseline_safetycost - safetycost);

        baseline_compcost   = baseline_compcost*  (1-EXE.compcost_rho)+ EXE.compcost_rho*compcost;      
        baseline_fallcost   = baseline_fallcost*  (1-EXE.fallcost_rho)+ EXE.fallcost_rho*fallcost;
        baseline_safetycost = baseline_safetycost*(1-EXE.safetycost_rho)+ EXE.safetycost_rho*safetycost; 
        
        % Begin IDYN learning
        for kkk=1:1
        if (perfectIDYN==0 && cond~= cond_unPERT) % && cond~= cond_CATCH)  % ERH
            exe_data.W = {W};  % the data returned is obtained with this W (weights or network structure see EXE.learning_type)
            exe_data.fb_MUL = fb_MUL;   % this was the fb_MUL used
            exe_data.isfalltrial = exe_data.zmp*0+failed;
            data = merge_data(data, exe_data);

           [err, maxerr] = learn_batch(1e-8, EXE.learning_type, body.par);  % produces global W
            errL(trialc) = sum(err);
            idyn_err = sum(err);
            track_err = sum(mean(abs(data.del_q)));
            my_err = track_err;   % track_err or idyn_err
            err_rho = 0.5;
            if (cit==1)
                %err_baseline = idyn_err;
                err_baseline = track_err;
            else
                DE = (my_err - err_baseline);
                err_baseline = err_baseline*(1-err_rho)+ err_rho*my_err;
                %DE = 2*fell-1;
                fb_MUL = min(max(0,fb_MUL*(1+sign(DE)*0.1)), MAX_fb_MUL);
            end
            if fb_MUL == MAX_fb_MUL
                fprintf('fb_MUL saturated at %d trial:%d\n',it);
            end
            if cond==cond_DEAD
                fb_MUL = 2000; %MAX_fb_MUL;  %ERH
            end
            %fprintf('Roll %d done (FBMUL:%3.2f).\n',trialc,fb_MUL);
        else
            exe_data.W = [];
            exe_data.fb_MUL = 0;  
            exe_data.isfalltrial = exe_data.zmp*0+failed;
            data = merge_data(data, exe_data);
        end
        end
        % End IDYN learning


        if (mod(trialc,50)==1 || failed==22 || trialc==MAXIT)
           plot_data(data, body.par, 120,sprintf('roll:%d ',trialc));
           drawnow;
        end

  
        % store data for plotting etc. 
        exper.h_a(it+1,:) = a;
        exper.h_sig(it+1) = sig;
        exper.h_costs(it+1,:) = [compcost, safetycost, fallcost, exe_data.fallcost_ba, exe_data.fallcost_fr];
        exper.h_failed(it+1) = failed;
        exper.h_costadvs(it+1,:) = [compcost_adv, safetycost_adv, fallcost_adv];
        [exper.h_TA(it+1), exper.h_ITA(it+1)] = signed_patharea(exe_data);
        if (exper.h_TA(it+1)<0) && ~failed
            checkfallcost = exper.h_costs(it+1,3)
            fprintf('********* it %d, TA Area negative (%3.2f) and not failed!!\n', it+1, exper.h_TA(it+1));
            cod=1;
        end
        
        exper.h_baseline_costs(it+1,:) = [baseline_compcost, baseline_safetycost,  baseline_fallcost];
        exper.h_fb_MUL(it+1) = fb_MUL;
        %exper.h_gradnorms(it+1,:)=[norm(delcc), norm(delfc), norm(delsc)]; 
        LO = -1;
        HI =  1;
        compcost_advclip   = max([min([HI,compcost_adv]), LO]);
        fallcost_advclip   = max([min([HI,fallcost_adv]), LO]);   % if not fall second part should be almost zero
        safetycost_advclip = max([min([HI,safetycost_adv]), LO]);
        if ~SILENT
            %COMP_FALL_SAFETY_advclips = [fallcost_advclip fallcost_advclip safetycost_advclip ]*100
        end
        

        if last_failed || failed || cit < body.SAFETY_SWITCH_succ
            netadv =failed*fallcost_advclip +  dangerCond*(1-failed)*safetycost_advclip + 0.0*effort_try*safetycost_advclip;  % ERH   1*(1-failed) was 0.75*(1-failed)
            effort_try = 0;
        else
            if cond~=cond_unPERT, fprintf('ACTING with CCOST adv:%1.7f\n',compcost_advclip); end
            netadv = (COST_grad_mul*compcost_advclip + 0*fallcost_advclip + 0*safetycost_advclip);
            effort_try=1;
        end
        mu = mu + lrate*delta*netadv;
  
        
        
        if (failed==1)
           fprintf('%s: %d/%d> (%s) FAILED  (succ:%d) FBG:%3.1f comp[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] safety[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] fall[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] [sig:%1.4f]\n', ...
                condtxt{cond}, it+1, MAXIT, EXE.fallstat_txt{4+fallstat}, succ(cond), fb_MUL, compcost,baseline_compcost,compcost_adv, safetycost,baseline_safetycost,safetycost_adv, fallcost,baseline_fallcost,fallcost_adv, sig);
   
        else   % NO FALL 
            succ(cond) = succ(cond) + 1;
            fprintf('%s: %d/%d> EXE OK (succ:%d) FBG:%3.1f comp[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] safety[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] fall[cost,baseln,adv:%3.3f, %3.3f, %+3.3f] [sig:%1.4f]\n', ...
                condtxt{cond}, it+1, MAXIT,  succ(cond), fb_MUL,compcost,baseline_compcost,compcost_adv, safetycost,baseline_safetycost,safetycost_adv, fallcost,baseline_fallcost,fallcost_adv, sig);
            %a_mu = [a, mu]'

            if (cond==cond_PERT) 
                if((gooda_safety_cost*1 >= 1*safetycost) || isempty(good_perta))
                    good_perta = a;
                    good_pertmu = mu;
                    if isempty(exe_data.W)
                        good_pertW = exe_data.W;
                    else
                        good_pertW = exe_data.W{1};
                    end
                    good_fb_MUL= exe_data.fb_MUL;
                    good_pert_it = it+1;
                    gooda_safety_cost = safetycost;
                    %fprintf('*** TRIAL :%d> New safepoint  safety:%2.5f ccost:%2.5f [ITA:%2.5f. TA:%2.5f]\n',good_pert_it,gooda_safety_cost,compcost, exper.h_ITA(it+1), exper.h_TA(it+1));
                    %execom
                    %fprintf('*** \n');
                end
            end
                
            impc = impc + 1;

            if cond==cond_PERT && succ(cond) == EXE.reqsucc
                success_cnt_reached = 1;
                exper.lastPERTix = trialc;
            end
        end
        if norm(delta)<1e-5
            fprintf('Last run must be exact the same as the latest safepoint!\n');
        end
        
        it = it + 1; 
        sig = sig*0.997; %0.997 for short experiments
        
        if (it>=MAXIT) 
            done=true; 
        elseif cond==cond_unPERT
            done = (cit == EXE.unPERT_cnt);
            lastunpertmu = a; % made is a
            if done,  exper.lastunPERTexe = exe_data; end
        elseif cond==cond_PERT
            done = success_cnt_reached;
            if done
                exper.success_cnt_reached = success_cnt_reached;
                exper.lastPERTexe = exe_data; 
                exper.condTrialCounts(cond_PERT) = cit;
                W = good_pertW;
                exper.pert_IDYN_W = good_pertW; %exe_data.W{1};
                exper.pert_fb_MUL = good_fb_MUL; %exe_data.fb_MUL;
                exper.pert_gooda = good_perta;
                [costLastPert, fallcostLastPert, exe_dataLastPert]= execute_movement(x0, make_trajSP(good_perta, EXE.endtime), good_fb_MUL,EXE.timesteps, EXE.endtime, environ, 1*1001);
                ylabel(sprintf('Last cache safe...(cache trial:%d)(%s)', good_pert_it,condtxt{cond}));
                %[costLastPert, fallcostLastPert, exe_dataLastPert]= execute_movement(x0, trajSP, fb_MUL,EXE.timesteps, EXE.endtime, environ, 1*1001);
                %ylabel(sprintf('Last PERT...(trial:%d)(%s)', it,condtxt{cond}));
            end
        elseif cond==cond_CATCH
            done = (cit == 1);
            exper.catchtrial_no = trialc;
            if done, exper.CATCHexe = exe_data; end
        elseif cond==cond_PERTx
            done = (cit == EXE.PERTx_cnt);
            if done, exper.lastPERTxexe = exe_data; end
        elseif cond==cond_DEAD
            if cit == 1, exper.firstDEADexe = exe_data; end
            done = (cit == EXE.DEAD_cnt);
            if done, exper.lastDEADexe = exe_data; end
        else
            done = true;
            fprintf('ABOARTING! No such condition!\n');
        end
    end  % while(~done) learning_loop , i.e. end of experiment trials

    exper.lasttrialno = trialc;

end %for

if ( exper.success_cnt_reached)
    fprintf('\nExperiment summary: O [Required Perterbud Condition Success Count is Reached]\n');
else
    fprintf('\nExperiment summary: X [Required Perterbud Condition Success Count is NOT Reached in %d iterations]\n',MAXIT);
end
for kk = [cond_unPERT cond_PERT cond_CATCH cond_PERTx cond_DEAD]  
    fprintf('     COND   %s : %d/%d (SUCC/TOT num of trials) \n', condtxt{kk}, succ(kk), exper.condTrialCounts(kk));
end
if exper.catchtrial_no == -1
    exper.catchtrial_no = exper.lastPERTix; 
    exper.nocatchtrial=true; 
else
    exper.nocatchtrial=false; 
end
exper.unPERTix  = 1 : EXE.unPERT_cnt;
exper.PERTix    = EXE.unPERT_cnt+1 : exper.lastPERTix;  %exper.catchtrial_no-1;
exper.PERTxix   = exper.catchtrial_no+1 : EXE.PERTx_cnt+exper.catchtrial_no;
exper.DEADix = exper.PERTxix(end)+1: EXE.DEAD_cnt+exper.PERTxix(end);

if (exper.lastPERTix==-1)
    fprintf('****** CANNOT make sufficient succesfull stand ups! Will not plot data ... *********\n');
    return; 
end

if MAKEMOVIE && ~SILENT 
    frames2avi(movieFrames, 'pertmovie.avi', movieRate);
end

figure(222); clf;
%set(gcf,'Position',[320,430,700,930]);
subplot(4,1,1); cla;
%plot(exper.h_cost, 'Color', [0.6 0.6 0.6],'Linewidth',2);  hold on; 

%bar((exper.h_fallcost(1:it)~=0)*max([0.1;exper.h_baseline_cost(1:it)+exper.h_adv(1:it)]),'FaceColor',[1 0.7 0.7],'EdgeColor','none'); hold on

%bar((exper.h_failed(1:it))*max([0.1;exper.h_baseline_costs(1:it,1)+exper.h_costadvs(1:it,1)]),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none'); hold on
plot(exper.h_baseline_costs(1:it,EXE.COMP), 'Color', [0.2 0.2 0.2],'Linewidth',3);  hold on;
%plot(exper.h_baseline_costs(1:it,COMP)+exper.h_costadvs(1:it,COMP),'.-','Color',[0.2 0.2 1],'Linewidth',1);
plot(exper.h_baseline_costs(1:it,EXE.SAFETY), 'Color', [0.2 0.6 0.2],'Linewidth',3);  hold on;
%plot(exper.h_baseline_costs(1:it,SAFETY)+exper.h_costadvs(1:it,SAFETY),'.-','Color',[0.2 0.6 1],'Linewidth',1);
plot(exper.h_baseline_costs(1:it,EXE.FALL), 'Color', [0.8 0.4 0.2],'Linewidth',3);  hold on;
%plot(exper.h_baseline_costs(1:it,FALL)+exper.h_costadvs(1:it,FALL),'.-','Color',[0.8 0.4 1],'Linewidth',1);
%plot(exper.h_baseline_cost+exper.h_adv,'.','Color',[0 0 0]);
plot(exper.h_fb_MUL(1:it)/50,'-','Linewidth',2, 'Color',[1 0.1 0.1]);
%ylim([0*min(exper.h_baseline_costs(:,1)+exper.h_costadvs(:,1)),max(exper.h_baseline_costs(:,1)+exper.h_costadvs(:,1))]);
grid on;
%xlabel('trial no');
legend('compcost','safetycost','fallcost','feedbMUL/50')
%legend('falls','compcost','+','safetycost','+','fallcost','+','feedbMUL/120')

subplot(4,1,2); cla;
hhh = max([0.1,max(exper.h_costadvs(1:it,EXE.COMP))]);
bar(hhh*exper.h_failed(1:it),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none'); hold on
bar(-hhh*exper.h_failed(1:it),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');  

plot(exper.h_costadvs(1:it,EXE.COMP),'.-','Color',[0.2 0.2 0.2],'Linewidth',3);
plot(exper.h_costadvs(1:it,EXE.SAFETY),'.-','Color',[0.2 0.6 0.2],'Linewidth',3);
plot(exper.h_costadvs(1:it,EXE.FALL),'.-','Color',[0.8 0.4 0.2],'Linewidth',3);
ylim([-hhh,hhh]);
grid on;
%xlabel('trial no');
legend('falls','','compcost_{adv}','safetycost_{adv}','fallcost_{adv}')

subplot(4,1,3); cla;
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

subplot(4,1,4); cla;
tmp1 = max([exper.h_costs(1:it,EXE.COMP);exper.h_costs(1:it,EXE.COMP)+exper.h_costadvs(1:it,EXE.COMP)]);
bar(exper.h_failed(1:it)*max([0.02; 1.1*exper.h_costs(1:it,EXE.FALL);tmp1]),'FaceColor',[1 0.7 0.7],'EdgeColor','none'); hold on;
bar(exper.catchtrial_no, max([0.02; 1.1*exper.h_costs(1:it, EXE.FALL);tmp1]),0.5,'FaceColor',[1 0 0]);
bar((1-exper.h_failed(1:it))*max([0.02;1.1*exper.h_costs(1:it, EXE.FALL);tmp1]),'FaceColor',[0.8 1 1],'EdgeColor','none'); hold on;
bar((1-exper.h_failed(exper.unPERTix))*max([0.02;1.1*exper.h_costs(1:it, EXE.FALL);tmp1]),'FaceColor',[0.6 0.9 0.6],'EdgeColor','none'); hold on;
bar(exper.DEADix,(1-exper.h_failed(exper.DEADix))*max([0.02;1.1*exper.h_costs(1:it,EXE.FALL);tmp1]),'FaceColor',[0.6 0.9 0.6],'EdgeColor','none'); hold on;

bar(exper.h_costs(1:it,EXE.FALL),0.5 ,'FaceColor',[0.6 0.6 0.6]); hold on;
bar(exper.DEADix, exper.h_costs(exper.DEADix,EXE.FALL),0.5 ,'FaceColor',[0.9 0.9 0.9]); hold on;
bar(exper.unPERTix, exper.h_costs(exper.unPERTix,EXE.FALL),0.5 ,'FaceColor',[0.9 0.9 0.9]); hold on;
plot(exper.h_baseline_costs(1:it,EXE.FALL),'y-','Linewidth',2)
plot(exper.h_costs(1:it,EXE.COMP), 'Color', [0 0.2 0.9],'Linewidth',2); 

xlim([0,exper.lasttrialno])
legend('fall','catch','success','unper_0','unper_1','fallcost','fallcost_{unp}','fallcost_{unp}','baselnfallcost','ccost'); xlabel('trial no');

          
 % DRAW the FINAL performance of learning...
 trajSP = make_trajSP(good_perta, EXE.endtime);  
 W = exper.pert_IDYN_W; 
 
 [compcost, fallcost, exe_data]= execute_movement(x0, trajSP, exper.pert_fb_MUL, EXE.timesteps, EXE.endtime, perturbed, 1010);
 ylabel(sprintf('Final (it:%d) performance (a)', it));
 drawnow;

 
 trajSP = make_trajSP(good_pertmu, EXE.endtime);  
 [compcost, fallcost, exe_data]= execute_movement(x0, trajSP, exper.pert_fb_MUL, EXE.timesteps, EXE.endtime, perturbed, 1011);
 ylabel(sprintf('Final (it:%d) performance (mu)', it));
 drawnow;


 

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

analyzeSFS(exper, 123, EXE);  % plot in figure 123 and using cost as COMP (composite)
figure(1311); clf;hold on
ha=exper.h_a(exper.succ_pert_ix,:);
if (~SPpar.radial), 
    if SPpar.NumCPm>2
        scatter(ha(:,[3]), ha(:,[3+SPpar.NumCPm])); 
        scatter(ha(:,[1]), ha(:,[1+SPpar.NumCPm])); 
        scatter(ha(:,[2]), ha(:,[2+SPpar.NumCPm]));  
    end
    cc = exper.h_costs(exper.succ_pert_ix,EXE.COMP);
    [ccval, ccix] = min(cc);
    plot( ha(ccix,1:SPpar.NumCPm), ha(ccix,SPpar.NumCPm+1:SPpar.NumCPm*2), 'k-', 'Linewidth',2);
    text( mean(ha(ccix,1:2)), mean(ha(ccix,SPpar.NumCPm+[1:2])), sprintf('%1.3f',cc(ccix)));
    text( mean(ha(:,1)), mean(ha(:,1+SPpar.NumCPm))+0.02,sprintf('(%1.3f, %1.3f)', mean(ha(:,1)), mean(ha(:,1+SPpar.NumCPm))));
    text( mean(ha(:,2)), mean(ha(:,2+SPpar.NumCPm))+0.04,sprintf('(%1.3f, %1.3f)', mean(ha(:,2)), mean(ha(:,2+SPpar.NumCPm))));
    %text( mean(ha(:,3)), mean(ha(:,3+1+SPpar.NumCPm))-0.02,sprintf('(%1.3f, %1.3f)', mean(ha(:,3)), mean(ha(:,6))));
    grid minor;
    title(sprintf('Successful trials (%d/%d)',length(cc),length(exper.PERTix)));
    drawnow;
    [~, sortix] = sort(cc);  oldk = 0;

    if (SILENT || false),  kkEnd = -1; else, kkEnd = length(cc); end
    for kk = 1:kkEnd
        figure(1311);
        k = sortix(kk);
        if oldk~=0
            plot( ha(oldk,1:SPpar.NumCPm), ha(oldk,SPpar.NumCPm+1:SPpar.NumCPm*2), 'y-', 'Linewidth',1.2);   
            %plot(oldpos(:,1),oldpos(:,2),'c.');
        end
        oldk = k;
        plot( ha(k,1:SPpar.NumCPm), ha(k,SPpar.NumCPm+1:SPpar.NumCPm*2), 'r-', 'Linewidth',1.2); hold on;
        %temptrajSP = make_trajSP(ha(k,:)', EXE.endtime);
        %[pos, state] = trajSP2cart(temptrajSP);
        %oldpos = pos;
        %plot(pos(:,1),pos(:,2),'m.');  
        %text( mean(ha(k,1:2)), mean(ha(k,4:5))+k*0.005, sprintf('%1.3f',cc(k)));
        xs = sprintf('%dth compcost:%2.5f\n',kk,cc(k));
        fprintf(xs);
        xlabel(xs);
        %drawnow
        c=' ';
        %c=input('press enter enter '' '' to quit'); 
        if c==' ', break; end
    end %for k
else
    scatter(ha(:,[1]), ha(:,[3])); 
    scatter(ha(:,[2]), ha(:,[4]));    
end

figure(1312); clf;
subplot(2,1,1);
plot(exper.h_costs(succ_pert_ix(sortix),EXE.COMP), exper.h_TA(succ_pert_ix(sortix)),'b-','Linewidth',2); hold on;
plot(exper.h_costs(succ_pert_ix(sortix),EXE.COMP), exper.h_TA(succ_pert_ix(sortix)),'o'); hold on;
grid minor
xlabel('compcost'); ylabel('TA');
subplot(2,1,2);
plot(exper.h_costs(succ_pert_ix(sortix),EXE.COMP), exper.h_ITA(succ_pert_ix(sortix)),'b-','Linewidth',2); hold on;
plot(exper.h_costs(succ_pert_ix(sortix),EXE.COMP), exper.h_ITA(succ_pert_ix(sortix)),'o'); hold on;
grid minor
xlabel('compcost'); ylabel('ITA');


show_last_trails_v1(exper);

%-- MAIN finished here -------------

 

function createOtherGlobals    
    global EXE
    global SPpar
    global body
    global RBUF RBUFc RBUFfull 
    global W dynnet first_roll
    
    body.SQUAT_STATE = [0.5 ; 0.95*pi/2]; 
    body.EXTENDED_STATE = [1.2 ; 1.02*pi/2]; %1.05*pi/2
    body.SQUAT_CART = forward_kin(body.SQUAT_STATE' );
    body.EXTENDED_CART = forward_kin(body.EXTENDED_STATE');
    body.SAFETY_SWITCH_succ_mu  = 50; 
    body.SAFETY_SWITCH_succ_sig = 0; 
    body.SAFETY_SWITCH_succ_sig =  body.SAFETY_SWITCH_succ_mu;
    body.cost_type     = 0;   % -1:effort 0: composite  1:safety 
    body.effortW       = 1;   % composite cost weight 1: full effort,  0:full safety  
             
    body.xnoiseon      = 0;        
    body.mass = 80;
    body.pert=3*body.mass; %perturbation coefficient
    body.par = [body.mass, body.pert];
    body.heel_dist =  0.13;
    body.toe_dist = 0.30;
    body.fallcost_threshold =  0.01 ; % (0.001)*100;  % (x) x is average zmp violation in meters
    
    SPpar = [];
    SPpar.k = 6;      % spline order (4 for cubic, 6 for quintic B-spline)
    SPpar.kb = 3;        % spline order at the boundaries (1-pos, 2-pos,vel, 3-pos,vel,acc)
    SPpar.radial = 0;
    SPpar.MAXRAD = 0.5;   % m
    SPpar.MINRAD = 0.05;   % m

    SPpar.NumCPm = 3;    % number of control points to be optimized

    SPpar.lb_cart=repmat([0.03    ; 0.5],1,SPpar.NumCPm);     % [minx miny]
    SPpar.ub_cart=repmat([0.22  ; 1.2],1,SPpar.NumCPm);    % [maxx maxy]
    SPpar.slb_cart = [reshape(SPpar.lb_cart',2*SPpar.NumCPm,1)];   % limits flattened
    SPpar.sub_cart = [reshape(SPpar.ub_cart',2*SPpar.NumCPm,1)];


    EXE.timesteps = 60;   % use 100 for more accurate sim
    EXE.endtime   = 1.2;
    EXE.numtrial  = 25;   % is this used??
        
    EXE.RBUFsize = EXE.timesteps*10;    % sample buffer for dynamic learning, if too small early samples will be lost. EXE.endtime*EXE.timesteps*EXE.numtrial samples are generated in total
    EXE.RBUFbatchsize = EXE.timesteps*4;% batch size from here.

    %EXE.safety_tho=1;

    testx = [0,0,0,0]; 
    hx = expandX(testx,[0,0]);
    EXE.INPDIM = length(hx);
    EXE.OUTDIM = 2;
    RBUF.input = zeros(EXE.RBUFsize, EXE.INPDIM);
    RBUF.output = zeros(EXE.RBUFsize, EXE.OUTDIM);
    RBUFc = 0;
    RBUFfull = 0;
    W =zeros(EXE.INPDIM, EXE.OUTDIM);
    dynnet = feedforwardnet(15); % hid unit size of the network to use for learning 2
    EXE.learning_type = 0;           % 2: matlab NN,  0: average batch pseuinv, 1:my grad desc. has problems
    EXE.offset_dyn = 1;              % 1:learn dyn as unperturbed + novel dyn   0: learn the usual way
    EXE.offset_SC  = 1*(1-EXE.offset_dyn)+EXE.offset_dyn*2000;           % hack to make matlab nn to progress in grad.descent 
                                                             % (but does not work well, if nn then use EXE.offset_dyn=0)
                                                             
    first_roll = 1;         % this is not an option. Must be always 1
end



function outdata = merge_data(data, exedata)
    global EXE RBUF RBUFc RBUFfull  
    global W
    if isempty(data)
        outdata=exedata;
    else
        outdata.des_xtraj   = [data.des_xtraj; exedata.des_xtraj];
        outdata.act_xtraj   = [data.act_xtraj; exedata.act_xtraj];
        outdata.des_cartpos = [data.des_cartpos; exedata.des_cartpos];
        outdata.act_cartpos = [data.act_cartpos; exedata.act_cartpos];
        
        outdata.ff_u  = [data.ff_u; exedata.ff_u];
        outdata.fb_u  = [data.fb_u; exedata.fb_u];
        outdata.net_u  = [data.net_u; exedata.net_u];
        outdata.del_q = [data.del_q; exedata.del_q];
        outdata.accel = [data.accel; exedata.accel];
        outdata.zmp   = [data.zmp; exedata.zmp];
        outdata.isfalltrial = [data.isfalltrial; exedata.isfalltrial];
        outdata.trajErr     = [data.trajErr; exedata.trajErr];
        outdata.idyn_nopert = [data.idyn_nopert; exedata.idyn_nopert];
        outdata.idyn_pert   = [data.idyn_pert; exedata.idyn_pert];
        outdata.W   = [data.W, exedata.W];
    end
    
end

function [traj_area, init_traj_area] = signed_patharea(data)
    global EXE
    ofx = data.act_cartpos(2:end,1) - data.act_cartpos(1,1);
    %ofx(ofx<0)=0;  % to consider only the positive side area
    traj_area = 10000*sum(ofx.*abs(data.act_cartpos(2:end,2)-data.act_cartpos(1:end-1,2)));  % in cm^2
    %T0_inms = 500;  INI_endtime = T0_inms/1000; t0_ix =  round((INI_endtime / EXE.endtime)*EXE.timesteps);
    t0_ix = min(find(data.act_cartpos(:,2) >= data.act_cartpos(1,2)+0.025*1));
    init_traj_area = 10000*sum(ofx(1:t0_ix-1,1).*(data.act_cartpos(2:t0_ix,2)-data.act_cartpos(1:t0_ix-1,2)));  %in cm^2
end
    
function plot_data(data, par, figno, str)
    global EXE perfectIDYN SILENT
    
    if SILENT return; end
    figure(102); clf;
    winsize = 500;
    N = size(data.des_xtraj,1);
    if N>winsize
        sel = N-winsize:N;
    else
        sel = 1:N;
    end
    
    subplot(2,1,1);
    des_q = data.des_xtraj(sel,[1,3]);
    act_q = data.act_xtraj(sel,[1,3]);
    tsel = sel*(EXE.endtime/EXE.timesteps);

    plot(tsel,des_q,'--', 'Linewidth',2); hold on;
    plot(tsel,act_q,'-',  'Linewidth',2); 
    legend('des_{q1}','des_{q2}','act_{q1}','act_{q2}');
    xlabel('time(s)');
    title(sprintf('tracking error (%s)',str));
    subplot(2,1,2);
    
    %isfall = data.isfalltrial(sel);
    %bar(tsel(isfall>0), isfall(isfall>0)); hold on;
    plot(tsel,data.del_q(sel,:),'-', 'Linewidth',2); hold on;
    %plot(tsel, data.isfalltrial(sel),'g-','Linewidth',5); 
    plot(tsel, data.zmp(sel),'k');
    legend('del_{q1}','del_{q2}');

    figure(figno); clf;
    subplot(2,1,1);
    if perfectIDYN==1
        title(sprintf('Using perfect inverse dynmics!'));
    else
        pred_torques = invdyn_model_pred(data.act_xtraj(sel,:), data.accel(sel,:), par)'; 
        plot(tsel, data.net_u(sel,:),'Linewidth',2); hold on;
        plot(tsel,pred_torques,'--','Linewidth',2);
        legend('torq_1','torq_2','pred-torq1','pred-torq2');
        xlabel('time(s)');
        subplot(2,1,2);
        err = data.net_u(sel,:) - pred_torques;
        plot(tsel,abs(err));
        merr = mean(abs(err));
        xlabel('time(s)');
        ylim([-3*max(merr), 3*max(merr)]);
        title(sprintf('InvDyn: recent abs err:%3.4f, %3.4f (%s)',merr(1), merr(2),str));
    end
end

function add_to_rbuf(x,accel, torque)
    global EXE RBUF RBUFc RBUFfull

    RBUFc = RBUFc+1;
    if RBUFc > EXE.RBUFsize
        RBUFc=1;
        RBUFfull=1;
    end

    RBUF.input(RBUFc,:) = expandX(x, accel);
    RBUF.output(RBUFc,:) = torque;
    RBUF.accel(RBUFc,:) = accel;
end
        
function [xin,yout,accel] = sample_rbuf
    global EXE RBUF RBUFc RBUFfull
    if RBUFfull
        sampsz = min(EXE.RBUFsize, EXE.RBUFbatchsize);
        lastD = EXE.RBUFsize;
    else
        sampsz = min(RBUFc, EXE.RBUFbatchsize);
        lastD = RBUFc;
    end
 
    ix = randperm(lastD, sampsz);
    xin = RBUF.input(ix,:);
    yout = RBUF.output(ix,:);
    accel = RBUF.accel(ix,:);
end

% Using the experience buffer learn the idynamicsof the system
function [err_out, err_max] = learn_batch(rate, ler_type, par)
    global EXE RBUF RBUFc RBUFfull 
    global W dynnet 

        
    [Xexpanded,Yactual,accel] = sample_rbuf;              % take a sample set from the big buffer
    Xorg = Xexpanded(:,1:4);                              % First 4 is the normal state 
    if EXE.offset_dyn==1
        unpert_torques = invdyn_true(Xorg, accel, [par(1),0])';
        Y = (Yactual - unpert_torques)*EXE.offset_SC;         % when learning unperturbed_known_dyn + perturbed_EXE.offset_dyn
    else
        Y = Yactual;
    end
    %err_in = mean(abs(Y - X*W));
    if ler_type==1                          % my simple grad descent (probably buggy)
        YY = Xexpanded*W;
        dEdW1 = (YY(:,1)-Y(:,1))'*Xexpanded;
        dEdW2 = (YY(:,2)-Y(:,2))'*Xexpanded;
        delW = [dEdW1; dEdW2]';
        W = W - 1*rate*delW;
        pred_torques = Xexpanded*W;
    end
    if ler_type==0                                 % averaged sample based pseudo inv solution
        beta = 0.5;
        W_samp = pinv(Xexpanded)*Y;
        delW = W - (W_samp*beta + (1-beta)*W);     % make an weight average of the current estimate to avoid large jumps. (beta can be tuned.)
        W = W - delW;   
        %W = W_samp*beta + (1-beta)*W;
        pred_torques = Xexpanded*W;
    end
    if ler_type==2                                 % use matlabs nn toolbox
        [dynnet,tr] = train(dynnet,Xexpanded', Y');
        pred_torques = dynnet(Xexpanded')';
    end
    
    errs = abs(Y - pred_torques)/EXE.offset_SC;
    err_out = mean(errs);
    err_max = max(errs);
    %fprintf('del Error: %1.7f\n',err_out - err_in);
end

function [W, err] = est_idyn(data, type, par)
    global dynnet
    global EXE
    X = expandX(data.act_xtraj, data.accel);
    if EXE.offset_dyn==1
        unpert_torques = invdyn_true(data.act_xtraj, data.accel, [par(1),0])';
        Y = (data.net_u - unpert_torques)*EXE.offset_SC;
    else
        Y = data.net_u;
    end
    
    if (type == 1) ||(type == 0)
        W = pinv(X)*Y;
        pred_torques = X*W;
        err = abs(Y - pred_torques);
    else
        [dynnet,tr] = train(dynnet,X',Y');
        pred_torques = dynnet(X');
        err = abs(Y' - pred_torques)';
        W = dynnet;
    end
    
end


function trajSP = make_trajSP(trajpar, endTime)
    global body SPpar

    % B-spline properties
    k = 6;      % spline order (4 for cubic, 6 for quintic B-spline)
    kb = 3;     % spline order at the boundaries (1-pos, 2-pos,vel, 3-pos,vel,acc)
    NumCPm = SPpar.NumCPm; % number of control points to be optimized
    % boundary conditions

    if size(trajpar,2)==1   % if it is given flat make 2x..
        trajpar = reshape(trajpar,length(trajpar)/2,2)';
    end
    
    % control points
    CPb = body.SQUAT_STATE;                % begin value (start state)
    CPe = body.EXTENDED_STATE;             % end value (end state)

 
    cart_CPm = trajpar; 

  
    CPm = inverse_kinM(cart_CPm);   % Convert to joint angles
    CP = [repmat(CPb,1,kb),CPm,repmat(CPe,1,kb)];     % all control points including start and end
    NumCP = size(CP,2);                               % total number of control points
    % knots
    knots = augknt(linspace(0,endTime,NumCP-k+2),k);      % knot generation
    % B-spline generation
    [q1s, q2s] = make_splines(knots, CP); 
    % We have N cartesian control points. I make them into states of 6-dim
    % vector
    N = NumCPm*2;
    mu = [reshape(cart_CPm',N,1)];                    % state are the knot configuration (CPm)

    CP = state2CPm(mu);
    CP = [repmat(CPb,1,kb),CPm,repmat(CPe,1,kb)]; 
    q1s.coefs = CP(1,:);
    q2s.coefs = CP(2,:);

    dq1s = fnder(q1s);
    dq2s = fnder(q2s);

    ddq1s = fnder(dq1s);
    ddq2s = fnder(dq2s);

    trajSP.q1s = q1s;
    trajSP.q2s = q2s;
    trajSP.dq1s = dq1s;
    trajSP.dq2s = dq2s;
    trajSP.ddq1s = ddq1s;
    trajSP.ddq2s = ddq2s;
    trajSP.CP = CP;
    trajSP.CPb = CPb;                % begin value
    trajSP.CPe = CPe;                  % end value
end

function [compcost, safetycost, fallcost, fallstat, data] =execute_movement(x0, trajSP, fb_MUL, timeSteps, endTime, par, showlip)
    global first_roll  body SILENT EXE cart_grad
    if ~exist('showlip','var')
        showlip = 0;
    end
    
    
    t = linspace(0,endTime,timeSteps)';
    dt = t(2)-t(1);
    x  = x0;

    [Jsafety, Jeffort,fallcost, p, q, data] = run_control(trajSP,  fb_MUL, x0, t, par);
    zmp = data.zmp;
    if body.cost_type == -1
        compcost = Jeffort;
    elseif body.cost_type == 1
        compcost = Jsafety;
    elseif body.cost_type == 0         % composite/mixed cost
       compcost = Jeffort*body.effortW + (1-body.effortW)*Jsafety;
    end
    safetycost = Jsafety;
    data.compcost = compcost;
    data.safetycost = Jsafety;
    
    failed = (abs(fallcost)>body.fallcost_threshold);  
    if (~failed) 
        fallstat = EXE.fallstat_SUCC;
    else
        if (data.fallcost_fr > data.fallcost_ba)
            fallstat = EXE.fallstat_FFALL;
        else
            fallstat = EXE.fallstat_BFALL;
        end
    end
    if fallcost<0  % external fall xnoiseon
        fallstat = -fallstat;
    end
    first_roll = 0; 
    %figno = max([showlip, (1-sign(showlip))*999]);  % open figno 999 if showlip==0
    %figure(figno);
    if ~SILENT
        if (showlip~=0)
            figure(showlip);
        else
            figure(999); cla; 
            axis([-0.4,0.4, 0, 1.3]);
        end
    end
    
   % set(gcf,'Position',[320,0,300,350])
    %grid on; grid minor;
    if showlip & ~SILENT >0
        clf;
        %plot(p(:,1), p(:,2))
        [maxfatx , maxfatix]  = draw_controlled_LIP(q(:,1),q(:,2),t,body,0.000);
        plot (zmp, p(:,2),'g-');
        plot(data.des_cartpos(:,1), data.des_cartpos(:,2),'k-','Linewidth', 1.5);
        plot(data.act_cartpos(:,1), data.act_cartpos(:,2),'m.','Linewidth', 2);
        
    end
    if ~SILENT,     
        show_CPgrad(trajSP.CP, cart_grad);
        if (failed)
            xlabel(sprintf('traj err:(%2.3f,  %2.3f)',data.trajErr(1),data.trajErr(2)));        
            title(sprintf('Jeff:%1.3f Jsaf:%1.3f COST:%1.3f [FAIL:%s!]', Jeffort, Jsafety, compcost, EXE.fallstat_txt{4+fallstat}),'FontSize', 9);
        else
            xlabel(sprintf('traj err:(%2.3f,  %2.3f)',data.trajErr(1),data.trajErr(2)));
            title(sprintf('Jeff:%1.3f Jsaf:%1.3f COST:%1.3f[SUCCESS]',Jeffort, Jsafety, compcost),'FontSize', 9);
        end
    end
end

function [Jsafety, Jeffort, fallcost, p, q, data] = run_control(trajSP,  fb_MUL, x0,t, par)
    global W trialc feedback_on fb_gainP_base fb_gainD_base fb_delay_ms body EXE

    N = length(t);
    dt = t(2)-t(1);
    zmp = zeros(N,1);
    p  = zeros(N,2);
    q   = zeros(N,2);
    dq   = zeros(N,2);

    Jsafety = 0;
    Jeffort = 0;
    fallcost_fr = 0;
    fallcost_ba = 0;
    fallcost_he = 0;
    trajErr = 0;
    x = x0;
    fb_u = zeros(2,1);
    data.des_xtraj = zeros(N,4);
    data.act_xtraj = zeros(N,4);
    data.des_cartpos = zeros(N,2);
    data.act_cartpos = zeros(N,2);
    data.ff_u =zeros(N,2);
    data.fb_u =zeros(N,2);
    data.net_u = zeros(N,2);
    data.del_q =zeros(N,2);
    data.accel = zeros(N,2);

    data.trajErr  = 0;
    fb_delay_sec=  fb_delay_ms/1000.0; %sec

    delayCycle = floor(fb_delay_sec/dt);
    makexnoise = 0;
    if body.xnoiseon && rand>0.2
        makexnoise = 1;
    end
     
    FALLING = false;
        
    for i = 1:N
        [des_x, des_acc] = spline2vec(trajSP.q1s,trajSP.q2s,trajSP.dq1s,trajSP.dq2s, trajSP.ddq1s, trajSP.ddq2s, t(i));
        
        data.des_xtraj(i,:) = des_x;
        data.act_xtraj(i,:) = x;
        data.des_cartpos(i,:) = forward_kin([des_x(:,1),des_x(:,3)]);
        data.act_cartpos(i,:) = forward_kin([x(:,1),x(:,3)]);
        
        delayed_i = max([i-delayCycle, 1]);
        delayed_x  = data.act_xtraj(delayed_i,:);
        
        delayed_x(1) = delayed_x(1) + delayed_x(2)*delayCycle*dt/2;   % linear prediction to compensate delay 
        delayed_x(2) = delayed_x(2) + delayed_x(3)*delayCycle*dt/2;
        
        del_q = des_x([1,3])-delayed_x([1,3]);                        % error in joint tracking
        del_dq = 0*des_x([2,4])-delayed_x([2,4]);                     % error in joint velocities (taking des vel=0 for stability?)
        
        data.trajErr = data.trajErr+abs(del_q);
        
        fb_u(1) = fb_MUL*fb_gainP_base(1)*del_q(1) + fb_MUL*fb_gainD_base(1)*del_dq(1);
        fb_u(2) = fb_MUL*fb_gainP_base(2)*del_q(2) + fb_MUL*fb_gainD_base(2)*del_dq(2);
        
        
        ff_in_x = des_x;

        
        torques_off =  invdyn_true(ff_in_x, des_acc, [par(1) 0*par(2)]); % use des_x or x ?
        torques_ok =  invdyn_true(ff_in_x, des_acc, [par(1) par(2)]); % use des_x or x ?

        ff_u = invdyn_model_pred(ff_in_x, des_acc, par);

        u = ff_u + fb_u;
        if (FALLING)
            u(2) = 0;
        end

        data.idyn_nopert(i,:) = torques_off;
        data.idyn_pert(i,:) = torques_ok;
        data.ff_u(i,:) = ff_u';
        data.fb_u(i,:) = fb_u';
        data.net_u(i,:) = u';
        data.del_q(i,:) = del_q;
 
        if makexnoise && i<N/3
            xnoise = 120;
        else
            xnoise =  0;
        end
    
        mot_noise = [0; xnoise]; 
        dx = forward_dyn_fn(x,u+mot_noise,par);  % simulation
        ddq1 = dx(2);
        ddq2 = dx(4);
        data.accel(i,:) = [ddq1, ddq2];
        
        add_to_rbuf(x,[ddq1, ddq2], u');   % store the ([x,accell], u) mapping for learning
        
        pos = forward_kin([x(1),x(3)]);

        x = x + dt*dx;                     % tick simulation (Euler integ.)
        x(1) = max(0.45,min(1.2, x(1)));   % linear motion limited bewtween 0.45 - 1.2 m
        ZMP = zmp_equation(x(1),x(2),ddq1,x(3),x(4),ddq2,par);
        F_net = netforce_equation(x(1),x(2),ddq1,x(3),x(4),ddq2,par);

        zmp(i) = ZMP;
        q(i,:) = [x(1), x(3)];
        dq(i,:) = [x(2), x(4)];
        p(i,:) = pos;

        count_err = pos(2)<1.1 && pos(2)>0.55;                          % note: avoiding error accumulation towards the end heig
        
        if ZMP>body.toe_dist && count_err
            if (~FALLING)      
                if (EXE.noCONTafterFALL==1) 
                    fallcost_fr = fallcost_fr  + (ZMP - body.toe_dist)*0 + (N-i)+0*(x(2))^2; 
                else   
                    fallcost_fr = fallcost_fr + (ZMP - body.toe_dist);
                    %fallcost_fr = fallcost_fr + (ZMP - body.toe_dist)*(N-i)/N; %*(N-i)/N;%*(N-i)*dt*1;
                    %fallcost_fr = fallcost_fr + abs(ZMP - 0.9*body.toe_dist)*(N-i)/N;  
                end
            end
            if EXE.noCONTafterFALL==1, FALLING = true; end
        end
        if ZMP<-body.heel_dist && count_err
            if (~FALLING)
                if (EXE.noCONTafterFALL==1) 
                    fallcost_ba = fallcost_ba  + (ZMP + body.heel_dist)*0 +(N-i)+ 0*abs(x(2))^2; 
                else 
                    fallcost_ba =  fallcost_ba - (ZMP + body.heel_dist) ; 
                end
                %fallcost_ba =  fallcost_ba - (ZMP + body.heel_dist)*(N-i)/N; %*(N-i)/N; %*abs(x(4));*(N-i)*dt*1;
                %fallcost_ba =  fallcost_ba + abs(ZMP - 0.9*body.toe_dist)*(N-i)/N; %*(N-i)/N; %*abs(x(4));*(N-i)*dt*1;
            end
            if EXE.noCONTafterFALL==1, FALLING = true; end
            %fprintf('Heel zmp bad: %1.4f height:%1.3f\n', ZMP,pos(2));
        end
         
        Jsafety = Jsafety + count_err*abs(ZMP - 0.5*body.toe_dist)  ;
        %Jsafety = max(Jsafety, 0.5*N*err_range*abs(ZMP - 0.95*body.toe_dist));
        %Jsafety = Jsafety + err_range*abs(pos(1) - 0.95*body.toe_dist);
        %Jsafety = max(Jsafety, 0.5*N*err_range*abs(pos(1) - 0.95*body.toe_dist));
        Jeffort = Jeffort + count_err*(abs(u(1)) + 6.3766*abs(u(2)))/1000; % 6.3766 normalizes the actuation based on maximal absolute torques
        %Jenergy = Jenergy + err_range*(u(1)^2 + 6.3766*u(2)^2)^0.5/1000;
        
    end
    
    data.zmp = zmp;
    data.tm  = t;
    % If you want to use max zmp as fallcost
    %frontfall = zmp(zmp>body.toe_dist & p(:,2)<1.15) - body.toe_dist;
    %backfall  = -zmp(zmp <-heel_dist & p(:,2)<1.15) + heel_dist;
    %fallcost  = 20*max([0;frontfall; backfall]);
    
    % if endpointis not close to 0,1.2 consider it as a failure
    if norm(pos-body.EXTENDED_CART) > 0.5
        fallcost_he = fallcost_he + norm(pos-body.EXTENDED_CART);
    end
    
    data.trajErr = 1000*data.trajErr/N;
    Jsafety = 10*Jsafety^2/N;
    Jeffort = 10*Jeffort/N;
    fallcost_fr =  50*fallcost_fr^2/N;
    fallcost_ba =  50*fallcost_ba^2/N;
    fallcost = (fallcost_fr + fallcost_ba + 0*fallcost_he)*(1-2*makexnoise) ;
    data.Jsafety = Jsafety;
    data.Jeffort = Jeffort;
    data.fallcost = fallcost;
    data.fallcost_fr = fallcost_fr;
    data.fallcost_ba = fallcost_ba;
    data.fallcost_he = fallcost_he;
end

function [q1s, q2s] = make_splines(knots, CP)
    q1s = spmak(knots,CP(1,:));
    q2s = spmak(knots,CP(2,:)); 
end

function  CP  = state2CPm(s)
    Dxy = reshape(s,length(s)/2,2)';
    CP = inverse_kinM(Dxy);
end


function [q1, q2] = inverse_kin(x,y)
    q2 = atan2(y,x);
    q1 = (x.^2 + y.^2).^0.5;
end
function Dlq = inverse_kinM(Dxy)
    x = Dxy(1,:)';
    y = Dxy(2,:)';
    [q1,q2] = inverse_kin(x,y);
     Dlq = [q1'; q2'];
end

function show_CP(CP)
    xx = cos(CP(2,3:7)).*CP(1,3:7); 
    yy = sin(CP(2,3:7)).*CP(1,3:7);
    plot(xx,yy,'ro'); hold on;
    plot(xx, yy,'b-'); 

    drawnow;
end

function show_CPgrad(CP, grad)
    xx = cos(CP(2,3:7)).*CP(1,3:7); 
    yy = sin(CP(2,3:7)).*CP(1,3:7);
    plot(xx,yy,'ro'); hold on;
    plot(xx, yy,'b-'); 
    dim = size(grad,1)/2;
    GRIX = 1;  %1:CompCOST 2:FailCOST 3:SafetyCOST
    grSC = 0.01;
    for k = 1:dim
        %plot([xx(k+1),xx(k+1)-grad(k, 1)*grSC],[xx(k+1),yy(k+1)-grad(k+dim,1)*grSC],'r--');  
        quiver(xx(k+1),yy(k+1),-grad(k, 1)*grSC,-grad(k+dim,1)*grSC,0,'r');
        quiver(xx(k+1),yy(k+1),-grad(k, 2)*grSC,-grad(k+dim,2)*grSC,0,'g');
        quiver(xx(k+1),yy(k+1),-grad(k, 3)*grSC,-grad(k+dim,3)*grSC,0,'b');
        %plot([xx(k+1),xx(k+1)-grad(k, 2)*grSC],[yy(k+1),yy(k+1)-grad(k+dim,2)*grSC],'g--');  
        %plot([xx(k+1),xx(k+1)-grad(k, 3)*grSC],[yy(k+1),yy(k+1)-grad(k+dim,3)*grSC],'k--');  
    end
    drawnow;
end

function show_state(s)
    c = state2CPm(s);
    xx = cos(c(2,:)).*c(1,:); 
    yy = sin(c(2,:)).*c(1,:);
    plot(xx,yy,'ro'); hold on;
    plot(xx, yy,'b-'); 
end

function p = forward_kin(q)
p = [ q(:,1).*cos(q(:,2)), q(:,1).*sin(q(:,2))];
end
function dx = forward_dyn_fn(x,T,par)
m = par(1);
K = par(2);
g = 9.81;
 
T1 = T(1);
T2 = T(2);
q1  = x(1);
dq1 = x(2);
q2  = x(3);
dq2 = x(4);


ddq1 = -(g*m*sin(q2) - dq2^2*m*q1 - T1 + K*dq2*q1*cos(q2)^2 + K*dq1*cos(q2)*sin(q2))/m;
ddq2 = (T2 - 2*dq1*dq2*m*q1 - g*m*q1*cos(q2) + K*dq1*q1*sin(q2)^2 + K*dq2*q1^2*cos(q2)*sin(q2))/(m*q1^2);
       
dx = [x(2), ddq1, x(4), ddq2];

end
function dx = forward_dyn_fn_nograv(x,T,par)
m = par(1);
K = par(2);
g = 0;
 
T1 = T(1);
T2 = T(2);
q1  = x(1);
dq1 = x(2);
q2  = x(3);
dq2 = x(4);


ddq1 = -(g*m*sin(q2) - dq2^2*m*q1 - T1 + K*dq2*q1*cos(q2)^2 + K*dq1*cos(q2)*sin(q2))/m;
ddq2 = (T2 - 2*dq1*dq2*m*q1 - g*m*q1*cos(q2) + K*dq1*q1*sin(q2)^2 + K*dq2*q1^2*cos(q2)*sin(q2))/(m*q1^2);
       
dx = [x(2), ddq1, x(4), ddq2];

end
function [T]=inv_dyn_fn(q1,dq1,ddq1,q2,dq2,ddq2,par)

g = 9.81;

m = par(1);
K = par(2);
 
T = [ m*ddq1 - m*(q1*dq2^2 - g*sin(q2)) + K*cos(q2)*(dq1*sin(q2) + dq2*q1*cos(q2)); m*q1^2*ddq2 + m*q1*(2*dq1*dq2 + g*cos(q2)) - K*q1*sin(q2)*(dq1*sin(q2) + dq2*q1*cos(q2))];
 
end


function [des_x, des_acc] = spline2vec(q1s, q2s, dq1s, dq2s, ddq1s, ddq2s, time)
    q1 = fnval(q1s,time);    
    q2 = fnval(q2s,time);

    dq1 = fnval(dq1s,time);
    dq2 = fnval(dq2s,time);

    ddq1 = fnval(ddq1s,time);
    ddq2 = fnval(ddq2s,time);

    des_x = [q1,dq1, q2, dq2];
    des_acc = [ddq1, ddq2];
end
        
function torques_u = invdyn_true(x, des_acc, par)
    torques_u = my_inv_dyn(x(:,1),x(:,2),des_acc(:,1),x(:,3),x(:,4),des_acc(:,2),par);
end

% make sure that first elements are always the original X!
function hX = expandX(Xin, accel)
    X = [Xin, accel];
    hX = [X, sin(X(:,2)),  X.^2, X.^3];
    %hX = X;
    %hX = [X, sin(X(:,2)), accel];
end

% function hX = expandX(X, accel)
function pred_torques = invdyn_model_pred(x, des_acc,par)
    global W dynnet EXE first_roll
    global perfectIDYN
    if perfectIDYN==1    %if using perfect idyn 
        pred_torques = invdyn_true(x, des_acc, par);
        return
    end
        
    pred_torques = 0;
    if EXE.offset_dyn == 1
        unpert_torques = invdyn_true(x, des_acc, [par(1),0]);
        divider = EXE.offset_SC;    
    else
        unpert_torques = 0;
        divider = 1; 
    end
    
    if first_roll==1    %if nothing learned yet return the unperturbed dynamics
        pred_torques = invdyn_true(x, des_acc, [par(1),0]);
    else
        X = expandX(x, des_acc);
        if (EXE.learning_type==0) || (EXE.learning_type==1)
            pred_torques = (X*W)'/divider+unpert_torques;
        end
        if (EXE.learning_type==2)
            pred_torques = dynnet(X')/divider+unpert_torques;
        end
        
    end
end

% function err = invdyn_model_learn(x, u, acc, time, par)
%     % learn data point (x,acc) -> u
%     return;
% end

function [T]= my_inv_dyn(q1,dq1,ddq1,q2,dq2,ddq2,par)

g = 9.81;

m = par(1);
K = par(2);
 
T = [ m*ddq1 - m*(q1.*dq2.^2 - g*sin(q2)) + K*cos(q2).*(dq1.*sin(q2) + dq2.*q1.*cos(q2)),...
      m*q1.^2.*ddq2 + m*q1.*(2*dq1.*dq2 + g*cos(q2)) - K*q1.*sin(q2).*(dq1.*sin(q2) + dq2.*q1.*cos(q2))]';
 
end

function [delcc, delfc , delsc] = numeric_cost_grads(a,fb_MUL, x0, timeSteps, endTime, par)
global first_roll  body SILENT EXE  perfectIDYN
    orgperfectIDYN = perfectIDYN;
    orgbodyxnoiseon = body.xnoiseon
    
    body.xnoiseon = 0;
    perfectIDYN = 1;
    
    myeps = 0.001;
    t = linspace(0,endTime,timeSteps)';
    dt = t(2)-t(1);
    
    [cost0, fallcost0,safetycost0, ~,~,~] = objfun(a, fb_MUL, x0, t, par);

    delcc = a*0; delfc = delcc; delsc = delcc;
    for k=1:length(a)
        a(k) = a(k) + myeps;
        [cost1, fallcost1,safetycost1, ~,~,~] = objfun(a, fb_MUL, x0, t, par);
        a(k) = a(k) - myeps;
        delcc(k) = (cost1-cost0)/myeps;
        delfc(k) = (fallcost1-fallcost0)/myeps;
        delsc(k) = (safetycost1-safetycost0)/myeps;
    end
    if (sum(isnan(delcc))+sum(isnan(delfc))+sum(isnan(delsc)))>0,
        fprintf('numeric deriv has NAN!\n');
    end
    perfectIDYN = orgperfectIDYN;
    body.xnoiseon = orgbodyxnoiseon;
end



function [compcost, fallcost, safetycost, zmp, xy, T_his] = objfun(a, fb_MUL, x0, t, par)
global EXE body
    trajSP = make_trajSP(a, EXE.endtime);
    dt = t(2)-t(1);
    [Jsafety, Jeffort,fallcost, p, q, data] = run_control(trajSP,  fb_MUL, x0, t, par);
    if body.cost_type == -1
        compcost = Jeffort;
    elseif body.cost_type == 1
        compcost = Jsafety;
    elseif body.cost_type == 0         % composite/mixed cost
       compcost = Jeffort*body.effortW + (1-body.effortW)*Jsafety;
    end 
    safetycost = Jsafety;
    zmp = data.zmp;
    xy = data.act_cartpos;
    T_his = data.net_u;
    
end

function [obj] = fmin_objfun(a)
global EXE body myfminpars
    
    myfminpars.iter = myfminpars.iter+1;
    if (mod(myfminpars.iter,10)==0),
        fprintf('fminsearch %d > \n',myfminpars.iter);
    end
    trajSP = make_trajSP(a, EXE.endtime);
    %dt = t(2)-t(1);
    [Jsafety, Jeffort,fallcost, p, q, data] = run_control(trajSP,  myfminpars.fb_MUL , myfminpars.x0, myfminpars.t, myfminpars.par);
    if body.cost_type == -1
        compcost = Jeffort;
    elseif body.cost_type == 1
        compcost = Jsafety;
    elseif body.cost_type == 0         % composite/mixed cost
       compcost = Jeffort*body.effortW + (1-body.effortW)*Jsafety;
    end 
    %safetycost = Jsafety;
    %zmp = data.zmp;
    %xy = data.act_cartpos;
    %T_his = data.net_u;
    
    obj = 100*fallcost + compcost;
end

function [des_cartpos, des_xtraj] = trajSP2cart(trajSP)
    global EXE
    t = linspace(0,EXE.endtime,EXE.timesteps)';
    des_xtraj=[];
    des_cartpos=[];
    for i=1:length(t)
        [des_x, des_acc] = spline2vec(trajSP.q1s,trajSP.q2s,trajSP.dq1s,trajSP.dq2s, trajSP.ddq1s, trajSP.ddq2s, t(i));
        des_xtraj(i,:) = des_x;
        des_cartpos(i,:) = forward_kin([des_x(:,1),des_x(:,3)]);
    end
end


function frames2avi(F, fname, rate)    
    writerObj = VideoWriter(fname);
    writerObj.FrameRate = rate;
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F{i} ;    
        frame
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end  