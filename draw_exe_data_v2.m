function [maxfatx , maxfatix]  = draw_exe_data_v2(exe_data, bodypar, t_pause)
   q1 = exe_data.act_xtraj(:,1);
   q2 = exe_data.act_xtraj(:,3);
   tm = exe_data.tm;
   [maxfatx , maxfatix] = draw_controlled_LIP(q1,q2,tm, bodypar, t_pause);

   plot (exe_data.zmp, exe_data.act_cartpos(:,2),'g-');
   plot(exe_data.des_cartpos(:,1), exe_data.des_cartpos(:,2),'k-','Linewidth', 1.5);
   plot(exe_data.act_cartpos(:,1), exe_data.act_cartpos(:,2),'m.','Linewidth', 2);
      
      if (exe_data.fallcost>bodypar.fallcost_threshold)
        xlabel(sprintf('traj err:(%2.3f,  %2.3f)',exe_data.trajErr(1),exe_data.trajErr(2)));        
        title(sprintf('Jeff:%1.3f Jsaf:%1.3f COST:%1.3f [FAIL!]', exe_data.Jeffort, exe_data.Jsafety, exe_data.compcost),'FontSize', 9);
    else
        xlabel(sprintf('traj err:(%2.3f,  %2.3f)',exe_data.trajErr(1),exe_data.trajErr(2)));
        title(sprintf('Jeff:%1.3f Jsaf:%1.3f COST:%1.3f[SUCCESS]',exe_data.Jeffort, exe_data.Jsafety, exe_data.compcost),'FontSize', 9);
    end
end

   