
function [x, ddq1, ddq2] = state_from_spline(CPm,CPb,CPe,kb,q1s,q2s,time)
% Given the trajecotry splines  and time, it samples them at given time and
% returns the trajectory point a flat state vector.
CP = [repmat(CPb,1,kb),CPm,repmat(CPe,1,kb)]; 
q1s.coefs = CP(1,:);
q2s.coefs = CP(2,:);

dq1s = fnder(q1s);
dq2s = fnder(q2s);

ddq1s = fnder(dq1s);
ddq2s = fnder(dq2s);


  
    q1 = fnval(q1s,time);    
    q2 = fnval(q2s,time);
    
    dq1 = fnval(dq1s,time);
    dq2 = fnval(dq2s,time);
    
    ddq1 = fnval(ddq1s,time);
    ddq2 = fnval(ddq2s,time);
    
x = [q1 dq1 q2 dq2 ];

   