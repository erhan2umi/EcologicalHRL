function F_net=netforce_equation(q1,dq1,ddq1,q2,dq2,ddq2,par)

m = par(1);
K = par(2);

% COM postion calculations

R = [ q1*cos(q2), q1*sin(q2)];
dR = [ dq1*cos(q2) - dq2*q1*sin(q2), dq1*sin(q2) + dq2*q1*cos(q2)];
ddR = [ - q1*cos(q2)*dq2^2 - 2*dq1*sin(q2)*dq2 + ddq1*cos(q2) - ddq2*q1*sin(q2), - q1*sin(q2)*dq2^2 + 2*dq1*cos(q2)*dq2 + ddq1*sin(q2) + ddq2*q1*cos(q2)];

% Net Force calculation
F_net = -K*dR(2) + m*9.81*sin(2*q2)/2;


end

