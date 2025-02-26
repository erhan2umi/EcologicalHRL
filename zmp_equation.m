function [ZMP]=zmp_equation(q1,dq1,ddq1,q2,dq2,ddq2,par)

m = par(1);
K = par(2);
g = 9.81;

ZMP = (m*q1*sin(q2)*(q1*cos(q2)*dq2^2 + 2*dq1*sin(q2)*dq2 - ddq1*cos(q2) + ddq2*q1*sin(q2)) - K*q1*sin(q2)*(dq1*sin(q2) + dq2*q1*cos(q2)) + m*q1*cos(q2)*(- q1*sin(q2)*dq2^2 + 2*dq1*cos(q2)*dq2 + g + ddq1*sin(q2) + ddq2*q1*cos(q2)))/(m*(- q1*sin(q2)*dq2^2 + 2*dq1*cos(q2)*dq2 + g + ddq1*sin(q2) + ddq2*q1*cos(q2)));

end

