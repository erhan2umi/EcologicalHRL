function [maxfatx , maxfatix]  = draw_controlled_LIP(q1,q2,t,bodypars,T_pause)


N=length(t);
p = forward_kin([q1,q2]);

if T_pause>0
    for i=1:N
        cla;
        plot(p(1:i,1),p(1:i,2),'m.'); hold on;
        drawBody_v2([q1(i),q2(i)], bodypars);  
        text(p(i,1)+0.1,p(i,2),sprintf('%1.3fm',t(i)));   
        drawnow; 
    end
else
    i=N;
    cla;
    plot(p(1:i,1),p(1:i,2),'m.'); hold on;
    drawBody_v2([q1(i),q2(i)], bodypars);  
    text(p(i,1)+0.1,p(i,2),sprintf('%1.3fm',t(i)));   
    plot([p(1),p(1)],[0,1.2],'--','Color',[1 0.5 0]);
    drawnow; 
end

[maxfatx , maxfatix] = max(p(:,1));
plot([0, maxfatx],[p(maxfatix,2), p(maxfatix,2)],'--');
text(maxfatx+0.05,p(maxfatix,2),sprintf('%1.3f',maxfatx));

function p = forward_kin(q)
p = [ q(:,1).*cos(q(:,2)), q(:,1).*sin(q(:,2))];

function p = drawBody_v2(q, bodypars)
ankleoff = bodypars.heel_dist;
footlen  = bodypars.heel_dist+bodypars.toe_dist;

p = [ q(1)*cos(q(2)), q(1)*sin(q(2))];

%plot([0,p(1)], [0,p(2)],'b-','LineWidth',3);
plot([0,p(1)], [0,p(2)],'Color',[0.5,0.5,0.5],'LineWidth',3);
plot([-ankleoff,footlen-ankleoff],[0,0],'k-','Linewidth',4);
plot(p(1), p(2),'ro','LineWidth',4); 
plot(p(1)/2,0,'rx');

axis equal;
axis([-0.5,+0.5, 0, 1.3]);

function show_CP(CP)
    xx = cos(CP(2,:)).*CP(1,:);
    yy = sin(CP(2,:)).*CP(1,:);
    plot(xx,yy,'x'); hold on;
    plot(xx, yy,'y-','Linewidth',2); 

    drawnow;
