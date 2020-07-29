function metMalkin(timestep,bet,del,mu,thet0)
    format long
    F1 = @(t,V)[V(2);
        -((1-bet)/(1-2*bet))*((1+del)*V(1) + mu*((V(1)/thet0)^2-1)*V(2)) - ((bet)/(1-2*bet))*V(1)*(V(2)^2)];
    [IC,count] = fullMet1Limit(F1,timestep);
    time = 0:timestep:(count*timestep);
    
    %This is the Exact Jacobian
    F2 = @(t,V1,V2)-[0, ((bet-1)/(1-2*bet))*( 1+del + (2*mu*V1(1)*V1(2))/(thet0^2) ) - (bet*(V1(2)^2))/(1-2*bet);
                     1, ((bet-1)*mu*((V1(1)/thet0)^2-1))/(1-2*bet) - (2*bet*V1(1)*V1(2))/(1-2*bet)]*V2;
                 
    %This is the Jacobian found numerically
    %step = 0.000001;  %stepsize for finding jacobian
    %F2 = @(t,V1,V2)-[(F1(t,V1+[step;0])-F1(t,V1-[step;0]))/(step*2),(F1(t,V1+[0;step])-F1(t,V1-[0;step]))/(2*step)]'*V2;
    
    F = @(t,V)[F1(t,V(1:2));F2(t,V(1:2),V(3:4))];
    
    Qic = [0;1];
    norm = (Qic')*F1(0,IC);
    Qic = Qic/norm;
    ICV = [IC;Qic];
    XQ = rk42d(F,time,ICV);
    
    %Un comment the rest to plot graphs
    plot(XQ(1,:),XQ(2,:));
%     hold on
%     plot(time/10,XQ(2,:));
    hold off
    title('Limit Cycle');
    xlabel('\theta')
    ylabel('\omega')
    %legend({'\theta','\omega'},'location','southwest')
    
    figure(2)
    plot(XQ(3,:),XQ(4,:));
%     plot(time/10,XQ(3,:))
%     hold on
%     plot(time/10,XQ(4,:));
%     hold off
    title('Q')
    xlabel('Q_\theta')
    ylabel('Q_\omega')
    %legend({'Q_\theta','Q_\omega'},'location','northeast')
end