function metSolver2
    format long
    %all the constants
        bet = 0.15;
        del = 0;
        mu = 0.01;
        thet0 = 0.39;
    timestep = 0.005;
    time = 0:timestep:2000;
    M = @(t,V)[1-bet*(cos(V(2))^2), bet*cos(V(1))*cos(V(2));
               bet*cos(V(1))*cos(V(2)), 1-bet*(cos(V(1))^2)] * [(1+del)*sin(V(1)) + mu*((V(1)/thet0)^2-1)*V(3) + bet*cos(V(1))*((V(3)^2)*sin(V(1)) + (V(4)^2)*sin(V(2)));
                                                                (1-del)*sin(V(2)) + mu*((V(2)/thet0)^2-1)*V(4) + bet*cos(V(2))*((V(3)^2)*sin(V(1)) + (V(4)^2)*sin(V(2)))];
    F = @(t,V)[V(3);
               V(4);
               M(t,V)/(bet*(cos(V(1))^2 + cos(V(2))^2) - 1)];
           
    IC = [0.7;-0.69;0;0];
    X = rk42d(F,time,IC);
    
    plot(time/10,X(1,:));
    hold on
    plot(time/10,X(2,:));
    hold off
    title('Angles of Metronomes');
    xlabel('Time (s)');
    ylabel('Angle (radians)');
    legend({'\theta_1','\theta_2'},'location','southwest');
    
    figure(2)
    plot(time/10,X(2,:)-X(1,:),'g');
    hold on
    plot(time/10,X(2,:)+X(1,:),'m');
    hold off
    xlabel('Time (s)')
    ylabel('Angle (radians)')
    title('Sum and difference of angles');
    legend({'\theta_2 - \theta_1','\theta_1 + \theta_2'},'location','southwest');
end