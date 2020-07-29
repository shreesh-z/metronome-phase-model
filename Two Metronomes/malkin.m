function XQ = malkin(timestep)
    %For van der Pol oscillator
    F1 = @(t,V)[V(1)-V(1)^3-V(2); V(1)];
    F2 = @(t,V1,V2)[(1-3*(V1(1)^2))*V2(1)+V2(2); -V2(1)];
    
%     for Andronov Hopf oscillator
%     F1 = @(t,V)[(V(1)-V(2))-V(1)*(V(1)^2+V(2)^2);
%                 (V(1)+V(2))-V(2)*(V(1)^2+V(2)^2)];
%     F2 = @(t,V1,V2)[(1-V1(2)^2-3*V1(1)^2)*V2(1) + (1-2*V1(1)*V1(2))*V2(2);
%                    -(1+2*V1(1)*V1(2))*V2(1) + (1-V1(1)^2-3*V1(2)^2)*V2(2)];
    
    %finding a point on the limit cycle and 
    %how many timesteps are needed to evolve it
    [IC,count] = fullMet1Limit(F1,timestep);
    
    
    %step = 0.0000001;
    %F2 = @(t,V1,V2)[(F1(t,V1+[step;0])-F1(t,V1-[step;0]))/(2*step),(F1(t,V1+[0;step])-F1(t,V1-[0;step]))/(2*step)]'*V2;
    
    F = @(t,V)[F1(t,V(1:2));F2(t,V(1:2),V(3:4))];
    time = 0:timestep:(timestep*count);
    
    %for van der Pol oscillator
    Qic = [-0.371195702701162;0.264733342201252];
    norm = Qic'*F1(0,IC);
    Qic = Qic/norm;
    ICV = [IC;Qic];
    
    %for Andronov Hopf oscillator
    %ICV = [1;0;0;1];
    
    XQ = rk42d(F,time,ICV);
    XQ(3:4,:) = flip(XQ(3:4,:));
    
    plot(XQ(1,:),XQ(2,:))
    %hold on
    %plot(time,XQ(2,:));
    %hold off
    title('Limit Cycle');
    xlabel('X')
    ylabel('Y')
    %legend({'x','y'},'location','southeast')
    
    figure(2)
    plot(time,XQ(3,:));
    hold on
    plot(time,XQ(4,:));
    hold off
    title('Q')
    xlabel('Time')
    ylabel('Variables')
    legend({'Q_x','Q_y'},'location','southwest')
end