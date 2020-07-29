function metH(timestep)
    %all the constants
        bet = 0.011;
        del = 0;
        mu = 0.01;
        thet0 = 0.2;
    %XQ = metMalkin(timestep,bet,del,mu,thet0); %for small angle eqs
    XQ = fullMet1Malkin(timestep,bet,del,mu,thet0); %for full eqs
    sizXQ = size(XQ);
    count = sizXQ(2);
    timevect = 0:timestep:(count*timestep);
    func1 = zeros(count+1,1);
    H = zeros(count+1,1);
    H1 = zeros(count+1,1);
    
    %interaction terms for small angle approx
    %F = @(V)-((bet)/(1-2*bet))*(( (1+del)*V(1) + mu*((V(1)/thet0)^2-1)*V(2)) - V(1)*V(2)^2);
    
    %interaction terms for full set of eqs
    %D is the full set of equations for both the oscills
    %Subtracting D from the expression for a free running oscill
    %   gives us the interaction terms
    D = @(V1,V2)[1-bet*(cos(V2(1))^2), bet*cos(V1(1))*cos(V2(1))] * [(1+del)*sin(V1(1)) + mu*((V1(1)/thet0)^2-1)*V1(2) + bet*cos(V1(1))*((V1(2)^2)*sin(V1(1)) + (V2(2)^2)*sin(V2(1)));
                                                                     (1-del)*sin(V2(1)) + mu*((V2(1)/thet0)^2-1)*V2(2) + bet*cos(V2(1))*((V1(2)^2)*sin(V1(1)) + (V2(2)^2)*sin(V2(1)))];
    F = @(V1,V2)((1+del)*sin(V1(1)) + mu*((V1(1)/thet0)^2-1)*V1(2) - (bet/2)*sin(2*V1(1))*(V1(2)^2))/(1-bet*(cos(V1(1))^2)) - D(V1,V2)/(1-bet*(cos(V1(1))^2 + cos(V2(1))^2));

    for phasediff = 0:1:count
        for phase = 1:1:count
            offset = phase+phasediff;
            while(offset>count)
                offset = offset - count;
            end
			%for full equations
            func1(phase) = XQ(4,phase)*F(XQ(1:2,phase),XQ(1:2,offset));
            %for small angle eqs
			%func1(phase) = XQ(4,phase)*F(XQ(1:2,offset));
        end
        func1(count+1) = func1(1);
        sum1 = simpson(func1);
        H(phasediff+1) = sum1/(count);
    end
    H1 = flip(H);
    G = H1-H;
	
    figure(3)
    plot(timevect/10,H,'b');
    title('Phase Model')
    xlabel('Phase Difference')
    ylabel('H')
    legend('H_{12}(\phi_2 - \phi_1)','location','southwest')
    grid on
    
    figure(2)
    plot(timevect/10,G,'g');
    title('G')
    xlabel('Phase Difference')
    ylabel('G')
    legend('G(\chi)','location','southeast')
    grid on
end
function sum = simpson(func)
    count = length(func);
    sum = (func(1)+func(count));
    for time = 2:1:count-1
        if(mod(time,2)==0)
            sum = sum+4*func(time);
        else
            sum = sum+2*func(time);
        end
    end
    sum = sum/3;
end
function sum = trap(func)
    count = length(func);
    sum = (func(1)+func(count))/2;
    for time = 2:1:(count-1)
        sum = sum + func(time);
    end
end