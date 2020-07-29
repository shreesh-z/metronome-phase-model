function V = rk42d(F,t,IC)
    siz = size(IC);
    tstep = t(2)-t(1);
    if(siz(1)==1)
        V = zeros(length(t),length(IC));
        V(1,:) = IC;
        for i = 1:1:length(t)-1
            k1 = F(t(i),V(i,:));
            k2 = F(t(i)+tstep/2,V(i,:)+k1*tstep/2);
            k3 = F(t(i)+tstep/2,V(i,:)+k2*tstep/2);
            k4 = F(t(i)+tstep,V(i,:)+k3*tstep);

            V(i+1,:) = V(i,:) + (tstep/6)*(k1+2*(k2+k3)+k4);
        end
    else 
        if(siz(2)==1)
            V = zeros(length(IC),length(t));
            V(:,1) = IC;
            for i = 1:1:length(t)-1
                k1 = F(t(i),V(:,i));
                k2 = F(t(i)+tstep/2,V(:,i)+k1*tstep/2);
                k3 = F(t(i)+tstep/2,V(:,i)+k2*tstep/2);
                k4 = F(t(i)+tstep,V(:,i)+k3*tstep);

                V(:,i+1) = V(:,i) + (tstep/6)*(k1+2*(k2+k3)+k4);
            end
        end
    end
end