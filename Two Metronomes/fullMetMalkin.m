function XQ = fullMetMalkin(timestep,bet,del,mu,thet0)
	
	F1 = @(t,V)[V(2);
				((1+del)*sin(V(1)) + mu*((V(1)/thet0)^2-1)*V(2)
				- (bet/2)*sin(2*V(1))*(V(2)^2))/(bet*(cos(V(1))^2) - 1)];
	F2 = @(t,V)-[0, ((1+del)*((bet*(cos(V(1))^2)-1)*cos(V(1))
					+ bet*sin(V(1))*sin(2*V(1)))
					+ mu*V(2)*(2*V(1)*(bet*(cos(V(1))^2)-1)
					+ bet*(V(1)^2)*sin(2*V(1)))/(thet0^2)
					+ mu*V(2)*bet*sin(2*V(1))
					- (bet*(V(2)^2)/2)*(2*(bet*(cos(V(1))^2)-1)*cos(2*V(1))
					+ bet*(sin(2*V(1))^2)))/((bet*(cos(V(1))^2)-1)^2);
				1, (mu*((V(1)/thet0)^2 - 1)
					- bet*V(2)*sin(2*V(1)))/(bet*cos(V(1)^2)-1)]
				*V(3:4);
	
	F = @(t,V)[F1(t,V(1:2)); F2(t,V)];
	
	[IC,count] = findLimitCycle(F1,timestep);
	time = 0:timestep:(count*timestep);
	Qic = [1;1];
	norm = (Qic')*F1(0,IC);
	Qic = Qic/norm;
	ICV = [IC;Qic];
	XQ = rk42d(F,time,ICV);
	
	plot(time/10,XQ(1,:))
	hold on 
	plot(time/10,XQ(2,:));
	hold off 
	title('Limit Cycle');
	xlabel('Time')
	ylabel('Variables')
	legend({'\theta','\omega'},'location','southeast')
	
	figure(2)
	plot(time/10,XQ(3,:));
	hold on
	plot(time/10,XQ(4,:));
	hold off
	title('Q, \theta_0 = 0.59')
	xlabel('Time')
	ylabel('Variables')
	legend({'Q_{\theta}','Q_{\omega}'},'location','southwest')
	grid on
end