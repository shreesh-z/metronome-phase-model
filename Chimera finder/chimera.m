%all physical params

N = 15;
freqBpm = 160;
m = 0.028;
M = 2.31;
l = 0.15;
L = 0.22;
k = 68;
g = 9.81;

l_bob = 0.073 - 0.00022*freqBpm;
r_cm = abs(0.178*l_bob - 0.0121);
I = 0.0000129 + 0.005*(l_bob^2);

%derived params

x0 = (m*r_cm)/M;
omega = sqrt((m*g*r_cm)/I);
kappa = (k/M)*((l/L)^2);
Omega = sqrt(g/L);

%all the constants:

mu_m = 0.011;         %nonlinearity of metros
mu_s = 0.00016;      %damping of swings
thet0 = 0.33;

%all the params that can be customized

%kappa = 20;        #coupling of swings
%freq = freqBpm/60;       #freq of metros in Hz
%omega = 2*pi*freq;
%x0 /= 100; 


%nondimensional parameters
beta = (x0*(omega^2))/g;
omeg_r2 = (Omega/omega)^2;
chi = kappa/(omega^2);

disp("x0        : "+x0)
disp("omega     : "+omega)
disp("kappa     : "+kappa)
disp("Beta      : "+beta)
disp("Omega^2_r : "+omeg_r2)
disp("Chi       : "+chi)

timestep = 0.01;
maxTime = 50;
oscill = round(maxTime/(2*pi)) + 2;

disp("Num of oscills : "+maxTime/(2*pi));

timeVec = 0:timestep:maxTime;
timeCount = length(timeVec);
disp("Count is: "+timeCount);

V = zeros(2*N+2,2,timeCount);
R = zeros(2*N);
F = zeros(2*N+2,2,4);
V_ = zeros(2*N+2,2);
cos_Phi_Psi = zeros(N,1);

V(1:N,1,1) = 2*thet0;   %for IP IC
%V(N+1:2*N,1:2,1) = 2*thet0*rand(N,2) - thet0;
V(N+1:2*N,1,1) = [0.14157382; 0.2614215; -0.29608024; -0.24948547; -0.07842186; -0.31545333; 0.06066683; 0.27495881; 0.32353985; -0.19375554; 0.07076737; 0.12151782; -0.28229673; 0.14350184; 0.03906534];
V(N+1:2*N,2,1) = [0.04105613; -0.22835415;  0.05830795; -0.20626723; -0.15800162; -0.21524905; 0.18515062; 0.05557701; -0.13233109; -0.31297807; -0.01390028; 0.28235838; 0.0213475; -0.20703707; -0.04162385];

for t = 1:(timeCount-1)
	V_ = V(:,:,t);
	for rk = 1:4
		
		%function evaluation for given V_ and RK4 step
		R = sin(V_(1:2*N,1)) + mu_m.*(((V_(1:2*N,1)./thet0).^2) - 1).*V_(1:2*N,2);
		cos_Phi_Psi = cos(V_(1:N,1));
		sum1Phi = dot(cos_Phi_Psi, cos_Phi_Psi);
		cos_Phi_Psi = cos(V_(N+1:2*N,1));
		sum1Psi = dot(cos_Phi_Psi, cos_Phi_Psi);

		sum2Phi = sum(R(1:N).*cos(V_(1:N,1)) + sin(V_(1:N,1)).*(V_(1:N,2).^2));
		sum2Psi = sum(R(N+1:2*N).*cos(V_(N+1:2*N,1)) + sin(V_(N+1:2*N,1)).*(V_(N+1:2*N,2).^2));

		F(2*N+1,2,rk) = (-(chi + omeg_r2)*V_(2*N+1,1) - mu_s*V_(2*N+1,2) + chi*V_(2*N+2,1) + sum2Phi)/(1 - beta*sum1Phi);
		F(2*N+2,2,rk) = (-(chi + omeg_r2)*V_(2*N+2,1) - mu_s*V_(2*N+2,2) + chi*V_(2*N+1,1) + sum2Psi)/(1 - beta*sum1Psi);
	
		F(1:(2*N+2),1,rk) = V_(1:2*N+2,2);
		F(1:N,2,rk) = -(R(1:N) + beta.*cos(V_(1:N,1).*F(2*N+1,2,rk)));
		F(N+1:2*N,2,rk) = -(R(N+1:2*N) + beta.*cos(V_(N+1:2*N,1)).*F(2*N+2,2,rk));
		
		%executing the RK4 steps
		if rk<3
			V_ = V(:,:,t) + (timestep/2)*F(:,:,rk);
		elseif rk==3
			V_ = V(:,:,t) + timestep*F(:,:,rk);
		else
			V(:,:,t+1) = V(:,:,t) + (timestep/6)*(F(:,:,1) + 2*(F(:,:,2)+F(:,:,3)) + F(:,:,4));
		end
	end
end

%re-scaling the swing oscillations
V(2*N+1:2*N+2,:,:) = V(2*N+1:2*N+2,:,:)*(x0/L);

figure(1)
for x = 1:N
	plot(timeVec/omega,reshape(V(x,1,:),[1,timeCount]))
        hold on
end
hold off
title("All Phi's")
xlabel("Slow Time (s)")
ylabel("Metronome Angle (rad)")

figure(2)
for x = 1:N
	plot(timeVec/omega,reshape(V(N+x,1,:),[1,timeCount]))
        hold on
end
hold off
title("All Psi's")
xlabel("Slow Time (s)")
ylabel("Metronome Angle (rad)")

figure(3)
plot(timeVec/omega,reshape(V(2*N+1,1,:),[1,timeCount]))
hold on
plot(timeVec/omega,reshape(V(2*N+2,1,:),[1,timeCount]))
hold off
title("Swing angles")
xlabel("Slow Time (s)")
ylabel("Deviation from Equilibrium Swing Angle (rad)")