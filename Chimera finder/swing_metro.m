%all physical params

N = 2; %N can be kept 1 also
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

disp("x0        : "+x0)
disp("omega     : "+omega)
disp("kappa     : "+kappa)
disp("Beta      : "+beta)
disp("Omega^2_r : "+omeg_r2)

timestep = 0.01;
maxTime = 50;
oscill = round(maxTime/(2*pi)) + 2;

disp("Num of oscills : "+maxTime/(2*pi));

timeVec = 0:timestep:maxTime;
timeCount = length(timeVec);
disp("Count is: "+timeCount);

V = zeros(N+1,2,timeCount);
R = zeros(N);
F = zeros(N+1,2,4);
V_ = zeros(N+1,2);
cos_Phi_Psi = zeros(N,1);

V(1:N,1,1) = 2*thet0;   %for IP IC

for t = 1:(timeCount-1)
	V_ = V(:,:,t);
	for rk = 1:4
		
		%function evaluation
		R = sin(V_(1:N,1)) + mu_m.*(((V_(1:N,1)./thet0).^2) - 1).*V_(1:N,2);
		sum1Phi = sum(cos(V_(1:N,1)).^2);
		sum2Phi = sum(R.*cos(V_(1:N,1)) + sin(V_(1:N,1)).*(V_(1:N,2).^2));
		
		F(1:N+1,1,rk) = V_(1:N+1,2);
		F(N+1,2,rk) = (-omeg_r2*V_(N+1,1) - mu_s*V_(N+1,2) + sum2Phi)/(1 - beta*sum1Phi);
		F(1:N,2,rk) = -(R + beta.*cos(V_(1:N,1)).*F(N+1,2,rk));
		
		%executing RK4 steps
		if rk<3
			V_ = V(:,:,t) + (timestep/2)*F(:,:,rk);
		elseif rk==3
			V_ = V(:,:,t) + timestep*F(:,:,rk);
		else
			V(:,:,t+1) = V(:,:,t) + (timestep/6)*(F(:,:,1) + 2*(F(:,:,2)+F(:,:,3)) + F(:,:,4));
		end
	end
end

%re-scaling the swing angles
V(N+1,:,:) = V(N+1,:,:)*(x0/L);

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
plot(timeVec/omega,reshape(V(N+1,1,:),[1,timeCount]))
title("Swing angles")
xlabel("Slow Time (s)")
ylabel("Deviation from Equilibrium Swing Angle (rad)")