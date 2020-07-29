using LinearAlgebra
using PyPlot
#using BenchmarkTools

#all physical params

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

#derived params
#29.43722
x0 = (m*r_cm)/M;
omega = sqrt((m*g*r_cm)/I);
kappa = (k/M)*((l/L)^2);
Omega = sqrt(g/L);

#all the constants:

mu_m = 0.011;         #nonlinearity of metros
mu_s = 0.00016;      #damping of swings
thet0 = 0.33;

#all the params that can be customized

kappa = 20;        #coupling of swings
#freq = freqBpm/60;       #freq of metros in Hz
#omega = 2*pi*freq;
x0 /= 100; 


#nondimensional parameters
beta = (x0*(omega^2))/g;
omeg_r2 = (Omega/omega)^2;
chi = kappa/(omega^2);

println("x0        : ",x0)
println("omega     : ",omega)
println("kappa     : ",kappa)
println("Beta      : ",beta)
println("Omega^2_r : ",omeg_r2)
println("Chi       : ",chi)

timestep = 0.01;
maxTime = 2500;
oscill = round(maxTime/(2*pi)) + 2;

println("Num of oscills : ",maxTime/(2*pi));
#print("Expected time  : ",convert(Int8, 0.000225425*(maxTime/timestep)))

timeVec = collect(0:timestep:maxTime);
timeCount = length(timeVec);
println("Count is: ",timeCount);
V = zeros(Float64,2*N+2,2,timeCount);
R = zeros(Float64,2*N);
F = zeros(Float64,2*N+2,2,4);
V_ = zeros(Float64,2*N+2,2);
cos_Phi_Psi = zeros(Float64,(N,1))

V[1:N,1,1] .= 2*thet0;   #for AP IC
#V[N+1:2*N,1:2,1] = 2*thet0*rand(Float64,(N,2)) .- thet0;
V[N+1:2*N,1,1] = [0.14157382; 0.2614215; -0.29608024; -0.24948547; -0.07842186; -0.31545333; 0.06066683; 0.27495881; 0.32353985; -0.19375554; 0.07076737; 0.12151782; -0.28229673; 0.14350184; 0.03906534];
V[N+1:2*N,2,1] = [0.04105613; -0.22835415;  0.05830795; -0.20626723; -0.15800162; -0.21524905; 0.18515062; 0.05557701; -0.13233109; -0.31297807; -0.01390028; 0.28235838; 0.0213475; -0.20703707; -0.04162385]
#V[1:N,1,1] .= acos(sum(cos.(V[N+1:2*N,1,1]))/N);
#V[1:N,2,1] .= sqrt(sum(V[N+1:2*N,2,1].^2)/N);

function findFunct(V_,rk)
	R = sin.(V_[1:2*N,1]) .+ mu_m.*(((V_[1:2*N,1]./thet0).^2).-1).*V_[1:2*N,2]
	cos_Phi_Psi = cos.(V_[1:N,1])
	sum1Phi = dot(cos_Phi_Psi, cos_Phi_Psi)
	cos_Phi_Psi = cos.(V_[N+1:2*N,1])
	sum1Psi = dot(cos_Phi_Psi, cos_Phi_Psi)
	
	#sum2Phi = dot(oneArr, R[1:N].*cos.(V_[1:N,1]) .+ sin.(V_[1:N,1]).*(V_[1:N,2].^2))
	#sum2Psi = dot(oneArr, R[N+1:2*N].*cos.(V_[N+1:2*N,1]) .+ sin.(V_[N+1:2*N,1]).*(V_[N+1:2*N,2].^2))
	sum2Phi = sum(R[1:N].*cos.(V_[1:N,1]) .+ sin.(V_[1:N,1]).*(V_[1:N,2].^2))
	sum2Psi = sum(R[N+1:2*N].*cos.(V_[N+1:2*N,1]) .+ sin.(V_[N+1:2*N,1]).*(V_[N+1:2*N,2].^2))
	
	F[2*N+1,2,rk] = (-(chi + omeg_r2)*V_[2*N+1,1] - mu_s*V_[2*N+1,2] + chi*V_[2*N+2,1] + sum2Phi)/(1 - beta*sum1Phi)
	F[2*N+2,2,rk] = (-(chi + omeg_r2)*V_[2*N+2,1] - mu_s*V_[2*N+2,2] + chi*V_[2*N+1,1] + sum2Psi)/(1 - beta*sum1Psi)
	
	F[1:(2*N+2),1,rk] = V_[1:2*N+2,2]
	F[1:N,2,rk] = -(R[1:N] + beta.*cos.(V_[1:N,1]).*F[2*N+1,2,rk])
	F[N+1:2*N,2,rk] = -(R[N+1:2*N] + beta.*cos.(V_[N+1:2*N,1]).*F[2*N+2,2,rk])
	
	return nothing 
end	

for t = 1:(timeCount-1)
	#V_ = reshape(V[:,:,t],(2*N+2,2))
	V_ = V[:,:,t]
	for rk = 1:4
		findFunct(V_,rk)
		if rk<3
			V_ = V[:,:,t] + (timestep/2)*F[:,:,rk]
		elseif rk==3
			V_ = V[:,:,t] + timestep*F[:,:,rk]
		else
			V[:,:,t+1] = V[:,:,t] + (timestep/6)*(F[:,:,1] + 2*(F[:,:,2]+F[:,:,3]) + F[:,:,4])
		end
	end
end

V[2*N+1:2*N+2,:,:] = V[2*N+1:2*N+2,:,:].*(x0/L);
#V[2*N+1,:,:] = V[2*N+1,:,:] .+ 0.15;
#V[2*N+2,:,:] = V[2*N+2,:,:] .- 0.15;

println(sum(V[N+1:2*N,1,end])/N-V[1,1,end]);

#indexx = (timeCount-1)/2
#println("At maxtime: ",V[N+1:2*N,1,convert(Int64,timeCount)])
#println("At half of maxtime: ", V[N+1:2*N,1,convert(Int64,indexx)])
#println("At 7 steps: ", V[N+1:2*N,1,8],"\n")
#=X = zeros(Float64,N,timeCount);
for i = 1:N
	X[i,:] = V[1,1,:];
end
X = reshape(V[N+1:2*N,1,:],N,timeCount).-X
Z2 = abs.(sum(exp.(1im.*X),dims=1))/N;
Z2 = reshape(Z2,timeCount);=#

figure(1)
for x = 1:N
	plot(timeVec./omega,V[x,1,:])
end
title("All Phi's")
xlabel("Slow Time (s)")
ylabel("Metronome Angle (rad)")

figure(2)
for x = 1:N
	plot(timeVec./omega,V[N+x,1,:])
end
title("All Psi's")
xlabel("Slow Time (s)")
ylabel("Metronome Angle (rad)")

figure(3)
plot(timeVec./omega,V[2*N+1,1,:])
plot(timeVec./omega,V[2*N+2,1,:])
title("Swing angles")
xlabel("Slow Time (s)")
ylabel("Deviation from Equilibrium Swing Angle (rad)")

#figure(4)
#plot(timeVec/omega,Z2)
#title("Z param")

show()