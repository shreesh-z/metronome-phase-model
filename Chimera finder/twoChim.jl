using LinearAlgebra
using PyPlot

N = 2;
freqBpm = 200;
m = 0.028;
M = 2.31;
l = 0.15;
L = 0.22;
g = 9.81;

l_bob = 0.073 - 0.00022*freqBpm;
r_cm = abs(0.178*l_bob - 0.0121);
I = 0.0000129 + 0.005*(l_bob^2);

x0 = (m*r_cm)/M;
omega = sqrt((m*g*r_cm)/I);
Omega = sqrt(g/L);

mu_m = 0.011;         #nonlinearity of metros
mu_s = 0.00016;      #damping of swings
thet0 = 0.33;

#=kappa = 50;        #coupling of swings
freq = freqBpm/60;       #freq of metros in Hz
omega = 2*pi*freq;
x0 = 0.0000091;=#

beta = (x0*(omega^2))/g;
omeg_r2 = (Omega/omega)^2;

beta = 0.0005;
omeg_r2 = 0.6;
thet0 = 1;

println("x0        : ",x0)
println("omega     : ",omega)
println("Free TimP : ",(2*pi)/omega)
println("Beta      : ",beta)
println("Omega^2_r : ",omeg_r2)

timestep = 0.01;
maxTime = 100;
repe = 500;
oscill = round(maxTime/(2*pi)) + 2;
println("Num of oscills : ",maxTime/(2*pi));

timeVec = collect(0:timestep:maxTime);
timeCount = length(timeVec);
timeVec = timeVec .+ (repe*maxTime);
println("Count is: ",timeCount);
global V = zeros(Float64,N+1,2,timeCount);
R = zeros(Float64,N)
F = zeros(Float64,N+1,2,4);
global V_ = zeros(Float64,N+1,2);

#V[1,1,1] = 2*thet0;
#V[2,1,1] = -2*thet0
V[:,:,1] = [1.6109 1.6338;
			-1.8599 -1.8410;
			-1  -2]
#V[1:N,1,1] .= 2*thet0
#V[2,1,1] = 0.0036*(L/x0) 

#limit cycle points
#V[:,:,1] = [0.07051553525115307 -0.6454509062924059; -0.11848332166894515 1.0900217522545772];
#init = zeros(Float64,N+1,2);
#init = [-0.6546581804868434 -0.09211355761541334; -0.6546581804868434 -0.09211355761541334; 2.1217870007739483 0.2771296047953938];
#init = [0.43069649291201895 0.49291497643815513; 0.43069649291201895 0.49291497643815513; 0.43069649291201895 0.49291497643815513; -0.0008221847429328909 -0.0008926953491499909];
#init = [-0.11297146804226176 0.6399500440993542; -0.11297146804226176 0.6399500440993542; -0.11297146804226176 0.6399500440993542; 0.00019788893098617933 -0.0012310166614600238];
#init = [-0.5856060936808007 0.2969393436000161; -0.5856060936808007 0.2969393436000161; -0.5856060936808007 0.2969393436000161; 0.001097076044939959 -0.0005199115426878482];
#V[:,:,1] = init

function findFunct(V_,rk)
	R = sin.(V_[1:N,1]) .+ mu_m.*(((V_[1:N,1]./thet0).^2).-1).*V_[1:N,2];
	sum1Phi = sum(cos.(V_[1:N,1]).^2);
	sum2Phi = sum(R.*cos.(V_[1:N,1]) .+ sin.(V_[1:N,1]).*(V_[1:N,2].^2));
	
	F[1:N+1,1,rk] = V_[1:N+1,2];
	F[N+1,2,rk] = (-omeg_r2*V_[N+1,1] - mu_s*V_[N+1,2] + sum2Phi)/(1 - beta*sum1Phi);
	F[1:N,2,rk] = -(R .+ beta.*cos.(V_[1:N,1]).*F[N+1,2,rk]);
	
	return nothing; 
end	

#=global vecDiff1 = 0;
global vecDiff2 = 0;
global timestep2 = timestep^2
global limitCount = 0 =#
for i = 1:repe
	for t = 1:(timeCount-1)
		V_ = V[:,:,t]
		for rk = 1:4
			findFunct(V_,rk);
			if rk<3
				V_ = V[:,:,t] + (timestep/2)*F[:,:,rk]
			elseif rk==3
				V_ = V[:,:,t] + timestep*F[:,:,rk]
			else
				V[:,:,t+1] = V[:,:,t] + (timestep/6)*(F[:,:,1] + 2*(F[:,:,2]+F[:,:,3]) + F[:,:,4])
			end
		end
		#=global vecDiff1 = vecDiff2;
		global vecDiff2 = sum((V[1,:,t+1].-init[1,:]).^2);
		if(vecDiff2<=timestep2)
			if(vecDiff2<vecDiff1)
				global limitCount = t+1;
				break
			end
		end  =#
	end
	V[:,:,1] = V[:,:,end];
end

println(V[:,:,end])
#println(init)
#limitCount = 6400;
#println("Time period is: ",(limitCount*timestep)/omega)
#println("Limit Count is: ",limitCount)
#println("Time period in slow time: ",limitCount*timestep)
#[0.07051553525115307 -0.6454509062924059; -0.11848332166894515 1.0900217522545772]
#this lies on the limit cycle

V[N+1,:,:] = V[N+1,:,:]*(x0/L);

figure(1)
for i = 1:N
	plot(timeVec/omega,V[i,1,:])
	#plot(timeVec[1:limitCount]/omega,V[i,1,1:limitCount])
end
#plot(timeVec[1:limitCount]/omega,V[2,1,1:limitCount])
#plot(timeVec/omega,V[1,1,:])
#plot(timeVec/omega,V[2,1,:])
title("Metro oscills")

figure(2)
#plot(timeVec[1:limitCount]/omega,V[N+1,1,1:limitCount])
plot(timeVec/omega,V[N+1,1,:])
title("Swing oscills")

figure(3)
for i = 1:N
	plot(V[i,1,:],V[i,2,:])
	#plot(V[i,1,1:limitCount],V[i,2,1:limitCount])
end
#plot(V[1,1,:],V[1,2,:])
#plot(V[2,1,:],V[2,2,:])
title("Metro Phase space")

figure(4)
#plot(V[N+1,1,1:limitCount],V[N+1,2,1:limitCount])
plot(V[N+1,1,:],V[N+1,2,:])
title("Swing phase space") 
show()