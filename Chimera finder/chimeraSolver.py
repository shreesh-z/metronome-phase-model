import numpy as np
import matplotlib.pyplot as plt
import time as timeLib

#this is the chimera solver

randomICs = {
    1: np.array([[ 0.14157382,  0.2614215 , -0.29608024, -0.24948547, -0.07842186, -0.31545333,  0.06066683,  0.27495881,  0.32353985, -0.19375554, 0.07076737,  0.12151782, -0.28229673,  0.14350184,  0.03906534],[0.04105613, -0.22835415,  0.05830795, -0.20626723, -0.15800162, -0.21524905,  0.18515062,  0.05557701, -0.13233109, -0.31297807, -0.01390028,  0.28235838,  0.0213475 , -0.20703707, -0.04162385]]),
    2: np.array([[-0.15290063,  0.3241958 ,  0.1120836 , -0.11188795, -0.30505519, -0.28341886,  0.14095393, -0.11355424,  0.01730479, -0.09568901, -0.1053118 ,  0.22120414, -0.17687524, -0.2194714 , -0.26651887],[0.04002332, -0.31970233, -0.22865928, -0.03622754, -0.28424214, 0.26495288,  0.08682372, -0.27123935, -0.08321772,  0.03173398, -0.13748985, -0.12752627, -0.28184933, -0.19892204, -0.07679633]]),
    3: np.array([[ 0.21701744,  0.18755557,  0.09966056, -0.28854996, -0.17635286, -0.01643674,  0.1727202 ,  0.0847469 , -0.08738345, -0.03049813, -0.12964068, -0.15500856, -0.09896353, -0.01853356, -0.10522117],[-0.00344076, -0.0126174 , -0.18630836,  0.03061491, -0.12380928, -0.08615331,  0.24836217, -0.09624896,  0.06622524, -0.31083154, -0.25894798,  0.04698623,  0.2034848 , -0.30503579, -0.04303378]]),
    4: np.array([[-0.25608729,  0.09047053,  0.09671974,  0.03058346, -0.07007772, -0.2251715 ,  0.0924076 , -0.13879626, -0.31291416,  0.04675506, -0.21283704,  0.27935477,  0.23249934,  0.17652559, -0.15843773],[0.04660167,  0.20153208,  0.29906662, -0.09469651,  0.25854828, -0.2800025 , -0.13962833, -0.13480114,  0.17154096,  0.27852486, 0.26370911,  0.24858714, -0.32204741,  0.16938723, -0.05016757]]),
    5: np.array([[ 0.20390167,  0.24951419, -0.26826464,  0.2513322 , -0.06521685, -0.1472395 , -0.29960625,  0.15699906,  0.10392526, -0.02591696, 0.18069027,  0.10559009,  0.11846373, -0.14896994,  0.19258071],[0.01096508,  0.06217605,  0.29208454,  0.20445947,  0.26906981, 0.22516534,  0.30244717, -0.05703924,  0.08187052, -0.12740978, 0.08876159,  0.23601612,  0.03232899,  0.2630207 , -0.17012485]]),
    6: np.array([[-0.28236782,  0.02272068,  0.16583858,  0.20703668, -0.03638689, 0.31897969, -0.25031005, -0.20690587,  0.10063275, -0.08376842, 0.08703523,  0.28748882, -0.21245525, -0.02141477, -0.06271233],[0.27097263,  0.23329002, -0.05803481,  0.2827175 ,  0.23302832, 0.23946649,  0.06657793,  0.31238389, -0.07159775,  0.02263875, -0.01246499, -0.19855326,  0.06315608, -0.22265275,  0.11808502]]),
    7: np.array([[-0.05060081, -0.2054146 , -0.06674946,  0.03285679,  0.2026197 , -0.03184899, -0.27685172,  0.0690691 ,  0.19276587,  0.03763812, -0.16118089,  0.20861508,  0.16564133,  0.16393983, -0.23569675],[0.13071014, -0.19924909,  0.17357686, -0.22026212,  0.29915596, -0.32026914,  0.28966824,  0.24144697,  0.18708392,  0.14308114, -0.00368212, -0.28380399, -0.00924944,  0.00729014, -0.14110797]]),
    8: np.array([[ 0.1070623 , -0.03246503, -0.28466685, -0.13138537,  0.05092291, 0.21429089, -0.1394381 ,  0.04477643, -0.16074862, -0.05627805, -0.02095395, -0.1479139 ,  0.32072046, -0.19726982, -0.20698106],[-0.28029995, -0.29186344, -0.32068636,  0.03183383,  0.23346458, -0.15516512, -0.27783551, -0.28181846, -0.15027134, -0.18853984, -0.04635679, -0.00600303,  0.21760392, -0.04373139, -0.0743389]]),
    9: np.array([[ 0.30193901, -0.22791099,  0.25686854,  0.05015035,  0.07947688, 0.01271234, -0.00403621, -0.10561666,  0.13115458,  0.09444678, -0.16631024,  0.14466168,  0.27436442,  0.25167074, -0.30699128],[0.05734172,  0.28317239,  0.25334837, -0.17983101,  0.17483002, 0.31991728, -0.17518918, -0.16270593,  0.08318364,  0.2243455 , -0.07216568, -0.2861955 , -0.23550847, -0.07095959,  0.24982223]]),
    10: np.array([[ 0.04646099, -0.14375112, -0.22678913,  0.07848184, -0.13759088, -0.06037901,  0.10562991,  0.19844441, -0.10766317, -0.25007225, 0.25419647,  0.08081354,  0.29445089,  0.20517673,  0.27329605],[0.26712177,  0.32415907, -0.09147037, -0.12900987,  0.04657682, 0.19501185,  0.00856151,  0.30200478, -0.24431119, -0.02569469, -0.19062023, -0.2423523 , -0.31037919,  0.04990117,  0.08438664]])
}
"""
rando = 0
while(rando not in range(1,11)):
    rando = (int)(input("Which set of ICs? (1-10): "))
"""
#all physical params

N = 15
freqBpm = 160
m = 0.028
M = 2.31
l = 0.15
L = 0.22
k = 68
g = 9.81
pi = np.pi

l_bob = 0.073 - 0.00022*freqBpm
r_cm = abs(0.178*l_bob - 0.0121)
I = 0.0000129 + 0.005*(l_bob**2)

#derived params

x0 = (m*r_cm)/M
omega = np.sqrt((m*g*r_cm)/I)
kappa = (k/M)*((l/L)**2)
Omega = np.sqrt(g/L)
#Omega = 4.44

#all the constants:

mu_m = 0.011         #nonlinearity of metros
mu_s = 0.00016      #damping of swings
thet0 = 0.33

#all the params that can be customized
"""
kappa = 50        #coupling of swings
freq = freqBpm/60       #freq of metros in Hz
omega = 2*pi*freq
x0 = 0.0000091
"""

#nondimensional parameters
beta = (x0*(omega**2))/g
omeg_r2 = (Omega/omega)**2
chi = kappa/(omega**2)
#chi = 0.092
#omeg_r2 = 0.6
#beta = 0.0005
print("x0        : ",x0)
print("omega     : ",omega)
print("kappa     : ",kappa)
print("Beta      : ",beta)
print("Omega^2_r : ",omeg_r2)
print("Chi       : ",chi)
    
timestep = 0.01
maxTime = 1000
oscill = round(maxTime/(2*pi)) + 2
print("Num of oscills : ",(maxTime/(2*pi)))
print("Expected time  : ",(int)(0.000225425*(maxTime/timestep)),"seconds")
    
time = np.arange(0,maxTime,timestep)
count = time.shape[0]
print(count)
V = np.zeros((2*N+2,2,count))   #matrix with all variables (tamil tha)
R = np.zeros(2*N)               #matrix with re-used terms (tamil na) except PHI^.. and PSI^..
#F = np.zeros((2*N+2,2))         #matrix to hold function values
F_ = np.zeros((2*N+2,2,4))      #matrix used to hold intermediate RK4 values
V_ = np.zeros((2*N+2,2))        #used to hold vectors of intermediate RK4 steps

#flip = np.zeros((N,2))              #to check where flipping of metro angle has happened
#flipIndex = np.zeros((N,2*oscill))  #to store times at which flips have happened
#indices = np.zeros(N,'int')         #to keep track of which index of flipIndex is to be filled next
#omAvg = np.zeros((N,2*oscill-2))    #stores all avg frequencies

#flippy = np.zeros(2)               #
#flippyIndex = np.zeros(2*oscill)
#flippyIndexCount = 0
#omAvgPhi = np.zeros(2*oscill-2)
#P_ = np.zeros(2)

#initial conditions

#initPhi = 2*thet0
#V[0:N,0,0] = np.zeros(N) + initPhi     #for AP IC
#V[N:2*N,0:2,0] = np.transpose(randomICs.get(rando))[0:N,:]
V[N:2*N,0:2,0] = 2*thet0*np.random.rand(N,2) - thet0
avAng = np.sum(V[N:2*N,0,0])
avVel = np.sum(V[N:2*N,1,0])
V[0:N,0,0] = avAng
V[0:N,1,0] = avVel
#initPshi = np.array([0,0])
#V[2*N:2*N+2,0,0] = initPshi

print("Init conditions of psi are:")
print(V[N:2*N,0,0])

t_init = timeLib.time()

for t in range(count-1):
    V_ = V[:,:,t]
    for rk in range(4):
        #devanagari letters 'ka' and 'cha':
        R = np.sin(V_[0:2*N,0]) + mu_m*((V_[0:2*N,0]/thet0)**2 - 1)*V_[0:2*N,1]
        
        cos_Phi = np.cos(V_[0:N,0])
        sum1Phi = cos_Phi.dot(cos_Phi)  #the devanagari letter 'a'
        cos_Psi = np.cos(V_[N:2*N,0])
        sum1Psi = cos_Psi.dot(cos_Psi)  #the devanagari letter 'i'
        
        oneArr = np.ones((1,N))
        sum2Phi = oneArr.dot(R[0:N]*np.cos(V_[0:N,0]) + np.sin(V_[0:N,0])*(V_[0:N,1]**2))[0]
        sum2Psi = oneArr.dot(R[N:2*N]*np.cos(V_[N:2*N,0]) + np.sin(V_[N:2*N,0])*(V_[N:2*N,1]**2))[0]
        
        F_[2*N,1,rk] = ( -(chi + omeg_r2)*V_[2*N,0] - mu_s*V_[2*N,1] + chi*V_[2*N+1,0] + sum2Phi)/(1 - beta*sum1Phi)
        F_[2*N+1,1,rk] = ( -(chi + omeg_r2)*V_[2*N+1,0] - mu_s*V_[2*N+1,1] + chi*V_[2*N,0] + sum2Psi)/(1 - beta*sum1Psi)
        
        F_[0:2*N+2,0,rk] = V_[0:2*N+2,1]
        F_[0:N,1,rk] = -(R[0:N] + beta*np.cos(V_[0:N,0])*F_[2*N,1,rk])
        F_[N:2*N,1,rk] = -(R[N:2*N] + beta*np.cos(V_[N:2*N,0])*F_[2*N+1,1,rk])
        if rk<2:
            V_ = V[:,:,t] + (timestep/2)*F_[:,:,rk]
        elif rk==2:
            V_ = V[:,:,t] + timestep*F_[:,:,rk]
    
    V[:,:,t+1] = V[:,:,t] + (timestep/6)*(F_[:,:,0] + 2*(F_[:,:,1]+F_[:,:,2]) + F_[:,:,3])      #RK4 step 
    
#checking for flips in phi population (since theyre in sync I only check first phi metro
#flippy = V[0,0,t:t+2]
#flippy = np.heaviside(flippy,0.5)
#flippy[1] = np.absolute(flippy[1] - flippy[0])
#if(flippy[1] != 0):
#    flippyIndex[flippyIndexCount] = (t*timestep)/omega
#    flippyIndexCount += 1

#checking for flips in psi population
#flip[:,0:2] = V[N:2*N,0,t:t+2]
#flip = np.heaviside(flip,0.5)
#flip[:,1] = np.absolute(flip[:,1] - flip[:,0])
#for i in range(N):
#    if(flip[i,1] != 0):
#        hindex = indices[i]
#        flipIndex[i,hindex] = (t*timestep)/omega
#        indices[i] += 1
    
V[2*N:2*N+2,:,:] = V[2*N:2*N+2,:,:]*(x0/L)

t_elapsed = timeLib.time() - t_init
print("Elapsed time   : ",t_elapsed)

#print(" ")
#print("At maxtime: ",V[N:2*N,0,count-1])
#print("At half of maxtime: ", V[N:2*N,0,round(count/2)-1])
#print("At 7 steps: ", V[N:2*N,0,7],"\n")

#POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING#
#__#POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING#
#____#POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING#
#______#POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING____POST_PROCESSING#

#finding average frequencies
#avIndex = 0
#while(avIndex < (2*oscill)-2):
#    omAvg[:,avIndex] = 1/(flipIndex[:,avIndex+2] - flipIndex[:,avIndex])
#    omAvgPhi[avIndex] = 1/(flippyIndex[avIndex+2] - flippyIndex[avIndex])
#    avIndex += 1

#finding angle difference between consecutive metronomes
#diff = np.zeros((N-1,count))
#for i in range(N-1):
#    diff[i,:] = V[N+i+1,0,:] - V[N+i,0,:]

#calculating the Z parameter
#since the \phi metronomes are always in sync
#the theta_sync is taken to be angle of first phi metro (refer to definition of Z)
"""
Z1 = np.zeros(count, 'complex128')
Z2 = np.zeros(count, 'complex128')
thetSync = np.zeros(count)
thetSync = V[0,0,:]
for i in range(N):
    Z1 += np.exp(1j*(V[i,0,:] + thetSync))
    Z2 += np.exp(1j*(V[N+i,0,:] + thetSync))

Z1 = Z1/N
Z2 = Z2/N

#Z1 = np.sum(np.exp(1j*(V[0:N,0,:] + V[0,0,:])),axis=0)/N
#Z2 = np.zeros(count)
#Z2 = np.sum(np.exp(1j*(V[N:2*N,0,:] + V[0,0,:])),axis=0)/N

V_avg = np.zeros((2,count))   #to plot average angles of the two populations
V_avg[0,:] = np.sum(V[0:N,0,:],axis=0)
V_avg[1,:] = np.sum(V[N:2*N,0,:],axis=0)
V_avg/=N
#for x in range(N):
#    V_avg[0,:] += V[x,0,:]
#    V_avg[1,:] += V[N+x,0,:]
    
#V_avg/=N
"""
#plotting all right metronome angles
plt.figure(1)
for x in range(N):
    plt.plot(time/omega,V[x,0,:])
#plt.plot(time/omega,V_avg[0,:])
plt.title("All \phi")

#plotting all left metronome angles    
plt.figure(2)
for x in range(N):
    plt.plot(time/omega,V[N+x,0,:])
#plt.plot(time/omega,V[N,0,:])
plt.grid()
#plt.plot(time/omega,V_avg[1,:])
plt.title("All \psi")

#plotting swing angles
plt.figure(3)
plt.plot(time/omega,V[2*N,0,:])
plt.plot(time/omega,V[2*N+1,0,:])
plt.title("Swing angles")
"""
plt.figure(4)
plt.plot(time/omega,V_avg[1,:]-V_avg[0,:])
plt.title("Average \psi - Average \phi")

plt.figure(5)
plt.plot(time/omega,np.absolute(Z1))
plt.plot(time/omega,np.absolute(Z2))
#plt.plot(Z2.real,Z2.imag)
#plt.plot(Z1.real,Z1.imag)
plt.title("mod(Zp) parameter")
"""
#plt.figure(6)
#for i in range(N-1):
#    plt.plot(time/omega,diff[i,:])
#plt.title("Diff of consecutive angles")

#plt.figure(6)
#for i in range(N):
#    plt.plot(flipIndex[i,0:-3]/omega,omAvg[i,0:-1]/omAvgPhi[0:-1])
#plt.plot(flipIndex[0,0:-3]/omega,omAvg[0,0:-1])
#plt.title("Average angular frequencies")
plt.show()
exit()