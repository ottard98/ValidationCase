import numpy as np
import openfoamparser as ofpp
import matplotlib.pyplot as plt

'''
This script returns Nusselts numbers, pressures and check mass conservation.
'''

L   = 0.06706622
dz  = 0.1*L
t_max = 20
p0      = 101325 #Pa
T0      = 600
Tstar   = 273
R       = 287
Pr      = 0.71
gamma   = 1.4
mustar  = 1.68e-05
S       = 110.5 #Sutherland
eps     = 0.6
g       = 9.81
Ra      = 1e6
Nuref   = 8.85978
pref    = 0.856338

rho0    = p0/(R*T0)
cp      = gamma*R/(gamma - 1)
mu0     = (T0/Tstar)**(3/2)*(Tstar + S)/(T0 + S)*mustar
Th      = T0*(1 + eps)
Tc      = T0*(1 - eps)
k0      = cp*mu0/Pr
dT      = Th-Tc
V       = L**2*dz

pressureField = ofpp.parse_internal_field('%s/p'%t_max)
pMax = max(pressureField)
pMin = min(pressureField)
pAvg = np.average(pressureField)

heatFlux = ofpp.parse_boundary_field('%s/wallHeatFlux'%t_max)
HFhot    = heatFlux[b'hotWall'][b'value']
HFcold   = heatFlux[b'coldWall'][b'value']

rhoAvg_vec = [rho0]
for t in range(t_max):
    t = t+1
    rhoField = ofpp.parse_internal_field('%i/rho' %t)
    rhoAvg_vec.append(np.average(rhoField))

mass_vec= [rho*V for rho in rhoAvg_vec]
massDiff= abs(mass_vec[-1] - mass_vec[0])
Nuh_vec = abs(L*HFhot/(k0*dT))
Nuc_vec = abs(L*HFcold/(k0*dT))
Nuh     = np.average(Nuh_vec)
Nuc     = np.average(Nuc_vec)
diffh   = (Nuh - Nuref)/Nuref*100
diffc   = (Nuc - Nuref)/Nuref*100
pMaxp0  = pMax/p0
pMinp0  = pMin/p0
pAvgp0  = pAvg/p0
pAvgDiff= (pAvgp0-pref)/pref*100
pMaxDiff= (pMaxp0-pref)/pref*100
pMinDiff= (pMinp0-pref)/pref*100

print("mu0    = %f" %mu0)
print("mustar = %f\n" %mustar)
def Ra_calc(Length):
    Ra = Pr*g*rho0**2*(Th - Tc)*Length**3/(T0*mu0**2)
    return Ra

def len_calc(Ra):
    Length = (Ra*T0*mu0**2/(Pr*g*rho0**2*(Th - Tc)))**(1/3)
    return Length

print("For Ra = %.0f --> L = %.8f" %(Ra,len_calc(Ra)))
print("For L  = %f --> Ra = %.0f\n" %(L,Ra_calc(L)))

print("Nusselts referance   = %.6f" %Nuref)
print("Nusselts hotWall     = %.5f, diff    = %.4f%% " %(Nuh,diffh))
print("Nusselts coldWall    = %.5f, diff    = %.4f%%" %(Nuc,diffc))
print("Nusselts hotWall max = %.5f" %max(Nuh_vec))
print("Nusselts hotWall min = %.5f" %min(Nuh_vec))
print("Nusselts coldWall max= %.5f" %max(Nuc_vec))
print("Nusselts coldWall min= %.5f\n" %min(Nuc_vec))

print("Pavg/P0              = %.6f, diff   = %.4f%%" %(pAvgp0,pAvgDiff))
print("Pmax/p0              = %.6f, diff   = %.4f%%" %(pMaxp0,pMaxDiff))
print("Pmin/p0              = %.6f, diff   = %.4f%%\n" %(pMinp0,pMinDiff))

print("Difference in total mass = %s"%massDiff)




mass_diff_tol_10 = [0.0, -3.037426267536031e-16, -1.946244543565151e-16, 7.079840186330344e-16, 1.4279145899541434e-15, 7.555059551057897e-16, 4.205078125985029e-16, -4.969305132315749e-16, -3.963233278884981e-16, -2.9896874906287785e-17, -1.1322153880676672e-15, -1.1418613992709992e-15, -1.7381014433704572e-15, -1.3622627602783571e-15, -1.5073595041430188e-15, -2.153635478502683e-15, -9.03509716045428e-16, -1.7056499170952505e-15, -2.35059774379362e-15, -1.8116611726418097e-15]
mass_diff_tol_16 = [0.0, 3.8031440518539184e-16, 5.216435465006664e-16, 9.932647152682827e-16, 1.320066966977937e-15, 2.2926437495424806e-15, 3.046343830404724e-15, 3.178480970176395e-15, 3.4996724756434366e-15, 3.4935026876556363e-15, 4.414823812515928e-15, 3.177525517011892e-15, 3.185562165615441e-15, 3.97271665976244e-15, 4.706602945922511e-15, 4.998819148329878e-15, 4.387854283475351e-15, 3.482247313852521e-15, 2.5819529348748704e-15, 2.6745268597461874e-15]
mass_diff_tol_16_50 =[0.0, 1.119305589511975e-14, 1.1689996572746691e-14, 6.015428091624489e-15, -6.099013302859543e-15, -6.067588380516409e-15, -8.633478182579549e-15, -5.327640726585786e-15, 3.788490381866419e-15, 9.118377439828323e-15, 1.2284644194905311e-14, 5.138850635169959e-15, 4.251919047968192e-15, 3.620710095674287e-15, 6.856345460999613e-15, 9.678340756862774e-15, 6.783571778303313e-15, 7.165583637515002e-15, 4.916389290034878e-15, 6.129679283681938e-15]
mass_diff_tol_16_600 = [0.0, 1.8922038415303266e-16, 3.1380537816698417e-16, 1.4393461466102875e-16, 1.5493249044817858e-16, 2.1377756335982934e-16, 5.02270208931066e-16, 2.384329983885075e-16, 6.015052008995908e-16, 6.919073332941478e-16, 6.083424508498275e-16, 1.621830924766754e-16, 2.2027600013116433e-16, 3.1799649718999845e-16, 5.0757602331266694e-17, 1.9484807105459023e-16, -1.0805091088354757e-16, -2.9093210045932905e-16, -5.236188273336634e-16, -7.622822186838241e-16]

t_vec = np.linspace(0,t_max,t_max+1)
mass_difference = [mass_vec[i+1] - mass_vec[1] for i in range(len(mass_vec) - 1)]
plt.plot(t_vec[1:],mass_diff_tol_16_600,color = 'b', label= '600x600, tol = 1e-16')
plt.plot(t_vec[1:],mass_diff_tol_16, color = 'r', label = '300x300, tol = 1e-16')
plt.plot(t_vec[1:],mass_diff_tol_10, color = 'r', label= '300x300, tol = 1e-10', linestyle = '--')

plt.plot(t_vec[1:], mass_diff_tol_16_50,color = 'orange', label= '50x50,     tol = 1e-16')

plt.ylabel('Mass (kg)')
plt.xlabel('Time (s)')
plt.grid()

plt.legend()
plt.show()

print(mass_vec[-1] - mass_vec[1])
