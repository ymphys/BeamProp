me = 0.511
mp = 938
mH = 938.5
c = 2.998e8

def calc_v(m,Ek):
    """
    calculate the velocity of a particle with mass m (in MeV/c^2) and kinetic energy Ek (in MeV)
    """
    gamma = Ek/m + 1
    beta = (1 - 1/gamma**2)**0.5
    return gamma,beta

import matplotlib.pyplot as plt
import numpy as np
Ek_lst = np.linspace(50,500,450)
beta_lst = []
t_lst = []
for Ek in Ek_lst:
    _,beta = calc_v(mp,Ek)
    t = 1e3/(beta*c)*1e6 # in us
    beta_lst.append(beta)
    t_lst.append(t)
plt.plot(Ek_lst,beta_lst)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Beta (v/c)')
plt.title('Proton Velocity vs Kinetic Energy')
plt.legend()
plt.savefig("proton_v_Ek.png")
plt.close()

plt.plot(Ek_lst,t_lst)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('prop 1km time used (us)')
plt.title('Proton Prop. time in 1 km vs Kinetic Energy')
plt.legend()
plt.savefig("proton_t_Ek.png")
plt.close()