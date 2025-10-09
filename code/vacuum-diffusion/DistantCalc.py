import numpy as np
from math import pi

# Define the beam paramters
# Take input from the users:
I = input("请输入电流强度，单位是安培：")
R = input("请输入束流半径，单位是厘米：")
E0 = input("请输入束流平均能量，单位是MeV：") 

I, R, E0 = float(I), 0.01*float(R), float(E0)

# Define physical constant
c, eps0, me, e0 = 3.0e8, 8.85e-12, 9.11e-31, 1.6e-19
MeV2J = 1.6e-13
# c: speed of light in the vacuum[m/s], eps0: vacuum dielectric constant[C^2/(N*m^2)], me: electron mass[kg], e0: elementary charge[C]
# MeV2J: 1MeV = 1.6e-13 J

# Calculate velocity, charge density and inherent relaxation time
gamma = E0*MeV2J/me/c**2
beta = np.sqrt(gamma**2-1)/gamma
v = beta*c
rho = I/(pi*R**2*v)
rho0 = rho/gamma
tau0 = np.sqrt(2*eps0*me/e0/rho0)

def distance(r):
    rR = r/R
    tau = tau0*(-2.81e-4*rR**3+3.27e-2*rR**2+1.77*rR-1.18)
    L = gamma*v*tau
    return L/1000

r = input("请输入扩散后的束流半径，单位是厘米：")
r = 0.01*float(r)
print(f"可传输距离为 {distance(r):.2f} km")

