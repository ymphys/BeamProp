# Simulate again for KE = 30 MeV but extend propagation distance to 20 m.
import numpy as np
from math import sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Physical constants
eps0 = 8.8541878128e-12
c = 299792458.0
qe = 1.602176634e-19
me = 9.10938356e-31
rc = qe**2/(4*np.pi*eps0*me*c**2)  # classical radius (electron)

# Beam / simulation parameters (defaults)
N = 1e9             # number of electrons in the bunch
KE_MeV = 30.0       # kinetic energy in MeV
gamma0 = 1.0 + KE_MeV/0.511   # relativistic factor
beta0 = np.sqrt(1.0 - 1.0/gamma0**2)

# Geometry / initial beam sizes
a0 = 1.0e-2         # initial transverse semi-axis (m)
zm0 = 5.0e-2        # initial half-length (m)
Ib= 3/4 * N * qe * beta0 * c / zm0   # beam current in Amperes
# print(f"IB ={Ib:.2f}")
# Normalized emittances (same as before)
eps_nx = 1.0e-6     # normalized transverse emittance (m·rad)
eps_nz = 1.0e-6     # normalized longitudinal normalized emittance

# Convert to un-normalized emittances used in dynamic equations:
eps_x = eps_nx/(beta0*gamma0)
eps_zzp = eps_nz/(beta0*(gamma0**3))

# Geometry factor g
# g = 1.0
xi = np.sqrt(1-a0**2/zm0**2)
g = 2/xi**2*(1/2/xi*np.log((1+xi)/(1-xi)) - 1) 

# Space-charge prefactors for transverse and longitudinal ODEs (from equations)
prefac_common = 3.0/2.0 * N * rc / (beta0**2 * gamma0**3)
prefac_long = 3.0/2.0 * N * rc / (beta0**2 * gamma0**5)

# ODE system
def deriv(s, y):
    a, ap, zm, zmp = y
    a = max(a, 1e-12)
    zm = max(zm, 1e-12)
    sc_trans = prefac_common * (1.0/(a*zm)) * (1.0 - (g/2.0)*(a**2/(gamma0**2 * zm**2)))
    sc_long = prefac_long * (g / (zm**2))
    app = sc_trans + (eps_x**2)/(a**3)
    zmpp = sc_long + (eps_zzp**2)/(zm**3)
    return [ap, app, zmp, zmpp]

# Integration range
s_max = 20000.0   # meters (extended)
s_eval = np.linspace(0, s_max, 4000000)
y0 = [a0, 0.0, zm0, 0.0]

sol = solve_ivp(deriv, [0, s_max], y0, t_eval=s_eval, method='RK45', rtol=1e-6, atol=1e-9)

a_sol = sol.y[0]
zm_sol = sol.y[2]

plt.figure(figsize=(8,4))
plt.plot(sol.t, a_sol)
plt.xlabel("Propagation distance s (m)")
plt.ylabel("Transverse semi-axis a(s) (m)")
plt.title(f"Transverse expansion a(s) — KE = {KE_MeV} MeV, s up to {s_max} m")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"vacuum-diffusion/transverse_expansion_{KE_MeV}MeV_{s_max}m.png")
plt.close()
plt.figure(figsize=(8,4))
plt.plot(sol.t, zm_sol)
plt.ylim(zm0 * 0.9, zm0 * 1.1)  # Limit y-axis to avoid excessive range
plt.plot(sol.t, zm_sol)
plt.ylim(zm0 * 0.9, zm0 * 1.1)  # Limit y-axis to avoid excessive range
plt.xlabel("Propagation distance s (m)")
plt.ylabel("Longitudinal half-length z_m(s) (m)")
plt.title(f"Longitudinal expansion z_m(s) — KE = {KE_MeV} MeV, s up to {s_max} m")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"vacuum-diffusion/longitudinal_expansion_{KE_MeV}MeV_{s_max}m.png")
plt.close()

final_a = a_sol[-1]
final_zm = zm_sol[-1]
print(f"Parameters used: Ib={Ib:.3g} A, KE={KE_MeV} MeV")
print(f"Initial a0={a0:.3e} m, zm0={zm0:.3e} m, eps_nx={eps_nx:.1e}, eps_nz={eps_nz:.1e}")
print(f"Final a(s={s_max} m) = {final_a:.3e} m")
print(f"Final zm(s={s_max} m) = {final_zm:.3e} m")
