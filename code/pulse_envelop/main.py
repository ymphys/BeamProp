import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# === Physical Constants (SI Units) ===
EPSILON_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
M_E = 9.10938356e-31         # Electron mass (kg)
Q_E = 1.602176634e-19        # Elementary charge (C)
C = 299792458                # Speed of light (m/s)

# === Inputs (SI Units) ===
gamma = 100
Ib = 1.8 # current, (A)
epsilon_n = 1e-6 # normalized emmitance (m-rad) 
theta0 = 5e-6 # beam initial angle (rad)
n_i = 1e5 # background ion density (m^-3)
zmax = 1e5 # total propagation distance (m)

# === Relativistic beta ===
def beta(gamma: float) -> float:
    """
    Calculate relativistic velocity factor beta = v/c.
    Parameters:
        gamma (float): Lorentz factor (must be >= 1)
    Returns:
        float: beta
    """
    if gamma < 1:
        raise ValueError("gamma must be >= 1")
    return np.sqrt(1 - 1 / gamma ** 2)

# === Alfven current ===
def alfven_current(gamma: float) -> float:
    """
    Calculate Alfven current (A).
    Parameters:
        gamma (float): Lorentz factor
    Returns:
        float: Alfven current (A)
    """
    b = beta(gamma)
    return 4 * np.pi * EPSILON_0 * M_E * C ** 3 * b * gamma / Q_E

# === Charge Neutralization Factor ===
def f_e(gamma: float, R: float, n_i: float, Ib: float) -> float:
    """
    Calculate charge neutralization factor f_e.
    Parameters:
        gamma (float): Lorentz factor
        R (float): Beam envelope radius (m)
        n_i (float): Ion density (m^-3)
        Ib (float): Beam current (A)
    Returns:
        float: Charge neutralization factor
    Note: n_i should be in m^-3 (not cm^-3)
    """
    b = beta(gamma)
    return Q_E * b * C * np.pi * R ** 2 * n_i / Ib

# === x factor ===
def x_factor(gamma: float, R: float, n_i: float, Ib: float) -> float:
    """
    Calculate x factor for beam-plasma system.
    Parameters:
        gamma (float): Lorentz factor
        R (float): Beam envelope radius (m)
        n_i (float): Ion density (m^-3)
        Ib (float): Beam current (A)
    Returns:
        float: x factor
    """
    fe = f_e(gamma, R, n_i, Ib)
    return (fe * gamma ** 2 - 1) / (gamma ** 2 - 1)

# === initial beam radius ===
def R_0(epsilon_n: float, theta0: float) -> float:
    """
    calculate initial beam radius for given emmitance and beam initial angle
    Parameters:
        epsilon_n (float): normalized emmitance
        theta0 (float): Beam initial angle (rad)
        n_i (float): Ion density (m^-3)
    Returns:
        float: R_0
    """
    return epsilon_n/beta(gamma)/gamma/theta0

def envelope_ode_1(z, y):
    """
    二阶束流包络微分方程转化为一阶方程组
    y[0]: R(z)
    y[1]: R'(z)
    返回: [R'(z), R''(z)]
    """
    # gamma = params['gamma']
    # n_i = params['n_i']
    # Ib = params['Ib']
    R = y[0]
    R_prime = y[1]
    I_A = alfven_current(gamma)
    xfac = x_factor(gamma, R, n_i, Ib)
    eps = epsilon_n/beta(gamma)/gamma
    R_2prime = -2 * xfac * Ib / (I_A * R) + eps**2 / R**3
    # R_2prime = eps**2 / R**3
    return [R_prime, R_2prime]

if __name__ == "__main__":
    
    R0 = R_0(epsilon_n,theta0)
    y0 = [R0, theta0]
    z_span = (0, zmax)
    z_eval = np.linspace(z_span[0], z_span[1], 2000)

    # 求解二阶微分方程
    sol1 = solve_ivp(lambda z, y: envelope_ode_1(z, y), z_span, y0, t_eval=z_eval)
    Rt1 = sol1.y[0][-1]

    # print out import variables
    print(f"R0={R0:2e}, theta0={theta0:2e},Rt1={Rt1:2e}")

    # 结果可视化
    plt.plot(sol1.t, sol1.y[0])
    plt.xlabel('z (m)')
    plt.ylabel('Envelope radius R (m)')
    plt.title('Beam Envelope Evolution')
    plt.grid()
    plt.show()