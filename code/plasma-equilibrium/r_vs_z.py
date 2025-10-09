import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# === Physical Constants (SI Units) ===
EPSILON_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
M_E = 9.10938356e-31         # Electron mass (kg)
Q_E = 1.602176634e-19        # Elementary charge (C)
C = 299792458                # Speed of light (m/s)

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

# === Emittance ===
def epsilon_rms(Ib: float, gamma: float, em_type: str) -> float:
    """
    Calculate RMS emittance (m·rad).
    Parameters:
        Ib (float): Beam current (A)
        gamma (float): Lorentz factor
        em_type (str): 'Ti', 'C', or 'LP'
    Returns:
        float: RMS emittance (m·rad)
    """
    if em_type == 'Ti':
        return 0.01 * 0.9 / gamma
    elif em_type == 'C':
        return 0.01 * 0.04 / gamma
    elif em_type == 'LP':
        return 0.01 * 0.005 * np.sqrt(Ib) / gamma
    else:
        raise ValueError("em_type must be 'Ti', 'C', or 'LP'")

# === Charge Neutralization Factor ===
def charge_neutralization_factor(gamma: float, R: float, n_i: float, Ib: float) -> float:
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
    fe = charge_neutralization_factor(gamma, R, n_i, Ib)
    return (fe * gamma ** 2 - 1) / (gamma ** 2 - 1)

# === Plasma Equilibrium Radius ===
def plasma_equilibrium_radius(Ib: float, n_i: float, gamma: float, em_type: str) -> float:
    """
    Calculate plasma equilibrium radius (m).
    Parameters:
        Ib (float): Beam current (A)
        n_i (float): Ion density (m^-3)
        gamma (float): Lorentz factor
        em_type (str): 'Ti', 'C', or 'LP'
    Returns:
        float: Equilibrium radius (m)
    """
    b = beta(gamma)
    I_A = alfven_current(gamma)
    eps = epsilon_rms(Ib, gamma, em_type)
    # Avoid division by zero
    denom = 2 * Q_E * n_i * np.pi * b * C * gamma ** 2
    if denom == 0:
        raise ZeroDivisionError("Denominator in plasma_equilibrium_radius is zero.")
    sqrt_term = np.sqrt(Ib ** 2 + denom * (gamma ** 2 - 1) * I_A * eps**2)
    return np.sqrt((Ib + sqrt_term) / denom)

def envelope_ode_2nd(z, y, params):
    """
    二阶束流包络微分方程转化为一阶方程组
    y[0]: R(z)
    y[1]: R'(z)
    返回: [R'(z), R''(z)]
    """
    gamma = params['gamma']
    n_i = params['n_i']
    Ib = params['Ib']
    em_type = params['em_type']
    R = y[0]
    R_prime = y[1]
    I_A = alfven_current(gamma)
    xfac = x_factor(gamma, R, n_i, Ib)
    eps = epsilon_rms(Ib, gamma, em_type)
    # 你的二阶方程右侧
    R_2prime = -2 * xfac * Ib / (I_A * R) + eps**2 / R**3
    return [R_prime, R_2prime]

# === Example usage ===
# All densities should be in m^-3, not cm^-3. If you have cm^-3, convert: n_i_SI = n_i * 1e6
# gamma = 10
# Ib = 100.0  # A
# n_i = 1e14  # m^-3
# R = 0.01    # m
# em_type = 'Ti'
# print(plasma_equilibrium_radius(Ib, n_i, gamma, em_type))

if __name__ == "__main__":
    try:
        gamma = float(input("Enter Lorentz factor gamma (>=1): "))
        if gamma < 1:
            raise ValueError("gamma must be >= 1")
        n_i_cm3 = float(input("Enter ion density n_i (cm^-3): "))
        n_i = n_i_cm3 * 1e6  # convert to m^-3
        Ib = float(input("Enter beam current Ib (A): "))
        em_type = input("Enter emittance type ('Ti', 'C', or 'LP'): ").strip()
        R_eq = plasma_equilibrium_radius(Ib, n_i, gamma, em_type)
    except Exception as e:
        print(f"Error: {e}")
    params = {
        'gamma': gamma,
        'n_i': n_i,
        'Ib': Ib,
        'em_type': em_type
    }
    # 初始条件
    R0_cm = float(input("Enter initial envelope radius R0 (cm): "))
    R0 = R0_cm * 1e-2  # convert to meters
    R0_prime = float(input("Enter initial envelope slope R0' : "))  # 新增输入
    y0 = [R0, R0_prime]
    z_span = (0, float(input("Enter max z (m): ")))
    z_eval = np.linspace(z_span[0], z_span[1], 200)

    # 求解二阶微分方程
    sol = solve_ivp(lambda z, y: envelope_ode_2nd(z, y, params), z_span, y0, t_eval=z_eval)

    print(f"\nPlasma equilibrium radius R_eq = {R_eq:.6e} m")

    # 结果可视化（将R单位从m转换为cm）
    plt.plot(sol.t, sol.y[0] * 100)  # m -> cm
    plt.xlabel('z (m)')
    plt.ylabel('Envelope radius R (cm)')
    plt.title('Beam Envelope Evolution')
    plt.grid()
    plt.show()