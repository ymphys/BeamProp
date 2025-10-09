from modules.beam_params import BeamParams
import numpy as np
from modules import constants
from scipy.integrate import solve_ivp

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
def charge_neutralization_factor(b: float, R: float, n_i: float, Ib: float) -> float:
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
    return constants.Q_E * b * constants.C * np.pi * R ** 2 * n_i / Ib

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
def plasma_equilibrium_radius(params: BeamParams) -> float:
    eps = epsilon_rms(params.Ib, params.gamma, params.em_type)
    denom = 2 * constants.Q_E * params.n_i_m3 * np.pi * params.b * constants.C * params.gamma ** 2
    if denom == 0:
        raise ZeroDivisionError("Denominator in plasma_equilibrium_radius is zero.")
    sqrt_term = np.sqrt(params.Ib ** 2 + denom * (params.gamma ** 2 - 1) * params.I_A * eps**2)
    req = np.sqrt((params.Ib + sqrt_term) / denom)
    if req is not None and req > 0 and np.isreal(req):
        return req
    else:
        raise ValueError("当前输入参数下不存在束流平衡半径！")

def sol_envelope_ode_2nd(params: BeamParams):
    y0 = [params.R0_m, 0.0]
    z_span = (0, params.z_max)
    z_eval = np.linspace(z_span[0], z_span[1], 200)
    def envelope_ode_2nd(z, y):
        R = y[0]
        R_prime = y[1]
        xfac = x_factor(params.gamma, R, params.n_i_m3, params.Ib)
        eps = epsilon_rms(params.Ib, params.gamma, params.em_type)
        R_2prime = -2 * xfac * params.Ib / (params.I_A * R) + eps**2 / R**3
        return [R_prime, R_2prime]
    sol = solve_ivp(envelope_ode_2nd, z_span, y0, t_eval=z_eval)
    R_T = sol.y[0, -1]
    return sol, R_T


