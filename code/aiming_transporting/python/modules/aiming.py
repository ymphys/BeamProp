from modules.beam_params import BeamParams
import numpy as np
from scipy.optimize import root_scalar

# 判断是否超出射程
def check_range(params: BeamParams):
    if params.R_m > 2 * params.a0_max / np.sin(params.Psi):
        raise ValueError(f"目标距离超出射程！当前目标角度下射程需小于: {2*params.a0_max/np.sin(params.Psi)*1e-3:.1f}公里，请选择距离或角度更小的目标后重试")

# 求解发射角theta 
def solve_theta(params: BeamParams):
    def a0(theta):
        return params.a0_max * np.sin(theta)

    def theta_equation(theta):
        a0_val = a0(theta)
        return params.R_m * np.sin(params.Psi) - a0_val * np.sqrt(2 * (1 - np.cos(params.wB * params.R_m * np.cos(params.Psi) / (params.v0 * np.cos(theta)))))

    theta_min, theta_max = 0.01, np.pi/2
    theta_test = np.linspace(theta_min, theta_max, 1000)
    f_test = [theta_equation(theta) for theta in theta_test]
    bracket = None
    for i in range(len(f_test)-1):
        if f_test[i] * f_test[i+1] < 0:
            bracket = [theta_test[i], theta_test[i+1]]
            break
    if bracket is not None:
        sol = root_scalar(lambda theta: theta_equation(theta), bracket=bracket, method='brentq')
        theta_val = sol.root if sol.converged else None
        if theta_val is None:
            raise ValueError("未解出合适打击轨迹，请检查输入参数")
    else:
        raise ValueError("未找到合适的theta区间")
    tau_val = params.R_m * np.cos(params.Psi) / params.v0 / np.cos(theta_val)
    a0_val = params.a0_max * np.sin(theta_val)
    xT = a0_val * (1 - np.cos(params.wB * tau_val))
    yT = a0_val * np.sin(params.wB * tau_val)
    phi_val = np.arctan2(-xT, yT)
    return theta_val, tau_val, phi_val
