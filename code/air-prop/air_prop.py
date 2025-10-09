import numpy as np
from math import erf

# ======================== 常量 ========================
m_e = 0.511  # MeV, electron rest mass
rho_STP = 1.225e-3  # g/cm^3, air density at STP
X0_air = 36.66      # g/cm^2, radiation length of air
K = 0.307           # MeV cm^2/g, Bethe constant
e_charge = 1.602e-19  # C

# ======================== 简化空气密度模型 ========================
def air_density(h_km):
    rho0 = rho_STP
    H = 8.5  # km, scale height
    return rho0 * np.exp(-h_km / H)

# ======================== Bethe-Bloch 能量损失 ========================
def dE_dx(E, rho):
    gamma = 1 + E / m_e
    beta2 = 1 - 1 / gamma**2
    Wmax = 2 * m_e * beta2 * gamma**2
    I = 85.7e-6  # MeV
    return K * rho * (1 / beta2) * (0.5 * np.log(2 * m_e * beta2 * gamma**2 * Wmax / I**2) - beta2)

# ======================== Highland 散射角 ========================
def theta0(E, rho, L):
    gamma = 1 + E / m_e
    beta = np.sqrt(1 - 1 / gamma**2)
    p = np.sqrt(E**2 + 2 * E * m_e)  # MeV/c
    x = rho * L / X0_air
    return (13.6 / (beta * p)) * np.sqrt(x) * (1 + 0.038 * np.log(x + 1e-12))

# ======================== 传播计算 ========================
def propagate(E0, r0, L, rho, nsteps=1000):
    dz = L / nsteps
    E = E0
    r = r0
    for _ in range(nsteps):
        E -= dE_dx(E, rho) * dz
        if E <= 0:
            E = 0
            break
        th = theta0(E, rho, dz)
        r = np.sqrt(r**2 + (dz * th)**2)
    return E, r

# ======================== 轫致辐射能散导致的电流衰减 ========================
def bremsstrahlung_current_loss(E0, L, rho, I0):
    """返回剩余电流和平均能量"""
    x = rho * L  # g/cm^2
    E_mean = E0 * np.exp(-x / X0_air)
    sigma_E = E_mean * np.sqrt(x / X0_air)
    # 保留比例
    keep_fraction = 0.5 * (1 + erf((E_mean - 0.5 * E0) / (np.sqrt(2) * sigma_E)))
    I_final = I0 * keep_fraction
    return I_final, E_mean

# ======================== 用户输入 ========================
h_km = float(input("请输入高度 (km)："))
L_km = float(input("请输入传输距离 (km)："))
I0 = float(input("请输入初始电流 I0 (A)："))

# 初始条件
E0 = 35.0   # MeV
r0 = 0.15   # cm (1.5 mm)
rho = air_density(h_km)
L = L_km * 1e5  # cm

# 计算传播结果
E_final, r_final = propagate(E0, r0, L, rho)
I_final, E_mean = bremsstrahlung_current_loss(E0, L, rho, I0)

# 注量率 (cm^-2s^-1)
P = I_final /(np.pi*r_final**2) / 1.602e-19

print(f"高度 {h_km:.1f} km, 传输 {L_km:.2f} km 后：")
print(f"  最终束流能量 E = {E_final:.2f} MeV")
print(f"  束流 RMS 半径 r = {r_final:.2f} cm")
print(f"  剩余电流 I = {I_final:.2f} A (初始 I0 = {I0:.2f} A)")
print(f"  最终到靶注量率 P = {P:.3e} cm^-2s^-1")
