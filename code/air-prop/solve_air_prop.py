import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, alpha, physical_constants, Avogadro, pi

# ---------- 常数 ----------
m_e = physical_constants["electron mass energy equivalent in MeV"][0]  # MeV
r_e = physical_constants["classical electron radius"][0]              # m
I_A = 17045.0   # Alfven current [A]

rho = 1.225e-3  # g/cm^3 = 1.225 kg/m^3
A = 14.0        # 平均原子量 (空气 ~ N2)
Z = 7.0
I_Z = 85.7e-6   # MeV
NA = Avogadro

X0 = 300.5      # m, 辐射长度
Ez = 0.0        # 外场 [MeV/m]

# ---------- 束流初始条件 ----------
W0 = 30.0       # MeV
a0 = 1e-3       # m
ap0 = 0.0       # rad
eps_tr0 = 1e-6  # m·rad
I_E = 1       # A

# ---------- 动量函数 ----------
def beta_gamma(W):
    gamma = 1.0 + W/m_e
    beta = np.sqrt(1.0 - 1.0/gamma**2)
    return beta, gamma

def pc_from_W(W):
    beta, gamma = beta_gamma(W)
    return beta * gamma * m_e   # MeV/c

# ---------- 损失项 ----------
def S_of_beta(beta, gamma):
    term = np.log(2*beta**2*gamma**2*m_e/I_Z) - beta**2
    pref = 4*pi*r_e**2 * NA * rho/A * Z * m_e
    return pref * (1.0/beta**2) * term  # MeV/m

def dW_dz(W, z):
    if W <= 0:
        return 0.0
    beta, gamma = beta_gamma(W)
    S = S_of_beta(beta, gamma)
    return -S - W/X0 + Ez

# ---------- 积分参数 ----------
zmax = 50.0      # m
dz   = 0.01       # m
Nz   = int(zmax/dz)

# 数组
zgrid = np.zeros(Nz+1)
Wgrid = np.zeros(Nz+1)
agrid = np.zeros(Nz+1)
apgrid= np.zeros(Nz+1)

# 初始值
zgrid[0] = 0.0
Wgrid[0] = W0
agrid[0] = a0
apgrid[0]= ap0

# 累积积分 ∫ a^2 p^2 dψ^2
Iint_grid = np.zeros(Nz+1)
p0 = pc_from_W(W0)

# ---------- 主循环 ----------
for i in range(Nz):
    z = zgrid[i]
    W = Wgrid[i]
    a = agrid[i]
    ap= apgrid[i]

    # 能量保护
    if W <= 1e-6:
        zgrid = zgrid[:i+1]
        Wgrid = Wgrid[:i+1]
        agrid = agrid[:i+1]
        apgrid= apgrid[:i+1]
        Iint_grid = Iint_grid[:i+1]
        break

    # ---- 更新能量 (RK4) ----
    k1 = dW_dz(W, z)
    k2 = dW_dz(W + 0.5*dz*k1, z + 0.5*dz)
    k3 = dW_dz(W + 0.5*dz*k2, z + 0.5*dz)
    k4 = dW_dz(W + dz*k3,     z + dz)
    W_new = W + (dz/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    # ---- β, γ, p ----
    beta, gamma = beta_gamma(W)
    pc = pc_from_W(W)

    # ---- dψ² ----
    Es = np.sqrt(4*pi/alpha) * m_e  # ≈ 21.2 MeV
    dpsi2_dz = (Es/(beta*pc))**2 / X0

    # ---- 积分 Iint ----
    Iint_new = Iint_grid[i] + a**2 * pc**2 * dpsi2_dz * dz

    # ---- 包络方程 RHS ----
    def envelope_rhs(a, ap, W, Iint):
        beta, gamma = beta_gamma(W)
        pc = pc_from_W(W)
        Wp = dW_dz(W, z)
        term1 = (I_E/I_A) / max(a,1e-9)
        term2 = (ap*Wp) / max(beta**2*W, 1e-12)
        RHS   = (p0**2*eps_tr0**2 + Iint) / (a**3 * pc**2)
        a_ddot = RHS - term1 - term2
        return ap, a_ddot

    # ---- RK4 for (a, ap) ----
    k1a, k1ap = envelope_rhs(a, ap, W, Iint_new)
    k2a, k2ap = envelope_rhs(a + 0.5*dz*k1a, ap + 0.5*dz*k1ap, W, Iint_new)
    k3a, k3ap = envelope_rhs(a + 0.5*dz*k2a, ap + 0.5*dz*k2ap, W, Iint_new)
    k4a, k4ap = envelope_rhs(a + dz*k3a,     ap + dz*k3ap,     W, Iint_new)

    a_new  = a  + (dz/6.0)*(k1a + 2*k2a + 2*k3a + k4a)
    ap_new = ap + (dz/6.0)*(k1ap+ 2*k2ap+ 2*k3ap+ k4ap)

    # ---- 存储 ----
    zgrid[i+1] = z + dz
    Wgrid[i+1] = max(W_new, 0.0)
    agrid[i+1] = max(a_new, 1e-9)
    apgrid[i+1]= ap_new
    Iint_grid[i+1] = Iint_new

# ---------- 绘图 ----------
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(zgrid, Wgrid)
plt.xlabel("z [m]")
plt.ylabel("W [MeV]")
plt.title("Energy evolution")

plt.subplot(1,2,2)
plt.plot(zgrid, agrid*1e3)
plt.xlabel("z [m]")
plt.ylabel("a [mm]")
plt.title("Envelope evolution (RK4)")

plt.tight_layout()
plt.savefig("air_propagation.png", dpi=300)
