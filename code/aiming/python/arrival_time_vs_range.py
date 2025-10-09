import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import os

# 获取当前脚本文件的目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 常量定义
c = 2.998e8  # 光速 [m/s]
qe = 1.609e-19  # 元电荷 [C]
me = 9.11e-31  # 电子质量 [kg]
# 变量输入
gamma = 70 # 洛伦兹因子，无量纲
beta = np.sqrt(1 - 1/gamma**2) # 相对论速度，无量纲
B = 2.5e-5 # 磁场横向分量 [T]

# 变量计算
v0 = c * beta  # 电子束速度，单位 [m/s]
wB = qe * B / (gamma * me)  # 回旋频率
# Rc = 1.7e-3*beta*gamma/B # 垂直情况下回旋半径 [m]
a0_max = v0 / wB  # 最大回转半径 [m]
# 求解发射角theta   
def a0(theta):
    return a0_max * np.sin(theta)

def theta_equation(theta, R, Psi):
    a0_val = a0(theta)
    return R * np.sin(Psi) - a0_val * np.sqrt(2 * (1 - np.cos(wB * R * np.cos(Psi) / (v0 * np.cos(theta)))))

def tau(R, theta,Psi):
    return R * np.cos(Psi) / v0 / np.cos(theta)

# 参数组
R1 = np.arange(1e3, 11.1e3, 1e3)
Psi1 = 0.5

tau_list = []

def find_bracket(func, R, Psi, theta_min=0.01, theta_max=np.pi/2, num=1000):
    theta_test = np.linspace(theta_min, theta_max, num)
    f_test = [func(theta, R, Psi) for theta in theta_test]
    for i in range(len(f_test)-1):
        if f_test[i] * f_test[i+1] < 0:
            return [theta_test[i], theta_test[i+1]]
    return None

for R_val in R1:
    bracket = find_bracket(theta_equation, R_val, Psi1)
    if bracket is not None:
        sol = root_scalar(theta_equation, args=(R_val, Psi1), bracket=bracket, method='brentq')
        theta_val = sol.root if sol.converged else None
        tau_list.append(tau(R_val,theta_val,Psi1)*1e6) # convert to microsecond
    else:
        tau_list.append(np.nan)

# 线性拟合 tau_list vs R1
coeff = np.polyfit(R1, tau_list, 1)
fit_line = np.polyval(coeff, R1)
print("拟合参数为",coeff)

# 绘图
plt.figure(figsize=(8,6))
plt.scatter(R1/1e3, tau_list, color='r', label=r'$\Psi=0.5$')
plt.plot(R1/1e3, tau_list, color='r')
plt.plot(R1/1e3, fit_line, 'k--', label='Linear fit')
plt.xlabel('R (km)')
plt.ylabel(r'$\tau$ (ms)')
plt.title('Arrival time vs Range')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'arrival_time_vs_range.png'), dpi=300)
plt.show()