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

# 参数组
R1 = np.arange(1e3, 15.1e3, 1e3)
Psi1 = 0.1
R2 = np.arange(1e3, 11.1e3, 1e3)
Psi2 = 0.5
R3 = np.arange(1e3, 9.1e3, 1e3)
Psi3 = 1.0

theta_list1 = []
theta_list2 = []
theta_list3 = []

def find_bracket(func, R, Psi, theta_min=0.01, theta_max=np.pi/2, num=1000):
    theta_test = np.linspace(theta_min, theta_max, num)
    f_test = [func(theta, R, Psi) for theta in theta_test]
    for i in range(len(f_test)-1):
        if f_test[i] * f_test[i+1] < 0:
            return [theta_test[i], theta_test[i+1]]
    return None

for R in R1:
    bracket = find_bracket(theta_equation, R, Psi1)
    if bracket is not None:
        sol = root_scalar(theta_equation, args=(R, Psi1), bracket=bracket, method='brentq')
        theta_list1.append(np.degrees(sol.root) if sol.converged else np.nan)
    else:
        theta_list1.append(np.nan)

for R in R2:
    bracket = find_bracket(theta_equation, R, Psi2)
    if bracket is not None:
        sol = root_scalar(theta_equation, args=(R, Psi2), bracket=bracket, method='brentq')
        theta_list2.append(np.degrees(sol.root) if sol.converged else np.nan)
    else:
        theta_list2.append(np.nan)

for R in R3:
    bracket = find_bracket(theta_equation, R, Psi3)
    if bracket is not None:
        sol = root_scalar(theta_equation, args=(R, Psi3), bracket=bracket, method='brentq')
        theta_list3.append(np.degrees(sol.root) if sol.converged else np.nan)
    else:
        theta_list3.append(np.nan)
# 绘图
plt.figure(figsize=(8,6))
plt.scatter(R1/1e3, theta_list1, color='b', label=r'$\Psi=0.1$')
plt.plot(R1/1e3, theta_list1, color='b')
plt.scatter(R2/1e3, theta_list2, color='g', label=r'$\Psi=0.5$')
plt.plot(R2/1e3, theta_list2, color='g')
plt.scatter(R3/1e3, theta_list3, color='r', label=r'$\Psi=1.0$')
plt.plot(R3/1e3, theta_list3, color='r')
plt.xlabel('R (km)')
plt.ylabel(r'$\theta$ (deg)')
plt.title('Aiming angle vs Range')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'aiming_angle_vs_range.png'), dpi=300)
plt.show()