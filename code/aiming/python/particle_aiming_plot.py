import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import root_scalar
from matplotlib.widgets import TextBox, Button
import os
import csv

# 常量定义
c = 2.998e8  # 光速 [m/s]
qe = 1.609e-19  # 元电荷 [C]
me = 9.11e-31  # 电子质量 [kg]

# 创建matplotlib输入界面
matplotlib.rcParams['font.sans-serif'] = ['Heiti TC', 'STHeiti', 'SimHei', 'Arial Unicode MS']
matplotlib.rcParams['axes.unicode_minus'] = False
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.4)
plt.title("电子束瞄准参数输入")
plt.axis('off')

# 创建输入框
plt.subplots_adjust(bottom=0.2)  # 调整底部边距
E_box = TextBox(plt.axes([0.3, 0.7, 0.4, 0.05]), '电子束能量:', initial="35")
B_box = TextBox(plt.axes([0.3, 0.6, 0.4, 0.05]), '磁场 B [g]:', initial="0.25")
R_box = TextBox(plt.axes([0.3, 0.5, 0.4, 0.05]), '目标距离 R [km]:', initial="15")
Psi_box = TextBox(plt.axes([0.3, 0.4, 0.4, 0.05]), '目标角度 Ψ [rad]:', initial="0.5")

# 确认按钮和回车键处理
def submit(event):
    global gamma, B, R_test, Psi_test
    try:
        gamma = float(E_box.text) / 0.511
        B = float(B_box.text) *1e-4
        R_test = float(R_box.text) *1e3
        Psi_test = float(Psi_box.text)
        plt.close()
    except ValueError:
        print("请输入有效的数值参数")

button_ax = plt.axes([0.4, 0.3, 0.2, 0.05])  # 调整按钮位置到输入框下方
button = Button(button_ax, '确认')
button.on_clicked(submit)

# 绑定回车键到submit函数
def on_key(event):
    if event.key == 'enter':
        submit(event)

fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()

# 变量计算
beta = np.sqrt(1 - 1/gamma**2) # 相对论速度，无量纲
v0 = c * beta  # 电子束速度，单位 [m/s]
wB = qe * B / (gamma * me)  # 回旋频率
# Rc = 1.7e-3*beta*gamma/B # 垂直情况下回旋半径 [m]
a0_max = v0 / wB  # 最大回转半径 [m]

# 判断是否超出射程
if R_test > 2*a0_max / np.sin(Psi_test):
    raise ValueError("目标距离超出射程！当前目标角度下射程需小于: {:.1f}公里".format(2*a0_max / np.sin(Psi_test)*1e-3)+"请选择距离或角度更小的目标后重试")

# 求解发射角theta   
def a0(theta):
    return a0_max * np.sin(theta)

def theta_equation(theta, R, Psi):
    a0_val = a0(theta)
    return R * np.sin(Psi) - a0_val * np.sqrt(2 * (1 - np.cos(wB * R * np.cos(Psi) / (v0 * np.cos(theta)))))

def tau(R, theta,Psi):
    return R * np.cos(Psi) / v0 / np.cos(theta)

def find_bracket(func, R, Psi, theta_min=0.01, theta_max=np.pi/2, num=1000):
    theta_test = np.linspace(theta_min, theta_max, num)
    f_test = [func(theta, R, Psi) for theta in theta_test]
    for i in range(len(f_test)-1):
        if f_test[i] * f_test[i+1] < 0:
            return [theta_test[i], theta_test[i+1]]
    return None

bracket = find_bracket(theta_equation, R_test, Psi_test)
if bracket is not None:
    sol = root_scalar(theta_equation, args=(R_test, Psi_test), bracket=bracket, method='brentq')
    theta_val = sol.root if sol.converged else None
    if theta_val is None:
        raise ValueError("未解出合适打击轨迹，请检查输入参数")
        
    # 判断是否超出射程
    if R_test > 2*a0_max / np.sin(Psi_test):
        raise ValueError("目标超出射程！最大射程为: {:.1f}米".format(2*a0_max / np.sin(Psi_test)))
        
    tau_val = tau(R_test,theta_val,Psi_test)
    
# 电子轨迹3D图
a0_val = a0(theta_val)
# 时间参数
t = np.linspace(0, tau_val, 1000)
# 轨迹方程
x = a0_val * (1 - np.cos(wB * t))
y = a0_val * np.sin(wB * t)
z = v0 * np.cos(theta_val) * t
# 靶目标位置
xT = a0_val * (1 - np.cos(wB * tau_val))
yT = a0_val * np.sin(wB * tau_val)
zT = v0 * np.cos(theta_val) * tau_val
# 坐标变换，y-y'夹角
phi = np.arctan2(-xT,yT)
xp = np.cos(phi) * x + np.sin(phi) * y
yp = - np.sin(phi) * x + np.cos(phi) * y
xTp = np.cos(phi) * xT + np.sin(phi) * yT
yTp = - np.sin(phi) * xT + np.cos(phi) * yT

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot(xp, yp, z)
ax.plot([0, xTp], [0, yTp], [0, zT])
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
ax.set_title('三维电子打击靶目标示意图')

plt.tight_layout()
plt.savefig(f'电子束瞄准_R{(R_test)}_Psi{(Psi_test)}.jpg', dpi=300, bbox_inches='tight')
plt.show()

# CSV文件头
header = [
    "电子束能量(MeV)", "磁场(g)", "目标距离(km)",
    "目标角度(rad)", "极角(deg)", "方位角(deg)", "打击时间(μs)"
]

# 准备数据行
data = [
    round(float(E_box.text), 1),
    round(float(B_box.text), 2),
    round(float(R_box.text), 1),
    round(Psi_test, 2),
    round(np.degrees(theta_val), 1) if 'theta_val' in globals() else 'nan',
    round(np.rad2deg(np.pi/2-phi), 1) if 'phi' in globals() else 'nan',
    round(tau_val*1e6, 1) if 'tau_val' in globals() else 'nan'
]

# 写入CSV文件
import sys
from pathlib import Path

# 获取正确的输出目录
if getattr(sys, 'frozen', False):
    # 如果是打包后的EXE
    base_dir = Path(sys.executable).parent
else:
    # 如果是脚本模式
    base_dir = Path(__file__).parent

csv_path = base_dir / 'parameters.csv'
file_exists = csv_path.exists()

with open(csv_path, 'a', newline='', encoding='utf_8_sig') as f:
    writer = csv.writer(f)
    if not file_exists or os.path.getsize(csv_path) == 0:
        writer.writerow(header)
    writer.writerow(data)