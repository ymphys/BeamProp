import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def beam_path(R, psi, gamma, B):
    # 常量定义
    c = 2.998e8  # 光速 [m/s]
    qe = 1.609e-19  # 元电荷 [C]
    me = 9.11e-31  # 电子质量 [kg]
    
    v0 = c * np.sqrt(1 - 1/gamma**2)  # 电子束速度
    wB = qe * B / (gamma * me)  # 回旋频率
    
    # 生成 theta 的测试区间，计算目标函数
    theta_test = np.linspace(0, np.pi/2, 100000)
    fn_test = (R * np.sin(psi))**2 - ((v0 * np.sin(theta_test) / wB)**2) * (2 - 2 * np.cos(wB * R * np.cos(psi) / (v0 * np.cos(theta_test))))
    
    # 找到第一个使目标函数小于0的 theta 区间
    neg_indices = np.where(fn_test < 0)[0]
    if len(neg_indices) == 0:
        print("Error: No valid theta found. Check your input parameters.")
        return
    index_test = neg_indices[0]
    
    r0 = R * np.sin(psi)  # 目标点在 y 方向的投影
    r0max = 2 * v0 / wB   # 最大可达半径
    
    if r0 > r0max:
        print("Error: Target out of range.")
    
    # 绘制目标函数随 theta 的变化曲线
    plt.plot(theta_test, fn_test)
    plt.grid(True)
    plt.show()
    
    # 定义用于求根的函数
    def root_fn(theta):
        return (R * np.sin(psi))**2 - (v0 * np.sin(theta) / wB)**2 * (2 - 2 * np.cos(wB * R * np.cos(psi) / (v0 * np.cos(theta))))
    
    # 用 fsolve 在找到的区间内求解 theta
    theta = fsolve(root_fn, [theta_test[index_test-1], theta_test[index_test]])[0]
    print('theta=',np.rad2deg(theta), 'deg')
    tmax = R * np.cos(psi) / (v0 * np.cos(theta))  # 到达目标点所需时间
    a0 = np.abs(v0 * np.sin(theta) / wB)  # 回旋半径
    
    phi = np.arctan2(1 - np.cos(wB * tmax), np.sin(wB * tmax))  # 轨道相位
    
    t = np.linspace(0, tmax, 1000)  # 时间序列
    x = a0 * (np.cos(wB * t) - 1)   # x 方向轨迹
    y = a0 * np.sin(wB * t)         # y 方向轨迹
    z = v0 * np.cos(theta) * t      # z 方向轨迹
    
    # 坐标旋转，得到实际空间轨迹
    xp = np.cos(phi) * x + np.sin(phi) * y
    yp = np.cos(phi) * y - np.sin(phi) * x
    
    xT = 0
    yT = R * np.sin(psi)
    zT = R * np.cos(psi)
    
    # 3D 绘图，显示粒子轨迹和目标点连线
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xp, yp, z)
    ax.plot([0, xT], [0, yT], [0, zT])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    
    # 计算瞄准角度
    xp = np.sin(theta) * np.cos(phi) * np.cos(psi) - np.cos(theta) * np.sin(psi)
    yp = -np.sin(theta) * np.sin(phi)
    zp = np.sin(theta) * np.cos(phi) * np.sin(psi) + np.cos(theta) * np.cos(psi)
    
    thetaT = np.arctan(xp / zp)  # 极角
    phiT = np.arctan(yp / zp)    # 仰角
    
    print(f"Targeting angles are {thetaT * 180 / np.pi} deg. polar, {phiT * 180 / np.pi} deg. elevation.")

# 示例调用
beam_path(1e3, 0.1, 70, 2.5e-5)