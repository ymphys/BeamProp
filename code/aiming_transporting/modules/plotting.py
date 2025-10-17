from modules.beam_params import BeamParams
import numpy as np
import matplotlib.pyplot as plt

def plot_trajectory(theta_val, tau_val, params: BeamParams):
    a0_val = params.a0_max * np.sin(theta_val)
    t = np.linspace(0, tau_val, 1000)
    x = a0_val * (1 - np.cos(params.wB * t))
    y = a0_val * np.sin(params.wB * t)
    z = params.v0 * np.cos(theta_val) * t
    xT = a0_val * (1 - np.cos(params.wB * tau_val))
    yT = a0_val * np.sin(params.wB * tau_val)
    zT = params.v0 * np.cos(theta_val) * tau_val
    phi = np.arctan2(-xT, yT)
    xp = np.cos(phi) * x + np.sin(phi) * y
    yp = -np.sin(phi) * x + np.cos(phi) * y
    xTp = np.cos(phi) * xT + np.sin(phi) * yT
    yTp = -np.sin(phi) * xT + np.cos(phi) * yT
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xp, yp, z)
    ax.plot([0, xTp], [0, yTp], [0, zT])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.set_title('电子束打击靶目标轨迹示意图')
    info = f'目标距离: {params.R_m*1e-3:.1f} km, 目标角度: {params.Psi:.2f} rad, 打击用时: {tau_val*1e6:.2f} $\\mu$ s'
    plt.annotate(info, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, color='black', verticalalignment='top', bbox=dict(boxstyle='round', fc='w', alpha=0.7))
    plt.tight_layout()
    plt.savefig(f'aiming_R{(params.R_m)}_Psi{(params.Psi)}.jpg', dpi=300, bbox_inches='tight')
    plt.show()

def plot_envelope(sol, params: BeamParams):
    plt.plot(sol.t * 1e-3, sol.y[0] * 100)  # x轴: m->km, y轴: m->cm
    plt.xlabel('传输距离(km)')
    plt.ylabel('束流包络半径(cm)')
    plt.title('束流包络半径随传输距离的变化')
    info = f'目标距离: {params.R} km, 电子束电流: {params.Ib} A, 背景离子密度: {params.n_i:.2e} cm^-3'
    plt.annotate(info, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, color='black', verticalalignment='top', bbox=dict(boxstyle='round', fc='w', alpha=0.7))
    plt.tight_layout()
    plt.savefig(f'envelop_R{(params.R)}_R_0{(params.R0)}.jpg', dpi=300, bbox_inches='tight')
    plt.grid()
    plt.show()