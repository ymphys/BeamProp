import numpy as np
# me = 0.511
# mp = 938
# mH = 938.5
c = 2.998e8
qe = 1.609e-19  # 元电荷 [C]
me = 9.11e-31  # 电子质量 [kg]
mp = 938/0.511*me
# 变量输入
gamma = 70 # 洛伦兹因子，无量纲
beta = np.sqrt(1 - 1/gamma**2) # 相对论速度，无量纲
B = 1.0e-7 # 磁场横向分量 [T]

# 变量计算
v0 = c * beta  # 电子束速度，单位 [m/s]
wB = qe * B / (gamma * me)  # 回旋频率
# Rc = 1.7e-3*beta*gamma/B # 垂直情况下回旋半径 [m]
a0_max = v0 / wB  # 最大回转半径 [m]

print(2*a0_max)

# def calc_v(m,Ek):
#     """
#     calculate the velocity of a particle with mass m (in MeV/c^2) and kinetic energy Ek (in MeV)
#     """
#     gamma = Ek/m + 1
#     beta = (1 - 1/gamma**2)**0.5
#     return gamma,beta

# import matplotlib.pyplot as plt
# Ek_lst = np.linspace(50,500,450)
# beta_lst = []
# t_lst = []
# for Ek in Ek_lst:
#     _,beta = calc_v(mp,Ek)
#     t = 1e3/(beta*c)*1e6 # in us
#     beta_lst.append(beta)
#     t_lst.append(t)
# plt.plot(Ek_lst,beta_lst)
# plt.xlabel('Kinetic Energy (MeV)')
# plt.ylabel('Beta (v/c)')
# plt.title('Proton Velocity vs Kinetic Energy')
# plt.legend()
# plt.savefig("proton_v_Ek.png")
# plt.close()

# plt.plot(Ek_lst,t_lst)
# plt.xlabel('Kinetic Energy (MeV)')
# plt.ylabel('prop 1km time used (us)')
# plt.title('Proton Prop. time in 1 km vs Kinetic Energy')
# plt.legend()
# plt.savefig("proton_t_Ek.png")
# plt.close()