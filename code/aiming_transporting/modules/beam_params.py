from dataclasses import dataclass

@dataclass
class BeamParams:
    E: float           # 电子束能量 (MeV)
    B: float           # 磁场 (g)
    R: float           # 目标距离 (km)
    Psi: float         # 目标角度 (rad)
    n_i: float         # 背景离子密度 (cm^-3)
    Ib: float          # 电子束电流 (A)
    em_type: str       # 发射度类型
    R0: float          # 初始束流包络半径 (cm)
    gamma: float = None    # 相对论因子
    B_T: float = None      # 磁场 (T)
    R_m: float = None      # 目标距离 (m)
    n_i_m3: float = None   # 背景离子密度 (m^-3)
    R0_m: float = None     # 初始包络半径 (m)
    b: float = None        # 相对论速度因子
    v0: float = None       # 束流速度 (m/s)
    wB: float = None       # 回旋频率
    a0_max: float = None   # 最大回转半径 (m)
    I_A: float = None      # Alfven current (A)
    z_max: float = None    # 实际传输距离 (m)
