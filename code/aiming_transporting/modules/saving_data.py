from modules.beam_params import BeamParams
import numpy as np
import os
import csv

def write_to_csv(params: BeamParams, theta, phi, tau, r_eq, R_T, fluence_rate, file_path='results.csv'):
    header = [
        "电子束能量(MeV)", "磁场(g)", "目标距离(km)",
        "目标角度(rad)", "极角(deg)", "方位角(deg)", "打击时间(μs)","初始束流半径(cm)","平衡束流半径(cm)",
        "背景离子密度(cm^-3)", "电子束电流(A)", "发射度类型",
        "到靶束流包络半径(cm)", "到靶注量率(cm^-2s^-1)"
    ]
    data = [
        params.E, params.B, params.R, params.Psi,
        round(np.degrees(theta), 1) if theta is not None else 'nan',
        round(np.rad2deg(np.pi/2-phi), 1) if phi is not None else 'nan',
        round(tau, 1) if tau is not None else 'nan',
        params.R0, round(r_eq,1) if r_eq is not None else 'nan',
        f'{params.n_i:.1e}' if params.n_i is not None else 'nan',
        params.Ib, params.em_type,
        round(R_T, 1) if R_T is not None else 'nan',
        f'{fluence_rate:.1e}' if fluence_rate is not None else 'nan'
    ]
    with open(file_path, 'a', newline='', encoding='utf_8_sig') as f:
        writer = csv.writer(f)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            writer.writerow(header)
        writer.writerow(data)