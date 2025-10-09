from modules import ui
from modules.beam_params import BeamParams
from modules import constants, aiming, envelop, plotting, saving_data
import pandas as pd


def prepare_params(raw_params):
    params = BeamParams(
        E=raw_params['E'],
        B=raw_params['B'],
        R=raw_params['R'],
        Psi=raw_params['Psi'],
        n_i=raw_params['n_i'],
        Ib=raw_params['Ib'],
        em_type=raw_params['em_type'],
        R0=raw_params['R0']
    )
    params.gamma = params.E / 0.511
    params.B_T = params.B * 1e-4
    params.R_m = params.R * 1e3
    params.n_i_m3 = params.n_i * 1e6
    params.R0_m = params.R0 * 1e-2
    params.b = (1 - 1/params.gamma**2) ** 0.5
    params.v0 = constants.C * params.b
    params.wB = constants.Q_E * params.B_T / (params.gamma * constants.M_E)
    params.a0_max = params.v0 / params.wB
    params.I_A = 4 * 3.141592653589793 * constants.EPSILON_0 * constants.M_E * constants.C ** 3 * params.b * params.gamma / constants.Q_E
    return params


def run_simulation(params: BeamParams):
    aiming.check_range(params)
    theta_val, tau_val, phi_val = aiming.solve_theta(params)
    r_eq = envelop.plasma_equilibrium_radius(params)
    params.z_max = tau_val * params.v0
    sol, R_T = envelop.sol_envelope_ode_2nd(params)
    fluence_rate = params.Ib / (constants.Q_E * 3.141592653589793 * R_T**2)
    plotting.plot_trajectory(theta_val, tau_val, params)
    plotting.plot_envelope(sol, params)
    saving_data.write_to_csv(
        params, theta_val, phi_val, tau_val*1e6, r_eq*1e2, R_T*1e2, fluence_rate*1e-4, 'results.csv')


def run_batch_simulation(csv_path):
    df = pd.read_csv(csv_path)
    for idx, row in df.iterrows():
        try:
            params = prepare_params(row)
            run_simulation(params)
        except Exception as e:
            print(f"\n第{idx+1}组参数仿真失败:")
            print(f"参数: {dict(row)}")
            print(f"错误详情: {str(e)}")
            import traceback
            traceback.print_exc()
            print("-"*50)


if __name__ == '__main__':
    user_inputs = ui.get_user_inputs()
    if user_inputs is None:
        print("未获取到输入参数，程序终止。")
    elif user_inputs.get('mode') == 'batch':
        run_batch_simulation(user_inputs['csv_path'])
    else:
        params = prepare_params(user_inputs)
        run_simulation(params)