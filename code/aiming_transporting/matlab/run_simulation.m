function run_simulation(params)
% 执行单次电子束仿真

% 调用各计算模块
aiming.check_range(params);
[theta_val, tau_val, phi_val] = aiming.solve_theta(params);
r_eq = envelop.plasma_equilibrium_radius(params);
params.z_max = tau_val * params.v0;

% 求解束流包络ODE
[sol, R_T] = envelop.sol_envelope_ode_2nd(params);

% 计算注量率
fluence_rate = params.Ib / (1.6e-19 * pi * R_T^2);

% 绘制结果
plotting.plot_trajectory(theta_val, tau_val, params);
plotting.plot_envelope(sol, params);

% 保存结果
saving_data.write_to_csv(params, theta_val, phi_val, tau_val*1e6, ...
    r_eq*1e2, R_T*1e-2, fluence_rate*1e-4, 'results.csv');
end