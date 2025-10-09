function write_to_csv(params, theta, phi, tau_us, r_eq_cm, R_T_cm, fluence_rate, filename)
% 保存仿真结果到CSV文件

% 创建结果表格
result_table = table(...
    params.E, params.B, params.R, params.Psi, ...
    params.n_i, params.Ib, params.em_type, params.R0, ...
    theta, phi, tau_us, r_eq_cm, R_T_cm, fluence_rate, ...
    'VariableNames', {...
    'E_MeV', 'B_g', 'R_km', 'Psi_rad', ...
    'n_i_cm3', 'Ib_A', 'em_type', 'R0_cm', ...
    'theta_rad', 'phi_rad', 'tau_us', 'r_eq_cm', 'R_T_cm', 'fluence_rate'});

% 写入文件
writetable(result_table, filename);

% 添加单位说明
fid = fopen(filename, 'a');
fprintf(fid, '\n\nUnits:\n');
fprintf(fid, 'E: MeV\nB: gauss\nR: km\nPsi: rad\n');
fprintf(fid, 'n_i: cm^-3\nIb: A\nR0: cm\n');
fprintf(fid, 'theta: rad\nphi: rad\ntau: μs\n');
fprintf(fid, 'r_eq: cm\nR_T: cm\nfluence_rate: cm^-2s^-1');
fclose(fid);
end