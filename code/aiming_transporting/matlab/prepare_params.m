function params = prepare_params(raw_params)
% 准备仿真参数结构体
params = struct();

% 基本参数
params.E = raw_params.E;  % 电子束能量 [MeV]
params.B = raw_params.B;  % 磁场 [g]
params.R = raw_params.R;  % 目标距离 [km] 
params.Psi = raw_params.Psi; % 目标角度 [rad]
params.n_i = raw_params.n_i; % 背景离子密度 [cm^-3]
params.Ib = raw_params.Ib;   % 电子束电流 [A]
params.em_type = raw_params.em_type; % 发射度类型
params.R0 = raw_params.R0;   % 初始束流半径 [cm]

% 计算衍生参数
params.gamma = params.E / 0.511;
params.B_T = params.B * 1e-4;  % 转换为特斯拉
params.R_m = params.R * 1e3;    % 转换为米
params.n_i_m3 = params.n_i * 1e6; % 转换为m^-3
params.R0_m = params.R0 * 1e-2;  % 转换为米

% 运动学参数
params.b = (1 - 1/params.gamma^2)^0.5;
params.v0 = 3e8 * params.b;  % 光速c=3e8 m/s
params.wB = 1.6e-19 * params.B_T / (params.gamma * 9.1e-31);
params.a0_max = params.v0 / params.wB;
params.I_A = 4 * pi * 8.85e-12 * 9.1e-31 * (3e8)^3 * params.b * params.gamma / 1.6e-19;
end