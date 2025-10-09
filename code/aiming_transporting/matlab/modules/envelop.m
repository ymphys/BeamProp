function [sol, R_T] = sol_envelope_ode_2nd(params)
% 求解束流包络二阶ODE

% 初始条件
y0 = [params.R0_m; 0];
z_span = [0, params.z_max];
z_eval = linspace(z_span(1), z_span(2), 200);

% 定义ODE函数
ode_fun = @(z,y) [
    y(2); 
    -2 * x_factor(params.gamma, y(1), params.n_i_m3, params.Ib) * params.Ib / (params.I_A * y(1)) + ...
    epsilon_rms(params.Ib, params.gamma, params.em_type)^2 / y(1)^3
];

% 求解ODE
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
sol = ode45(ode_fun, z_span, y0, options);

% 获取最终半径
R_T = deval(sol, z_span(2));
R_T = R_T(1);
end

function eps = epsilon_rms(Ib, gamma, em_type)
% 计算RMS发射度(m·rad)
switch em_type
    case 'Ti'
        eps = 0.01 * 0.9 / gamma;
    case 'C'
        eps = 0.01 * 0.04 / gamma;
    case 'LP'
        eps = 0.01 * 0.005 * sqrt(Ib) / gamma;
    otherwise
        error('em_type must be ''Ti'', ''C'', or ''LP''');
end
end

function fe = charge_neutralization_factor(b, R, n_i, Ib)
% 计算电荷中和因子
fe = 1.6e-19 * b * 3e8 * pi * R^2 * n_i / Ib;
end

function xfac = x_factor(gamma, R, n_i, Ib)
% 计算x因子
fe = charge_neutralization_factor(gamma, R, n_i, Ib);
xfac = (fe * gamma^2 - 1) / (gamma^2 - 1);
end

function req = plasma_equilibrium_radius(params)
% 计算等离子体平衡半径
eps = epsilon_rms(params.Ib, params.gamma, params.em_type);
denom = 2 * 1.6e-19 * params.n_i_m3 * pi * params.b * 3e8 * params.gamma^2;

if denom == 0
    error('Denominator in plasma_equilibrium_radius is zero.');
end

sqrt_term = sqrt(params.Ib^2 + denom * (params.gamma^2 - 1) * params.I_A * eps^2);
req = sqrt((params.Ib + sqrt_term) / denom);

if isempty(req) || req <= 0 || ~isreal(req)
    error('当前输入参数下不存在束流平衡半径！');
end
end