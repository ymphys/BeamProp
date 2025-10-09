function [theta_val, tau_val, phi_val] = solve_theta(params)
% 求解电子束发射角度和飞行参数

% 定义辅助函数
a0 = @(theta) params.a0_max * sin(theta);

theta_equation = @(theta) params.R_m * sin(params.Psi) - ...
    a0(theta) * sqrt(2 * (1 - cos(params.wB * params.R_m * cos(params.Psi) / ...
    (params.v0 * cos(theta))));

% 寻找解区间
theta_min = 0.01;
theta_max = pi/2;
theta_test = linspace(theta_min, theta_max, 1000);
f_test = arrayfun(theta_equation, theta_test);

bracket = [];
for i = 1:length(f_test)-1
    if f_test(i) * f_test(i+1) < 0
        bracket = [theta_test(i), theta_test(i+1)];
        break;
    end
end

if ~isempty(bracket)
    options = optimset('Display','off');
    [theta_val, ~, exitflag] = fzero(theta_equation, bracket, options);
    
    if exitflag ~= 1
        error('未解出合适打击轨迹，请检查输入参数');
    end
else
    error('未找到合适的theta区间');
end

% 计算其他参数
tau_val = params.R_m * cos(params.Psi) / params.v0 / cos(theta_val);
a0_val = a0(theta_val);
xT = a0_val * (1 - cos(params.wB * tau_val));
yT = a0_val * sin(params.wB * tau_val);
phi_val = atan2(-xT, yT);
end

function check_range(params)
% 检查目标是否在射程内
max_range = 2 * params.a0_max / sin(params.Psi);
if params.R_m > max_range
    error('目标距离超出射程！当前目标角度下射程需小于: %.1f公里', max_range*1e-3);
end
end