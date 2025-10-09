% aiming_angle_vs_range.m
clear; clc;

% 常量定义
c = 2.998e8;      % 光速 [m/s]
qe = 1.609e-19;   % 元电荷 [C]
me = 9.11e-31;    % 电子质量 [kg]
gamma = 70;       % 洛伦兹因子
beta = sqrt(1 - 1/gamma^2);
B = 2.5e-5;       % 磁场 [T]

% 变量计算
v0 = c * beta;
wB = qe * B / (gamma * me);
a0_max = v0 / wB;

% 参数组
R1 = 1e3:1e3:15e3; Psi1 = 0.1;
R2 = 1e3:1e3:11e3; Psi2 = 0.5;
R3 = 1e3:1e3:9e3;  Psi3 = 1.0;

theta_list1 = nan(size(R1));
theta_list2 = nan(size(R2));
theta_list3 = nan(size(R3));

% 定义目标函数
a0 = @(theta) a0_max * sin(theta);
theta_equation = @(theta, R, Psi) R .* sin(Psi) - a0(theta) .* sqrt(2 * (1 - cos(wB .* R .* cos(Psi) ./ (v0 .* cos(theta)))));

% 求解函数
find_bracket = @(f, R, Psi) ...
    find_bracket_helper(f, R, Psi, 0.01, pi/2, 1000);

for i = 1:length(R1)
    bracket = find_bracket(theta_equation, R1(i), Psi1);
    if ~isempty(bracket)
        theta_list1(i) = rad2deg(fzero(@(theta) theta_equation(theta, R1(i), Psi1), bracket));
    end
end

for i = 1:length(R2)
    bracket = find_bracket(theta_equation, R2(i), Psi2);
    if ~isempty(bracket)
        theta_list2(i) = rad2deg(fzero(@(theta) theta_equation(theta, R2(i), Psi2), bracket));
    end
end

for i = 1:length(R3)
    bracket = find_bracket(theta_equation, R3(i), Psi3);
    if ~isempty(bracket)
        theta_list3(i) = rad2deg(fzero(@(theta) theta_equation(theta, R3(i), Psi3), bracket));
    end
end

% 绘图
figure;
hold on;
scatter(R1/1e3, theta_list1, 'b', 'filled');
plot(R1/1e3, theta_list1, 'b');
scatter(R2/1e3, theta_list2, 'g', 'filled');
plot(R2/1e3, theta_list2, 'g');
scatter(R3/1e3, theta_list3, 'r', 'filled');
plot(R3/1e3, theta_list3, 'r');
xlabel('R (km)');
ylabel('\theta (deg)');
title('Aiming angle vs Range');
legend('\Psi=0.1', '\Psi=0.5', '\Psi=1.0');
grid on;
hold off;

% --- 辅助函数 ---
function bracket = find_bracket_helper(func, R, Psi, theta_min, theta_max, num)
    theta_test = linspace(theta_min, theta_max, num);
    f_test = arrayfun(@(theta) func(theta, R, Psi), theta_test);
    idx = find(f_test(1:end-1).*f_test(2:end) < 0, 1);
    if ~isempty(idx)
        bracket = [theta_test(idx), theta_test(idx+1)];
    else
        bracket = [];
    end
end