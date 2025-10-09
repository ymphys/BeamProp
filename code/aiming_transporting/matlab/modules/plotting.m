function plot_trajectory(theta_val, tau_val, params)
% 绘制电子束轨迹
figure;
hold on;

% 计算轨迹点
t = linspace(0, tau_val, 100);
x = params.a0_max * sin(theta_val) * (1 - cos(params.wB * t));
y = params.a0_max * sin(theta_val) * sin(params.wB * t);
z = params.v0 * cos(theta_val) * t;

% 绘制3D轨迹
plot3(z*1e-3, x*1e-3, y*1e-3, 'LineWidth', 2);
plot3([0 params.R_m*cos(params.Psi)*1e-3], [0 0], [0 params.R_m*sin(params.Psi)*1e-3], 'r--');

% 标注目标点
plot3(params.R_m*cos(params.Psi)*1e-3, 0, params.R_m*sin(params.Psi)*1e-3, 'ro', 'MarkerSize', 10);

% 设置图形属性
xlabel('距离 (km)');
ylabel('横向位移 (km)');
zlabel('高度 (km)');
title(sprintf('电子束轨迹 (θ=%.2f rad)', theta_val));
grid on;
axis equal;
view(30, 30);
end

function plot_envelope(sol, params)
% 绘制束流包络
figure;
hold on;

% 获取解数据
z = sol.x;
R = sol.y(1,:);

% 绘制包络
plot(z*1e-3, R*1e2, 'b-', 'LineWidth', 2);
plot(z*1e-3, -R*1e2, 'b-', 'LineWidth', 2);

% 设置图形属性
xlabel('传播距离 (km)');
ylabel('束流半径 (cm)');
title('束流包络演化');
grid on;
end