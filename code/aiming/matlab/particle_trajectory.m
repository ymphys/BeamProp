% 常量定义
c = 2.998e8;      % 光速 [m/s]
qe = 1.609e-19;   % 元电荷 [C]
me = 9.11e-31;    % 电子质量 [kg]

% 变量输入
gamma = 35/0.511;           % 洛伦兹因子
beta = sqrt(1 - 1/gamma^2); % 相对论速度
B = 2.5e-5;                 % 磁场 [T]
theta = deg2rad(30);        % 速度与磁场夹角

% 变量计算
v0 = c * beta;
wB = qe * B / (gamma * me);
a0 = v0 * sin(theta) / wB;

% 时间参数
t = linspace(0, 2*pi/wB*5, 1000); % 画5圈

% 轨迹方程
x = a0 * (1 - cos(wB * t));
y = a0 * sin(wB * t);
z = v0 * cos(theta) * t;

% 绘图
figure('Position',[100 100 800 600]);
plot3(x, y, z, 'LineWidth', 2);
hold on;
xlim([0, 2*a0]);
ylim([-a0, a0]);
zlim([0, max(z)]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title('三维电子回旋轨迹图','FontSize',14);

% B场箭头（+z方向）
arrow_len = 0.5 * max(z);
quiver3(0, 0, 0, 0, 0, arrow_len, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(0, 0, arrow_len*1.1, '\bf\itB', 'Color', 'b', 'FontSize', 14);

% 初速度箭头（y-z平面，夹角30度）
v_arrow_len = 0.25 * max(z);
vy = v_arrow_len * sin(theta);
vz = v_arrow_len * cos(theta);
% 箭头
quiver3(0, 0, 0, 0, vy, vz, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% 标签放在箭头终点稍远处
text(0, vy*0.2, vz*0.2, '\bf\itv_0', 'Color', 'r', 'FontSize', 14);

grid on;
view(45, 25);
hold off;