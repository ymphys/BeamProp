function particle_aiming_plot()
% 电子束瞄准参数计算和可视化工具 - MATLAB版本
% 常量定义
c = 2.998e8;  % 光速 [m/s]
qe = 1.609e-19;  % 元电荷 [C]
me = 9.11e-31;  % 电子质量 [kg]

% 创建主GUI界面
mainFig = figure('Name', '电子束瞄准参数输入', 'NumberTitle', 'off', ...
                'Position', [100, 100, 400, 400], 'Resize', 'off');
movegui(mainFig, 'center');

% 创建输入框和标签
uicontrol('Style', 'text', 'Position', [50, 350, 100, 20], ...
         'String', '电子束能量:', 'HorizontalAlignment', 'left');
E_box = uicontrol('Style', 'edit', 'Position', [150, 350, 200, 20], ...
                'String', '35', 'Tag', 'E_box');

uicontrol('Style', 'text', 'Position', [50, 300, 100, 20], ...
         'String', '磁场 B [g]:', 'HorizontalAlignment', 'left');
B_box = uicontrol('Style', 'edit', 'Position', [150, 300, 200, 20], ...
                'String', '0.25', 'Tag', 'B_box');

uicontrol('Style', 'text', 'Position', [50, 250, 100, 20], ...
         'String', '目标距离 R [km]:', 'HorizontalAlignment', 'left');
R_box = uicontrol('Style', 'edit', 'Position', [150, 250, 200, 20], ...
                'String', '15', 'Tag', 'R_box');

uicontrol('Style', 'text', 'Position', [50, 200, 100, 20], ...
         'String', '目标角度 Ψ [rad]:', 'HorizontalAlignment', 'left');
Psi_box = uicontrol('Style', 'edit', 'Position', [150, 200, 200, 20], ...
                  'String', '0.5', 'Tag', 'Psi_box');

% 确认按钮
uicontrol('Style', 'pushbutton', 'Position', [150, 150, 100, 30], ...
         'String', '确认', 'Callback', @submitCallback);

% 等待用户输入
uiwait(mainFig);

% 回调函数
function submitCallback(~, ~)
    % 获取输入框对象
    E_box = findobj(mainFig, 'Tag', 'E_box');
    B_box = findobj(mainFig, 'Tag', 'B_box');
    R_box = findobj(mainFig, 'Tag', 'R_box');
    Psi_box = findobj(mainFig, 'Tag', 'Psi_box');
    
    try
        % 获取输入值并计算
        E_val = str2double(get(E_box, 'String'));
        B_val = str2double(get(B_box, 'String'));
        R_val = str2double(get(R_box, 'String'));
        Psi_val = str2double(get(Psi_box, 'String'));
        
        gamma = E_val / 0.511;
        B = B_val * 1e-4;
        R_test = R_val * 1e3;
        Psi_test = Psi_val;
        
        % 关闭输入窗口
        uiresume(mainFig);
        close(mainFig);
        
        % 执行计算和绘图
        calculateAndPlot(gamma, B, R_test, Psi_test, E_val, B_val, R_val, Psi_val);
    catch ME
        errordlg(['输入错误: ' ME.message], '错误');
    end
end

function bracket = find_bracket(func, R, Psi, theta_min, theta_max, num)
    % 类似Python的find_bracket函数
    theta_test = linspace(theta_min, theta_max, num);
    f_test = arrayfun(@(theta) func(theta, R, Psi), theta_test);
    
    for i = 1:length(f_test)-1
        if f_test(i) * f_test(i+1) < 0
            bracket = [theta_test(i), theta_test(i+1)];
            return;
        end
    end
    bracket = [];
end

function calculateAndPlot(gamma, B, R_test, Psi_test, E_val, B_val, R_val, Psi_val)
    % 变量计算
    beta = sqrt(1 - 1/gamma^2); % 相对论速度，无量纲
    v0 = c * beta;  % 电子束速度，单位 [m/s]
    wB = qe * B / (gamma * me);  % 回旋频率
    a0_max = v0 / wB;  % 最大回转半径 [m]

    % 判断是否超出射程
    if R_test > 2*a0_max / sin(Psi_test)
        error("目标距离超出射程！当前目标角度下射程需小于: %.1f公里\n请选择距离或角度更小的目标后重试", 2*a0_max / sin(Psi_test)*1e-3);
    end

    % 求解发射角theta   
    a0_func = @(theta) a0_max * sin(theta);
    theta_equation = @(theta, R, Psi) R * sin(Psi) - a0_func(theta) * sqrt(2 * (1 - cos(wB * R * cos(Psi) / (v0 * cos(theta)))));
    
    % 使用find_bracket查找合适的区间
    bracket = find_bracket(theta_equation, R_test, Psi_test, 0.01, pi/2, 1000);
    
    if isempty(bracket)
        error("未找到合适的解区间，请检查输入参数");
    end
    
    % 使用fzero求解方程
    options = optimset('Display', 'off');
    theta_val = fzero(@(theta) theta_equation(theta, R_test, Psi_test), bracket, options);

    % 判断是否超出射程
    if R_test > 2*a0_max / sin(Psi_test)
        error("目标超出射程！最大射程为: %.1f米", 2*a0_max / sin(Psi_test));
    end

    tau_val = R_test * cos(Psi_test) / v0 / cos(theta_val);

    % 电子轨迹3D图
    a0_val = a0_func(theta_val);
    % 时间参数
    t = linspace(0, tau_val, 1000);
    % 轨迹方程
    x = a0_val * (1 - cos(wB * t));
    y = a0_val * sin(wB * t);
    z = v0 * cos(theta_val) * t;
    % 靶目标位置
    xT = a0_val * (1 - cos(wB * tau_val));
    yT = a0_val * sin(wB * tau_val);
    zT = v0 * cos(theta_val) * tau_val;
    % 坐标变换，y-y'夹角
    phi = atan2(-xT, yT);
    xp = cos(phi) .* x + sin(phi) .* y;
    yp = -sin(phi) .* x + cos(phi) .* y;
    xTp = cos(phi) * xT + sin(phi) * yT;
    yTp = -sin(phi) * xT + cos(phi) * yT;

    % 绘制3D图
    fig = figure;
    plot3(xp, yp, z, 'LineWidth', 2);
    hold on;
    plot3([0, xTp], [0, yTp], [0, zT], 'r--', 'LineWidth', 2);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title('三维电子打击靶目标示意图');
    grid on;
    view(3);

    % 保存图像
    saveas(fig, sprintf('电子束瞄准_R%.0f_Psi%.2f.jpg', R_test, Psi_test));

    % 准备数据并写入CSV文件
    header = {'电子束能量(MeV)', '磁场(g)', '目标距离(km)', ...
              '目标角度(rad)', '极角(deg)', '方位角(deg)', '打击时间(μs)'};
    
    % 四舍五入到1位小数
    data = {round(E_val, 1), round(B_val, 1), round(R_val, 1), ...
            round(Psi_test, 2), round(rad2deg(theta_val), 1), ...
            round(rad2deg(pi/2-phi), 1), round(tau_val*1e6, 1)};

    % 检查文件是否存在并写入
    filename = 'parameters.csv';
    if exist(filename, 'file')
        % 读取现有数据
        existingData = readtable(filename);
        % 追加新数据
        newRow = cell2table(data, 'VariableNames', header);
        updatedData = [existingData; newRow];
        writetable(updatedData, filename);
    else
        % 创建新文件
        writetable(cell2table(data, 'VariableNames', header), filename);
    end
end
end