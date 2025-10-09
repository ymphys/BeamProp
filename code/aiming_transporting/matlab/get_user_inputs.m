function user_inputs = get_user_inputs()
% MATLAB参数输入界面

% 创建主对话框
fig = figure('Name', '电子束参数输入', 'NumberTitle', 'off', ...
    'Position', [500, 500, 400, 500], 'MenuBar', 'none', 'ToolBar', 'none');

% 默认参数值
defaults = struct('E', '35', 'B', '0.25', 'R', '15', 'Psi', '0.5', ...
    'n_i', '1e5', 'Ib', '1', 'em_type', 'Ti', 'R0', '5');

% 创建输入控件
uicontrol('Style', 'text', 'Position', [20, 450, 150, 20], ...
    'String', '电子束能量 [MeV]:', 'HorizontalAlignment', 'left');
E_edit = uicontrol('Style', 'edit', 'Position', [200, 450, 150, 20], ...
    'String', defaults.E);

% 添加其他参数控件...(类似上面)

% 添加批量处理按钮
batch_btn = uicontrol('Style', 'pushbutton', 'Position', [150, 50, 100, 30], ...
    'String', '批量处理', 'Callback', @batch_callback);

% 添加确认按钮
ok_btn = uicontrol('Style', 'pushbutton', 'Position', [150, 10, 100, 30], ...
    'String', '确认', 'Callback', @ok_callback);

% 等待用户操作
uiwait(fig);

% 回调函数
    function batch_callback(~,~)
        [file, path] = uigetfile('*.csv', '选择批量参数文件');
        if file
            user_inputs.mode = 'batch';
            user_inputs.csv_path = fullfile(path, file);
            delete(fig);
        end
    end

    function ok_callback(~,~)
        % 参数验证和收集...(类似Python版本)
        delete(fig);
    end
end