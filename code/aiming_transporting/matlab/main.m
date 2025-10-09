%% 电子束瞄准与注量率估算主程序
clear; clc; close all;

% 获取用户输入
user_inputs = get_user_inputs();
if isempty(user_inputs)
    disp('未获取到输入参数，程序终止。');
    return;
end

if isfield(user_inputs, 'mode') && strcmp(user_inputs.mode, 'batch')
    run_batch_simulation(user_inputs.csv_path);
else
    params = prepare_params(user_inputs);
    run_simulation(params);
end