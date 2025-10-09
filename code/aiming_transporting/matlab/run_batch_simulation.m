function run_batch_simulation(csv_path)
% 批量执行电子束仿真

% 读取CSV文件
try
    params_table = readtable(csv_path);
catch ME
    error('无法读取CSV文件: %s', ME.message);
end

% 获取参数列名
param_names = params_table.Properties.VariableNames;

% 遍历每一组参数
for idx = 1:height(params_table)
    try
        % 准备参数结构体
        raw_params = table2struct(params_table(idx,:));
        params = prepare_params(raw_params);
        
        % 执行仿真
        run_simulation(params);
        
    catch ME
        fprintf('\n第%d组参数仿真失败:\n', idx);
        disp('参数:');
        disp(raw_params);
        fprintf('错误详情: %s\n', ME.message);
        fprintf('错误堆栈:\n');
        for k = 1:length(ME.stack)
            fprintf('  %s (行%d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('----------------------------------------\n');
    end
end
end