% MATLAB Ground Filtering Evaluation Visualization
clear; clc; close all;

% 1. 读取数据 (假设数据保存为 data.csv)
% 建议CSV格式：Stage, Name, Precision, Recall, IoU_p, F1, RNE_F1, VD_Plant

% CSV 文件路径
% csvFile = 'C:\Users\BaseBoot\Desktop\论文\小论文1\groundPlantSegRes.csv';
csvFile = 'E:\Datasets\3d\result\groundSegResult.csv';
%
opts = detectImportOptions(csvFile);
data = readtable(csvFile, opts);

% 提取核心分析指标
methods = {'Ours', 'CSF', 'SMRF'};
stages = {'V6', 'V8', 'V11'};
metrics = {'IoU_p', 'F1', 'RNE_F1', 'VD_p', 'Recall_p'};

%% 2. 数据归一化处理
normData = data;
for i = 1:length(metrics)
    mName = metrics{i};
    colData = data.(mName);
    if strcmp(mName, 'VD_p')
        % 负向指标归一化：(max - x) / (max - min)
        normData.(mName) = 1 ./ (1 + data.VD_p);
    end
end

%% 方案 A: 组合柱状图 (IoU_p 与 RNE_F1)
figure('Color', 'w', 'Name', 'Scheme A');
targetMetrics = {'Precision', 'Recall_p', 'F1', 'IoU_p'};
for i = 1:length(targetMetrics)
    subplot(1, length(targetMetrics), i);
    % 重组数据用于绘图 [Methods x Stages]
    plotVec = zeros(length(stages), length(methods)); 
    for s = 1:length(stages)
        for m = 1:length(methods)
            idx = strcmp(data.Stage, stages{s}) & strcmp(data.Method, methods{m});
            plotVec(m, s) = data.(targetMetrics{i})(idx);
        end
    end
    b = bar(plotVec');
    set(gca, 'XTickLabel', stages);
    title([targetMetrics{i}], 'Interpreter','none');
    ylabel('Value');
    legend(methods, 'Location', 'best');
    grid on;
    % 缩小Y轴范围以突出差异
    if i == 1
        ylim([0.98, 1]);
    else
        ylim([0.9 1]);
    end
end

%% 方案 B: 归一化后的综合得分 (箱线图修复版)
figure('Color', 'w', 'Name', 'Scheme B');
combinedScores = [];
groupLabels = {};

for m = 1:3
    mIdx = strcmp(data.Method, methods{m});
    % 提取该方法在所有指标(metrics)上的归一化分值
    mScores = data{mIdx, metrics}; 
    mScores_flat = mScores(:); % 转为一维列向量
    
    combinedScores = [combinedScores; mScores_flat];
    
    % 修复 repmat 逻辑：复制因子直接使用 [行数, 1]
    % 产生与分数向量等长的标签向量
    currentLabels = repmat(methods(m), size(mScores_flat, 1), 1);
    groupLabels = [groupLabels; currentLabels];
end

boxplot(combinedScores, groupLabels);
ylabel('Normalized Score (0-1)');
title('Algorithm Comprehensive Performance (Normalized)');
grid on;

%% 方案 C: 关键指标折线图 (RNE_F1 随生育期的变化)
figure('Color', 'w', 'Name', 'Scheme C');
hold on;
colors = {'#D95319', '#0072BD', '#EDB120'}; % 区分 Ours, CSF, SMRF
for m = 1:3
    yVal = [];
    for s = 1:3
        idx = strcmp(data.Stage, stages{s}) & strcmp(data.Method, methods{m});
        yVal(s) = data.RNE_F1(idx);
    end
    plot(1:3, yVal, '-o', 'LineWidth', 2, 'Color', colors{m}, 'MarkerFaceColor', colors{m});
end
set(gca, 'XTick', 1:3, 'XTickLabel', stages);
xlabel('Growth Stage'); ylabel('RNE\_F1');
title('Trend of RNE\_F1 over Stages');
legend(methods); grid on;

%% 方案 D: 雷达图 (以 V11 期为例)
figure('Color', 'w', 'Name', 'Scheme D', 'Position', [600, 350, 850, 550]);

s = 2;
v11Idx = strcmp(data.Stage, stages{s});
v11Data = data(v11Idx, :);
% 准备雷达图指标 (归一化后的数据)
radarMetrics = {'Precision', 'Recall_p', 'F1', 'IoU_p', 'RNE_F1', 'VD_p'};
P = [];
for m = 1:3
    mIdx = strcmp(v11Data.Method, methods{m});
    P(m, :) = v11Data{mIdx, radarMetrics};
end
% 使用第三方或自定义雷达图函数，此处用极坐标简易模拟
theta = linspace(0, 2*pi, length(radarMetrics) + 1);
for m = 1:3
    polarplot(theta, [P(m,:) P(m,1)], 'LineWidth', 2);
    hold on;
end

% 设置刻度和标签
% 1. 径向范围锁定 0.8~1
% set(gca, ...
%     'RLim', [0.8 1], ...          % 画布只留 0.8–1
%     'RTickLabel', string(0.8:0.05:1));

set(gca, ...
    'ThetaTick', rad2deg(theta(1:end-1)), ...
    'ThetaTickLabel', radarMetrics, ...
    'TickLabelInterpreter', 'none', ...
    'FontSize', 14);                      % ← 统一设置字体大小

title(sprintf('Performance Comparison at %s (Normalized)', stages{s}));
legend(methods);




%%
% MATLAB PLY Point Cloud Segmentation Matrix Visualization
clear; clc; close all;

% 1. 配置参数
stages = {'V6', 'V8', 'V11'};
methods = {'GT', 'Ours', 'CSF', 'SMRF'}; 

ground_color = [0.0 1.0 0.0]; % 地面：中灰色
plant_color = [0.0 0.7 0.0];  % 植株：翠绿色

base_dir = 'E:\Datasets\3d\groundPlantSegResultVisualize\';

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.85]);

% 2. 循环遍历并绘图
% num_method = length(methods);
num_method = 2;
for s = 1:length(stages)
    for m = 1:num_method
        % 动态构建文件名 (假设命名格式如: V11_Ours.ply)
        fileName = fullfile(base_dir, sprintf('%s_%s.ply', stages{s}, methods{m}));
        
        if exist(fileName, 'file')
            % 使用 pcread 读取 PLY 文件
            ptCloud = pcread(fileName);
            xyz = ptCloud.Location;
            
            % 关键步骤：提取 Intensity 属性作为分类标签
            % 假设 0 = Ground, 1 = Plant
            labels = ptCloud.Intensity;
            
            % 计算子图位置
            ax = subplot(length(stages), num_method, (s-1)*num_method + m);
            
            % 提取地面点和植株点索引
            idx0 = (labels == 0);
            idx1 = (labels == 1);
            
            % 绘制地面点
            scatter3(xyz(idx0,1), xyz(idx0,2), xyz(idx0,3), 2, ground_color, 'filled', 'MarkerEdgeAlpha', 0.2);
            hold on;
            % 绘制植株点 (点径略大，突出主体)
            scatter3(xyz(idx1,1), xyz(idx1,2), xyz(idx1,3), 6, plant_color, 'filled');
            
            % 视觉一致性优化
            if s == 1, title(methods{m}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none'); end
            
            % 2. 每一列最左侧注明生育期 (Row Header)
            if m == 1
                % 使用 text 函数在子图左侧添加垂直文本，避免与 ylabel 坐标轴重叠
                ylh = ylabel(stages{s}, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
                % 调整位置，使其更靠近左侧边缘
                set(ylh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]); 
            end
            
            % 视角锁定（这对 SCI 论文对比至关重要）
            view(45, 15); 
            axis equal;
            axis off; % 隐藏坐标轴增加美观度
            grid off;
        else
            fprintf('未找到文件: %s\n', fileName);
        end
    end
end

% 添加全局标题
sgtitle('Maize Segmentation Results (Grey: Ground, Green: Plant)', 'FontSize', 16);

% 添加统一图例
h_dummy1 = scatter3(nan, nan, nan, 20, ground_color, 'filled');
h_dummy2 = scatter3(nan, nan, nan, 20, plant_color, 'filled');
legend([h_dummy1, h_dummy2], {'Ground', 'Maize Plant'}, 'Location', 'southoutside', 'Orientation', 'horizontal');

%% 根部放大图
% MATLAB Maize Root Detail Comparison (PLY Version)
clear; clc; close all;

% 1. 配置
stages = 'V6'; 
methods = {'GT', 'Ours', 'CSF', 'SMRF'};
root_z_limit = 0.5; % 重要：根据您的点云尺度调整，只保留地面以上0.5米范围内的点
ground_color = [0.8 0.8 0.8]; % 地面淡灰色
plant_color = [0.0 0.8 0.0];  % 植株鲜绿色

slice_z_min = 0.04; % 离地 5cm
slice_z_max = 0.07; % 离地 15cm

base_dir = 'E:\Datasets\3d\groundPlantSegResultVisualize\';

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.1, 0.8, 0.4]);
hold on;

for m = 1:length(methods)
    fileName = fullfile(base_dir, sprintf('%s_%s.ply', stages, methods{m}));
    if ~exist(fileName, 'file'), continue; end
    
    ptCloud = pcread(fileName);
    xyz = ptCloud.Location;
    labels = ptCloud.Intensity;
    
    % 提取特定高度的植株点 (Label=1)
    % 这能展示茎秆被切除后的横截面完整度
    slice_idx = (xyz(:,3) >= slice_z_min & xyz(:,3) <= slice_z_max);
    s_xyz = xyz(slice_idx, :);
    
    subplot(1, 4, m);
    % 俯视图绘制
    scatter(s_xyz(:,1), s_xyz(:,2), 15, [0 0.6 0], 'filled', 'MarkerFaceAlpha', 0.5);
    
    title(methods{m}, 'FontSize', 15);
%     axis equal; 
    % 统一坐标轴范围，确保圆饼大小可比
%     xlim([-0.15, 0.15]); ylim([-0.15, 0.15]); 
    grid on; box on;
    xlabel('X (m)'); if m==1, ylabel('Y (m)'); end
end

sgtitle(['Root Cross-section Comparison (Z: ', num2str(slice_z_min), '-', num2str(slice_z_max), 'm)'], 'FontSize', 16);
hold off;


%% 
% MATLAB Vertical Point Distribution Analysis
clear; clc; close all;

stages = 'V6'; 
methods = {'GT', 'Ours', 'CSF', 'SMRF'};
colors = {'#000000', '#D95319', '#0072BD', '#EDB120'}; % 黑(GT), 橙(Ours), 蓝, 黄
z_step = 0.01; % 分层步长：1cm
slice_z_min = 0.0001; 
slice_z_max = 0.3; % 离地 15cm

base_dir = 'E:\Datasets\3d\groundPlantSegResultVisualize\';

figure('Color', 'w', 'Name', 'Vertical Distribution', 'Position', [280 410 1020 420]);
hold on;

for m = 1:length(methods)
     fileName = fullfile(base_dir, sprintf('%s_%s.ply', stages, methods{m}));
    if ~exist(fileName, 'file'), continue; end
    
    ptCloud = pcread(fileName);
    xyz = ptCloud.Location;
    labels = ptCloud.Intensity;
    
    % 只统计植株点 (Label == 1)
    p_xyz = xyz(labels == 1, :);
    
    % 分层统计
    z_bins = slice_z_min:z_step:slice_z_max;
    counts = zeros(length(z_bins)-1, 1);
    bin_centers = z_bins(1:end-1) + z_step/2;
    
    for b = 1:length(z_bins)-1
        counts(b) = sum(p_xyz(:,3) >= z_bins(b) & p_xyz(:,3) < z_bins(b+1));
    end
    
    % 绘制曲线 (横轴点数，纵轴高度)
    plot(bin_centers, counts, 'Color', colors{m}, 'LineWidth', 2, 'DisplayName', methods{m});
end

% 美化图表

% 1. X轴（高度Z）设置：去掉科学计数法，设置均匀刻度
ax = gca;
ax.XAxis.Exponent = 0; % 强制关闭科学计数法
xtickformat('%.2f');   % 设置显示两位小数 (例如 0.10, 0.20)

% 手动设置 X 轴均匀刻度（例如每 0.1m 一个刻度）
xticks(0:0.05:slice_z_max); 
xlim([0, slice_z_max]);

% 2. Y轴（点数）设置
% 如果点数太多导致显示科学计数法，同样关闭
ax.YAxis.Exponent = 0; 
ytickformat('%d');     % 整数显示点数

xlabel('Height Z (m)', 'FontSize', 14);
ylabel('Number of Plant Points', 'FontSize', 14);
title([stages, ': Vertical distribution of plant points (Z: 0 - ', num2str(slice_z_max), ' m)'], 'FontSize', 16);

legend('Location', 'northeast', ...
       'Units','normalized', ...     % 用归一化坐标，随窗口等比
       'Position',[0.73 0.68 0.09 0.22]); % [left bottom width height]

grid on;
set(gca, 'XScale', 'log'); % 如果点数差异巨大，建议使用对数坐标