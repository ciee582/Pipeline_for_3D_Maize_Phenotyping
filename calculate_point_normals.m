%% ==============================================
%  功能：批量读取点云 (XYZC)，计算方向向量 (DXDYDZ)，并保存为 (XYZDXDYDZC)
%  日期：2025-10-13
%  ==============================================
clc; clear vars;

basePath = 'E:\Datasets\3d\myMaize-leaf-stem\New\';
filePath = '1214';
inputFolder = fullfile(basePath, filePath);
outputFolder = fullfile(basePath, [filePath '+normal']);

% 邻域大小
k = 5;  

mode = 'principal'; % 'principal' 主方向 / 'normal' 法向方向

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

files = dir(fullfile(inputFolder, '*.txt'));
% files = dir(fullfile(inputFolder, '08-2_XYZC.txt'));

for f = 1:length(files)
    filePath = fullfile(files(f).folder, files(f).name);
    data = readmatrix(filePath);  % 每行格式: X Y Z C
    pts = data(:, 1:3);
    C = data(:, 4);

    % 构建KDTree
    Mdl = KDTreeSearcher(pts);
    dirs = zeros(size(pts));

%     fprintf('正在处理文件 %d/%d: %s ...\n', f, length(files), files(f).name);

    % --- 计算方向向量 ---
    for i = 1:size(pts, 1)
        idx = knnsearch(Mdl, pts(i,:), 'K', k);
        nbrs = pts(idx, :);

        % PCA特征分解
        COV = cov(nbrs);
        [V, D] = eig(COV);
        eigVals = diag(D);

        switch mode
            case 'principal'
                [~, maxIdx] = max(eigVals);
                dirVec = V(:, maxIdx);   % 主方向向量
            case 'normal'
                [~, minIdx] = min(eigVals);
                dirVec = V(:, minIdx);   % 法向量
            otherwise
                error('mode 参数只能为 "principal" 或 "normal"');
        end

        dirVec = dirVec / norm(dirVec); % 单位化
        dirs(i, :) = dirVec;
    end

    % --- 合并输出 ---
    outputData = [pts, dirs, C];
    outPath = fullfile(outputFolder, files(f).name);
    writematrix(outputData, outPath, 'Delimiter', ' ');
    fprintf('→ 已保存至 %s\n', outPath);
end

fprintf('\n所有文件已处理完成！\n');


%% --- 可视化 ---
filename = fullfile(outputFolder, '02-1_GT_XYZC_nobg.txt');

% 读取数据：每行 [X Y Z Nx Ny Nz Label]
data = load(filename);
pts   = data(:, 1:3);          % XYZ 坐标
norms = data(:, 4:6);          % 法向量 Nx Ny Nz
label = data(:, 7);            % Label (0~4)

% 归一化法向量
norms = norms ./ vecnorm(norms, 2, 2);  % 每行归一化为单位向量

% 根据 label 生成颜色（0~4 共5类）
colors = lines(5);             % 使用 MATLAB 内置 colormap
C = colors(label + 1, :);      % label 是 0~4，需 +1 转为 1~5 索引

% 创建图形窗口
figure('Name', sprintf('点云与法向量可视化 - %s', filename), 'Color', 'w');
hold on;

% 绘制点云（按标签着色）
scatter3(pts(:,1), pts(:,2), pts(:,3), 10, C, 'filled');

% 抽稀绘制法向量:为了清晰，只绘制部分法向量（例如每隔10个点）
step = 10;
idx = 1:step:size(pts,1);
% 绘制法向量（quiver3）
quiver3(pts(idx,1), pts(idx,2), pts(idx,3), ...
        norms(idx,1), norms(idx,2), norms(idx,3), ...
        0.5, 'k', 'LineWidth', 1);  % 黑色箭头，长度缩放0.5

% 设置坐标系
axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('点云 + 法向量 + 标签着色', 'FontSize', 12);

% 添加图例
% 生成每个 Label 的图例项
numLabels = 5;
h_legend = gobjects(numLabels, 1);  % 预分配句柄数组
for k = 1:numLabels
    h_legend(k) = scatter3(NaN, NaN, NaN, 10, colors(k,:), 'filled');
end

% 构造图例字符串：'Label 0', 'Label 1', ..., 'Label 4'
legendLabels = compose('Label %d', 0:(numLabels-1));

% 添加图例（仅点云类别，不包含法向量）
legend(h_legend, legendLabels, 'Location', 'bestoutside');

% 刷新图形
drawnow;
