% MATLAB代码：植株与地面点云分离
% 撰写日期：2025年3月12日
close all;
clear;
clc;

%% 1. 读取PLY点云数据
% pathname = 'E:\Datasets\3d\20241012\';
% filename = 'maize-60-Cloud-SOR-Octree12.ply';
% layerHeight = 0.02    茎秆层比例0.2     pKNN = 15   dbscan_minpts=minpts    stemDiameterThresh = 0.01;
% filename = 'maize-60-Cloud-SOR-plants222-optimized-Oct11.ply';
% filename = 'maize-60-Cloud-SOR-plants222-optimized-Oct11-Labeled-GT - Cloud.ply';
% layerHeight = 0.02    茎秆层比例0.2     pKNN = 25   dbscan_minpts=minpts   eps
% = avgNNDist * 1.9;    stemDiameterThresh = 0.5; cubeSide = 0.16;

% pathname = 'D:\MaizePCD\20250510a\';
% filename = '目标4-Noise-Octree11-NoChecker.ply';
% pKNN = 20;    dbscan_minpts=5

% pathname = 'E:\Datasets\3d\20241026\';
% filename = '目标42-Cloud-Cloud.ply';
% layerHeight = 0.02   茎秆层比例0.2    pKNN = 25;    dbscan_minpts=minpts   eps = avgNNDist * 0.9;
% filename = '目标42.las.segmented-Octree11-2.ply';
% pKNN = 60;    eps = avgNNDist * 1.9;  densityDiffThresh = 0.01 * globalDensity;
% OUTPUT '目标42.las.plants24152.ply'

pathname = 'E:\Datasets\3d\20241214\';
% filename = '目标4-Cloud-Octree10.ply';
filename = '目标4-Cloud-ClassOctree10-GT-labels01-align.ply';
% 目标4-Cloud-Noise+SOR_filter-Spatial0.01-Cloud2.ply
% layerHeight = 0.02   茎秆层比例0.1   pKNN = 25;    dbscan_minpts=1

% pathname = 'E:\Datasets\3d\20250422\';
% filename = '目标3-SOR+Clean-Cloud.ply';

ptCloud = pcread([pathname filename]); % 替换为你实际的点云文件路径
points = ptCloud.Location; % 获取点云坐标 [x, y, z]

%% 
plyFile = 'E:\Datasets\3d\20241026\目标42.las.segmented-Octree11-2.ply';
ptCloud_base = pcread(plyFile);

% Estimate Y shift
deltaY = estimateYShift_byHistogram(ptCloud_base, ptCloud, 0.05);

% Apply Y translation to ptCloud2
points = ptCloud.Location;
points(:,2) = points(:,2) + deltaY;

% 构造参数列表（按已有属性动态补齐）
args = {};

if ~isempty(ptCloud.Intensity)
    args = [args, {'Intensity', ptCloud.Intensity}];
end

if ~isempty(ptCloud.Color)
    args = [args, {'Color', ptCloud.Color}];
end

if ~isempty(ptCloud.Normal)
    args = [args, {'Normal', ptCloud.Normal}];
end

% 重新构造 pointCloud（属性完整保留）
ptCloud = pointCloud(points, args{:});


%% 绘制 intensity 直方图并寻找 Valley（最佳分割点）
ints = ptCloud.Intensity;

figure; 
histogram(ints, 500);
title('Intensity Histogram');
xlabel('Intensity');
ylabel('Count');
grid on;

%% 
% 1. intensity 转 [0,1]
% ints_norm = (ints - min(ints)) / (max(ints) - min(ints));
% 直方图均衡化
ints_norm = histeq(uint16(ints), 256);
ints_norm = double(ints_norm) / double(max(ints_norm));

% 2. 选择 colormap（jet/summer/parula/hot）
cmap = jet(256);

% 3. 将 intensity 映射到 RGB
idx = ceil(ints_norm * 255) + 1;
colors = cmap(idx, :);

% 4. 可视化
figure; 
scatter3(points(:, 1), points(:, 2), points(:, 3), 10, colors, 'filled');
title('Point Cloud Colored by Intensity');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% 2. 预处理：去噪和体素滤波
% % 统计滤波去噪
% meanK = 50; % 邻域点数
% stdThresh = 1; % 标准差阈值
% [~, inliers] = pcdenoise(ptCloud, 'NumNeighbors', meanK, 'Threshold', stdThresh);
% ptCloudClean = select(ptCloud, inliers);
% 
% % 体素滤波降低点云密度
% gridSize = 0.01; % 体素网格大小（单位：米）
% ptCloudDownsampled = pcdownsample(ptCloudClean, 'gridAverage', gridSize);
% points = ptCloudDownsampled.Location;
% 
% % 可视化
% figure;
% pcshow(points, [0.5 0.5 0.5], 'MarkerSize', 20); % 未分割点（灰色）

%% 3. 沿Z轴分层处理（优化方法1）
% 获取Z轴范围
zCoords = points(:, 3);
zMin = min(zCoords);
zMax = max(zCoords);
layerHeight = 0.02;     % 每层高度（单位：米，可调）
numLayers = ceil((zMax - zMin) / layerHeight);

% 分层并计算密度
layerPoints = cell(numLayers, 1);
layerDensities = zeros(numLayers, 1);
for i = 1:numLayers
    zLow = zMin + (i-1) * layerHeight;
    zHigh = zMin + i * layerHeight;
    idx = points(:, 3) >= zLow & points(:, 3) < zHigh;
    layerPoints{i} = points(idx, :);
    layerDensities(i) = size(layerPoints{i}, 1) / layerHeight;      % 密度 = 点数/层高
end

%% 绘制折线图
figure;
plot(layerDensities, 1:numLayers, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('层号 (layerNumber)');
ylabel('密度 (layerDensities)');
title('各层密度折线图 (Layer Density Line Chart)');
grid on;


%% 优化方法1：结合高度约束选择茎秆层
% stemHeightMid = zMin + (zMax - zMin) * 0.2;
% validLayerIdx = find(zMin + (0:numLayers-1) * layerHeight <= stemHeightMid);
% [~, minDensityIdx] = min(layerDensities(validLayerIdx));
% targetLayerIdx = validLayerIdx(minDensityIdx);
% targetLayerPoints = layerPoints{targetLayerIdx};
% 
% if isempty(targetLayerPoints)
%     error('目标层点云为空，请检查分层参数或点云数据');
% end

%% 
% ---- 1. 平滑密度曲线 ----
smoothedDensity = smooth(layerDensities, 3, 'moving'); % 移动平均平滑

[~, minimaLocs] = findpeaks(smoothedDensity);    % 波峰
if ~isempty(minimaLocs)
    groundEndIdx = minimaLocs(1); % 第一个波谷的位置
else
    % 回退策略：如取全局最小值或斜率平稳点
    [~, groundEndIdx] = min(smoothedDensity);
end

% % ---- 2. 计算一阶与二阶导数 ----
% d1 = diff(smoothedDensity);   % 一阶导，反映上升或下降趋势,[diff导致输出长度比输入少1，后续索引需对齐]
% d2 = diff(d1);                % 二阶导，反映曲率变化（拐点位置）[长度比d1再少1]
% 
% % ---- 3. 寻找显著下降后转平稳的拐点 ----
% % 下降段判定：d1 < 0
% % 拐点判定：二阶导数由负变正（即地面段结束）
% signChangeIdx = find(d2(1:end-1) < 0 & d2(2:end) > 0);
% 
% if ~isempty(signChangeIdx)
%     groundEndIdx = signChangeIdx(1) + 1; % +1 对齐原索引
% else
%     % 若未检测到明显拐点，取密度下降后首次稳定段
%     [~, groundEndIdx] = min(abs(d1)); 
% end

%
figure;
plot(layerDensities, 'r', 'LineWidth', 1.5); hold on;
plot(smoothedDensity, '--b', 'LineWidth', 2);
if exist('groundEndIdx', 'var')
    xline(groundEndIdx, '--k', 'LineWidth', 1);
    % 添加横向标签并调整样式
    text(groundEndIdx, layerDensities(groundEndIdx), '拐点', ...
        'HorizontalAlignment', 'center', ...
        'Rotation', 0, ...
        'FontSize', 14, ...
        'Color', 'g', ...
        'FontWeight', 'bold');
end
legend('原始数据', '平滑后', '拐点');
hold off;
xlabel('层号 (layerNumber)');
ylabel('密度 (layerDensities)');
title('各层密度折线图 (Layer Density Line Chart)');
grid on;



%% 

groundLayerPoints = layerPoints{groundEndIdx};

% 
% ---- 4. 确定茎秆基部层 ----
baseLayerIdx = groundEndIdx + 5;
baseLayerIdx = min(baseLayerIdx, numel(layerDensities));

% 获取相应层点云
targetLayerPoints = layerPoints{baseLayerIdx};

%% 可视化每一层的点云
figure;
hold on;
colors = jet(numLayers); % 使用不同颜色表示不同的层
for i = 1:numLayers
    scatter3(layerPoints{i}(:, 1), ...
             layerPoints{i}(:, 2), ...
             layerPoints{i}(:, 3), ...
             1, ...
             'MarkerFaceColor', colors(i,:), ...
             'MarkerEdgeColor', colors(i,:) ...
             );
end

% 在原始点云中，可视化提取出来的特殊层
scatter3(groundLayerPoints(:, 1), groundLayerPoints(:, 2), groundLayerPoints(:, 3), 15, 'filled');

scatter3(targetLayerPoints(:, 1), targetLayerPoints(:, 2), targetLayerPoints(:, 3), 15, 'filled');

title([num2str(numLayers) ' 层点云可视化' ]);
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off;


%% 4. 投影聚类定位茎秆（优化方法2）
% 投影到XY平面
targetLayerXY = targetLayerPoints(:, 1:2);

% 自适应DBSCAN参数
% 计算每个点到最近邻的距离，取平均值
pKNN = 60;
distances = pdist2(targetLayerXY, targetLayerXY, 'euclidean', 'Smallest', pKNN);   % 每个点的 pKNN 个最近邻距离
avgNNDist = mean(distances(2:pKNN, :), 'all');    % 只取最近邻（第1近是自己，距离为0）
eps = avgNNDist * 1.9;                           % 动态调整邻域半径
localDensities = sum(distances < eps, 2);
minPts = min(2, round(median(localDensities) * 0.03));     % 自适应minPts

% DBSCAN聚类, 第三个参数{V6=5, V8=10, V11=20}
[clusterIdx, C] = dbscan(targetLayerXY, 0.38, 10);

% 可视化结果
% figure;
% gscatter(targetLayerXY(:,1), targetLayerXY(:,2), clusterIdx);
% hold on;
% legend('Cluster 1', 'Cluster 2', 'Noise', 'Core points', 'Location', 'Best');
% title('DBSCAN Clustering');
% hold off;

% 移除噪声点
validIdx = clusterIdx ~= -1;
filteredXY = targetLayerXY(validIdx, :);
filteredLabels = clusterIdx(validIdx);

numClusters = max(filteredLabels);

% ==== 计算每簇中心点（向量化，避免循环） ====
clusterCenters = accumarray(filteredLabels, filteredXY(:,1), [], @mean);
clusterCenters(:,2) = accumarray(filteredLabels, filteredXY(:,2), [], @mean);

% 可视化 XY 平面投影和聚类结果
fig = figure(Position=[100, 100, 1200, 600]);
scatter(targetLayerXY(:, 1), targetLayerXY(:, 2), 5, 'k', 'filled');       % 绘制所有投影点（黑色）
hold on;
colors = lines(numClusters);        % 为每个簇分配不同颜色
for i = 1:numClusters
    clusterPoints = targetLayerXY(clusterIdx == i, :);
    scatter(clusterPoints(:, 1), clusterPoints(:, 2), 10, colors(i, :), 'filled'); % 绘制聚类点
    
    scatter(clusterCenters(i, 1), clusterCenters(i, 2), 50, 'r', 'x', 'LineWidth', 2); % 绘制聚类中心（红色X）
    % 显示标签信息
    text(clusterCenters(i, 1)+0.1, clusterCenters(i, 2), int2str(i), 'Color', 'red', 'FontSize', 14);
end

xlabel('X (地面东西方向)', 'FontSize', 14);
ylabel('Y (地面南北方向)', 'FontSize', 14);

title(['玉米群体点云XY平面投影与聚类结果：' num2str(numClusters) '株'], 'FontSize', 14);
grid on;
hold off;

%% 保存聚类中心
[~, name, ~] = fileparts(filename);

outputPath = fullfile(pathname, [name, '-clusterCenters.mat']);
save(outputPath, 'clusterCenters')

%% 5. 生长分割分离植株与地面（优化版）
% 初始化分割结果
labels = zeros(size(points, 1), 1); % 0=地面（默认），1=植株

% 目标层高度下界（从第3步分层逻辑）
targetLayerLowerBound = zMin + (baseLayerIdx - 1) * layerHeight; % 目标层下边界

% 将高于目标层下边界的点先标记为植株
plantMask = points(:, 3) > targetLayerLowerBound;
labels(plantMask) = 1;

% 初始化剩余点池
remainingPoints = points(~plantMask, :);
remainingIdx = find(~plantMask);

% 全局密度均值
globalDensity = mean(layerDensities);

% 定义立方体边长和密度差阈值
cubeSide = 0.4;            % 立方体边长（米，可调，控制地面附着点大小，越小附着点越小）
densityDiffThresh = 0.01 * globalDensity;   % 控制地面附着点大小，比例越小附着点越小，但根部越短，比例大根部分割的越长
stemDiameterThresh = 0.5;
radius = cubeSide / 4; % 搜索半径

for i = 1:numClusters
    % 从目标层过滤后的点中找到聚类中心的最近点
    seedXY = clusterCenters(i, :);
    seedIdxInLayer = knnsearch(targetLayerXY, seedXY, 'K', 1);
    seedPoint = mean(targetLayerPoints(seedIdxInLayer, :), 1);

    % 在剩余点中找到对应索引
    seedIdx = knnsearch(remainingPoints, seedPoint, 'K', 1);

    % 初始化生长队列
    queue = seedIdx;
    visited = false(size(remainingPoints, 1), 1);
    visited(seedIdx) = true;
    labels(remainingIdx(seedIdx)) = 1;
    
    prevDensity = 0; % 初始上一密度
    
    while ~isempty(queue)
        % 取当前生长点
        currentIdx = queue(1);
        queue(1) = []; % 出队

        currentPoints = remainingPoints(currentIdx, :);
        
        % 邻域查找（局部向量化）
        [neighborsCell, ~] = rangesearch(remainingPoints, currentPoints, radius);
        neighbors = unique([neighborsCell{:}]');    % 去重
        
        % 向下生长：Z 必须小于当前点，并且未访问过
        % 筛选Z值小于currentZ的点
        downwardMask = remainingPoints(neighbors, 3) < currentPoints(3) & ~visited(neighbors);
        newIdx = neighbors(downwardMask);

        if isempty(newIdx)
            break; % 无向下点可添加，停止生长
        end
        
        % 标记新点为植株
        visited(newIdx) = true;
        labels(remainingIdx(newIdx)) = 1;
        
        % --- 局部密度计算（稳健版） ---
        % 计算新添加点的局部密度
        newPoints = remainingPoints(newIdx, :);
        localRange = max(newPoints, [], 1) - min(newPoints, [], 1);

        % 避免过小体积，至少取 cubeSide 的体积
        localVolume = prod(max(localRange, cubeSide)); % 动态体积
        localDensity = size(newPoints, 1) / localVolume;
        
        % 判断密度变化        
        densityDiff = localDensity - prevDensity;
        if densityDiff > densityDiffThresh
            break;
        end

        % 横截面直径限制
        xyDiameter = max(pdist(newPoints(:,1:2)));
        if xyDiameter > stemDiameterThresh
            break;
        end

        prevDensity = localDensity;

        queue = [queue; newIdx(:)];
    end
end

%% 6. 可视化结果
groundCloud = select(ptCloud, find(labels == 0)).Location;
plantCloud = select(ptCloud, find(labels == 1)).Location;

figure('Position',[200 150 900 600]);
axNew = gca;
scatter3(groundCloud(:, 1), groundCloud(:, 2), groundCloud(:, 3), 3, [1 0 0]);   % 地面（绿色）
hold on;
scatter3(plantCloud(:, 1), plantCloud(:, 2), plantCloud(:, 3), 5, [0 1 0], "filled");   % 植株（红色）

% for i = 1:numClusters    
%     scatter(clusterCenters(i, 1), clusterCenters(i, 2), 50, 'b', 'x', 'LineWidth', 2); % 绘制聚类中心（红色X）
%     % 显示标签信息
% %     text(clusterCenters(i, 1), clusterCenters(i, 2), int2str(i), 'Color', 'red', 'FontSize', 14);
% end

% 载入参数并应用到新坐标区
load('myView.mat','viewParam');   % 结构体 viewParam 回到工作区

axNew.CameraPosition  = viewParam.CameraPosition;
axNew.CameraTarget    = viewParam.CameraTarget;
axNew.CameraUpVector  = viewParam.CameraUpVector;
axNew.XLim = viewParam.XLim;
axNew.YLim = [-2.5, 2.5];
axNew.ZLim = viewParam.ZLim;

xlabel('X(地面东西方向)', 'FontSize', 14); 
ylabel('Y(地面南北方向)', 'Rotation', -45, 'FontSize', 14); 
zlabel('Z(地上竖直方向)', 'FontSize', 14);

title('群体点云地面分割结果', 'FontSize', 14);
legend({'地面（红色）', '植株（绿色）'}, ...
       'Units','normalized', ...     % 用归一化坐标，随窗口等比
       'Position',[0.75 0.86 0.15 0.06]); % [left bottom width height]

%% 获取图形属性
% 把当前坐标区的视图参数打包存盘
viewParam.CameraPosition  = get(gca,'CameraPosition');
viewParam.CameraTarget    = get(gca,'CameraTarget');
viewParam.CameraUpVector  = get(gca,'CameraUpVector');
viewParam.XLim            = get(gca,'XLim');
viewParam.YLim            = get(gca,'YLim');
viewParam.ZLim            = get(gca,'ZLim');

% save('myView.mat','viewParam');   % 自动存到当前工作目录

%% 7. 保存分割结果
% 地面点 labels == 0
% 植株点 labels == 1

% 1. 将 labels 转换为 Intensity（确保 labels 是列向量）
if size(labels, 2) > 1
    labels = labels'; % 转置为列向量（N×1）
end

% 2. 创建新点云并赋值 Intensity（若原点云无 Intensity 字段）
ptCloudWithIntensity = pointCloud(ptCloud.Location, 'Intensity', labels);

% 3. 保存为 PLY 文件（包含 Intensity 属性）
[~, name, ~] = fileparts(filename);
pcwrite(ptCloudWithIntensity, [pathname name '-withPlantsLables-0.4.ply'], 'PLYFormat', 'binary');

% 可选：验证保存的属性
disp(['保存的点云包含 Intensity 属性: ', num2str(isfield(ptCloudWithIntensity.Intensity, 'Intensity'))]);

%%
plantCloud = select(ptCloud, find(labels == 1));
pcwrite(plantCloud, [pathname '目标42.las.plants24152.ply']);

%% 地面点
groundCloud = select(ptCloud, find(labels == 0));
pcwrite(groundCloud, [pathname '目标42.las.plants24152-ground.ply']);

%% 带标签保存
% ptCloudWithLabels = pointCloud(ptCloud.Location, 'Intensity', labels-1);
% pcwrite(ptCloudWithLabels, [pathname 'plants241-withLables.ply']);
% 此方法保存的标签点云数据精度较低
% 直接分开保存植株点云和地面点云，之后通过CloudCompare添加标签（注意命名）后合并，
% 之后打开PLY格式文件，将 scalar_intensity（或添加标签时的名称修改为 intensity，因为 pcread 函数无法读取非标准属性）

%% 辅助函数
function center = computeClusterCenters(XY, method)
% 计算点簇的中心点
% XY     : N x 2 点集
% method : 'mean' | 'median' | 'pca'
%
% 返回值:
%   center : 1x2 的中心坐标

    arguments
        XY (:,2) double
        method (1,:) char {mustBeMember(method,{'mean','median','pca'})} = 'mean'
    end

    switch method
        case 'mean'
            % ==== 方法 1: 均值中心 ====
            center = mean(XY, 1);

        case 'median'
            % ==== 方法 2: 中位数中心 (Weiszfeld 算法) ====
            tol = 1e-6;  
            maxIter = 100;
            x0 = median(XY,1);   % 初始点
            for iter = 1:maxIter
                dists = sqrt(sum((XY - x0).^2, 2));
                if any(dists == 0)   % 防止除零
                    x0 = XY(dists == 0, :);
                    break;
                end
                w = 1 ./ dists;
                x1 = sum(XY .* w, 1) / sum(w);
                if norm(x1 - x0) < tol
                    break;
                end
                x0 = x1;
            end
            center = x0;

        case 'pca'
            % ==== 方法 3: PCA 主轴中心 ====
            [coeff, score, ~] = pca(XY);
            mainAxis = coeff(:,1);
            proj = score(:,1);
            midProj = (min(proj) + max(proj)) / 2;
            center = mean(XY,1) + midProj * mainAxis';
    end
end

function [center, radius] = circfit(x,y)
    % CIRCFIT fits a circle to given x,y data in the least square sense.
    % [center, radius] = circfit(x,y)
    % x,y: column vectors of points (x(i),y(i))
    % center: (x0,y0) coordinate of the center
    % radius: radius of the fitted circle

    x = x(:);
    y = y(:);
    a = [x y ones(size(x))] \ (x.^2 + y.^2);
    center = [a(1)/2 a(2)/2];
    radius = sqrt((a(3) + center(1)^2 + center(2)^2));
end

function deltaY = estimateYShift_byHistogram(ptCloud1, ptCloud2, binWidth)
% Estimate Y-axis translation between two point clouds
% using 1D histogram cross-correlation
%
% ptCloud1, ptCloud2 : pointCloud objects
% binWidth           : histogram bin size (e.g. 0.02 or 0.05)
%
% deltaY > 0 means ptCloud2 should be shifted +deltaY in Y

if nargin < 3
    binWidth = 0.05;
end

Y1 = ptCloud1.Location(:,2);
Y2 = ptCloud2.Location(:,2);

% Use common Y range
yMin = min([Y1; Y2]);
yMax = max([Y1; Y2]);

edges = yMin:binWidth:yMax;

% 1D histograms (density normalized)
h1 = histcounts(Y1, edges, 'Normalization', 'pdf');
h2 = histcounts(Y2, edges, 'Normalization', 'pdf');

% Remove mean to emphasize shape
h1 = h1 - mean(h1);
h2 = h2 - mean(h2);

% Cross-correlation
[c, lags] = xcorr(h1, h2);

% Best lag
[~, idx] = max(c);
bestLag = lags(idx);

% Convert lag to physical shift
deltaY = bestLag * binWidth;
end


