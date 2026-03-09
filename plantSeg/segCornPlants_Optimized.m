clc;
clearvars;

%% 加载数据
% pathname = 'E:\Datasets\3d\20241012\';
% filename = 'plants223.ply';
% filename = 'plants240.ply';
% filename = 'plants241-withLables.ply';
% closedImg 5

pathname = 'E:\Datasets\3d\20241026\';
% filename = 'plants222.ply';
% filename = 'plants2414-optimized(131415-4).ply';
filename = '目标42.las.plants2415-labels.ply';

% pathname = 'E:\Datasets\3d\20241214\';
% filename = '目标4-Cloud-plants-01-withLables-SOR-segment-GT.ply';
% filename = '目标4-Cloud-plants-01-withLables-SOR-segment-GT - Figure.ply';

% 分组数
groupSize = 4;

PLYPoints = pcread([pathname filename]);

pNums = PLYPoints.Count;

%% 按 Intensity 属性重新读取并提取对应标签的点云数据
% 提取 Intensity 属性（即 labels）
ints = PLYPoints.Intensity;

% 逻辑索引,定义目标标签
mask = (ints == 1);

% 提取对应点云
targetPoints = select(PLYPoints, find(mask));

% 可视化
% figure;
% pcshow(targetPoints);

Points = double(targetPoints.Location);  % 将点云转换为双精度类型
ptCloud = Points(:, 1:3);  % 提取XYZ坐标

%% 投影点密度图像生成
% --- 优化 1: 定义分辨率参数 ---
GridResolution = 0.001;         % 定义栅格单元的物理尺寸 (例如 0.005 米/像素)

% 点云数据 X 和 Y 方向上的范围
x_range = [min(ptCloud(:, 1)), max(ptCloud(:, 1))];
y_range = [min(ptCloud(:, 2)), max(ptCloud(:, 2))];

x_len = x_range(2) - x_range(1);
y_len = y_range(2) - y_range(1);

% 计算所需的像素尺寸
% 根据点云的 X 和 Y 范围，动态计算二维投影图像的宽度和高度
% img_width = max(200, round((x_range(2)-x_range(1))*200)+100); % 分辨率
% img_height = max(200, round((y_range(2)-y_range(1))*200)+100);
% img_width = round((x_range(2)-x_range(1)) * 200) + 1; % 分辨率
% img_height = round((y_range(2)-y_range(1))* 200) + 1;
img_width = ceil(x_len / GridResolution) + 1;       % 向上取整以确保所有点都能包含
img_height = ceil(y_len / GridResolution) + 1;

% KDE
% 2. 使用 histcounts2 进行栅格计数
% histcounts2 性能优于 accumarray + sub2ind
% 创建用于 histcounts2 的 bin 边缘
x_edges = x_range(1):GridResolution:x_range(2) + GridResolution;
y_edges = y_range(1):GridResolution:y_range(2) + GridResolution;

% 统计，H 即为计数矩阵，其尺寸为 (length(y_edges)-1) x (length(x_edges)-1)
H = histcounts2(ptCloud(:, 2), ptCloud(:, 1), y_edges, x_edges);

figure;
imshow(H);
title('栅格计数点密度（直方）图');

%% 3. 应用高斯平滑来近似 KDE (离散“计数直方图”转换为连续“密度估计”)
% Bandwidth h (物理单位) = Sigma (像素单位) * GridResolution
Sigma_pixel =  0.001/ GridResolution; 

% H 矩阵的尺寸可能与 img_height/img_width 略有不同，但通常非常接近
% 由于 histcounts2 的输出通常比 ceil() 大的尺寸小 1，我们直接使用 H
img_kde = imgaussfilt(H, Sigma_pixel);

% 归一化 (可选，取决于可视化需求)
img_kde = img_kde / max(img_kde(:));

% 自适应阈值分割，确定高密度区域（通常对应植株）与低密度区域（通常对应地面）
binaryImg_KDE = imbinarize(img_kde, 'adaptive', 'Sensitivity', 0.5);

%% Visualize
figure;
% subplot(1, 2, 1);
% imagesc(img_kde);
% axis equal tight;
% title('高斯核平滑密度估计图');
% colorbar;
% 
% subplot(1, 2, 2);
imshow(binaryImg_KDE);
title('自适应阈值分割图');

%% ==== 多方向形态处理（关键步骤：弥补断裂叶片）====
% 超像素闭合： 有效地连接沿着线性结构元方向上的细小断裂和间隙，同时抑制与主要结构元方向不一致的噪声点和不规则突出。
angles = 0:15:179;   % 多方向结构元（可根据叶片方向性调整）
bw_clean = binaryImg_KDE;

for th = angles
    if th == 0
        se = strel('line', 15, th);         % 越小，小区域更容易被保留，但也更容易破碎；越大，会错误地连接和合并相邻区域，越容易形成一个大的连通区域，但分割精度会降低
    else
        se = strel('line', 5, th);
    end
%     marker = imerode(bw_clean, se);                 % 标记图像,腐蚀操作会切断细小的连接，只保留图像中比结构元 se 更大、更厚的部分。
%     bw_clean = imreconstruct(marker, bw_clean);     % 重构操作相当于用标记图像作为"种子"，在其上进行膨胀，但膨胀的范围被限制在掩模图像 bw_recon 的范围内。
    % 方向性形态学重建（如去除短线、保留长线），这正是 imreconstruct 的典型应用


% 多方向形态学梯度:
% 程序执行时间: 8.6855 秒
% N = 37
% 未被覆盖点数量：712
% 精度 95.7%
%     dilated = imdilate(bw_clean, se); % 前景区域“变胖”
%     eroded = imerode(bw_clean, se);   % “变瘦”
%     bw_clean = dilated - eroded;      % 形态学梯度 = 膨胀 − 腐蚀

    grad = imdilate(binaryImg_KDE, se) - imerode(binaryImg_KDE, se);
    bw_clean = bw_clean | grad;  % 或 max(edge_sum, grad)
end

figure;
imshow(bw_clean);
title('多方向闭运算');

%% 最终闭合（提供边界连贯性）
bw_close = imclose(bw_clean, strel('disk', 5));    % 越大，膨胀越明显，靠近的小区域会被合并，但是也更加容易与周边区域重叠
bw_close = imfill(bw_close,'holes');

figure;
imshow(~bw_close);
% title('闭合边界');
%% 基于距离或形状相似性的区域合并
% 开始计时
tic;

% 1: 识别所有区域和目标区域
% bw_close 是经过形态学处理后的二值图像
% AreaThreshold = 40000; % 设定区分大小区域的面积阈值（可调）
AreaThreshold = 90000; % 20241026
% AreaThreshold = 180000; % 20241214

% 获取所有连通区域的标记矩阵 L 和属性 stats
L = bwlabel(bw_close, 8); % 使用 bwlabel 确保连通性一致
stats = regionprops(L, 'Area', 'PixelIdxList'); % 用 PixelIdxList 更快读取点集

% 识别大区域和小区域的标签
large_labels = find([stats.Area] >= AreaThreshold);
small_labels = find([stats.Area] < AreaThreshold);

% 大区域掩模
BW_large = ismember(L, large_labels);

% 距离变换（到大区域的距离）
D = bwdist(BW_large);  % double matrix

% 初始结果：包含所有大区域（并保留原始大区域形状）
bw_final = BW_large;

% 预置邻域偏移（8邻域）
nb = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];

[rMax, cMax] = size(D);
maxStepsLimit = rMax*cMax; % 安全上限，防止无限循环

% 为 KD-tree 版（可选加速）：预取大区域像素坐标（用于可能查找目标 label）
[y_large, x_large] = find(BW_large); % 若需知道连到哪个大区域，可用此信息

% 逐个处理小区域
for idx = 1:length(small_labels)
    k = small_labels(idx);
    BW_small = (L == k);

    % 如果小区域已经与 BW_large 相交（意外情况），跳过
    if any(BW_small & BW_large, 'all')
        bw_final(BW_small) = 1;
        continue;
    end

    % 在小区域内找到距离 D 最小的位置（即离大区域最近的点）
    D_in_small = D;
    D_in_small(~BW_small) = Inf;
    [min_val, min_linidx] = min(D_in_small(:));

    if isinf(min_val)
        % 表示小区域完全没有到达大区域的路径（理论上不可能），跳过
        bw_final(BW_small) = 1;
        continue;
    end

    % 从该点出发，沿着 D 的局部最小方向移动，直到达到 D==0
    [r, c] = ind2sub(size(D), min_linidx);

    steps = 0;
    prev_rc = [-1 -1];
    while true
        steps = steps + 1;
        if steps > maxStepsLimit
            warning('Reached max steps for region %d; aborting path extension for this region.', k);
            break;
        end

        % 将当前位置设为前景（连接路径）
        bw_final(r,c) = 1;

        % 如果已经到达大区域边界（D == 0），结束
        if D(r,c) == 0
            break;
        end

        % 检查 8 邻域，找到邻域中 D 值最小的位置（贪心下降）
        minDnb = Inf;
        next_r = -1; next_c = -1;
        for t = 1:size(nb,1)
            rr = r + nb(t,1);
            cc = c + nb(t,2);
            if rr < 1 || rr > rMax || cc < 1 || cc > cMax
                continue;
            end
            val = D(rr,cc);
            if val < minDnb
                minDnb = val;
                next_r = rr; next_c = cc;
            end
        end

        % 如果找不到更小的邻域（被困在局部鞍点），尝试选择与当前点不同的邻域或跳出
        if next_r == -1
            % 尝试将相邻像素中任何一个设为前景并继续（容错策略）
            for t = 1:size(nb,1)
                rr = r + nb(t,1);
                cc = c + nb(t,2);
                if rr < 1 || rr > rMax || cc < 1 || cc > cMax
                    continue;
                end
                next_r = rr; next_c = cc;
                break;
            end
            if next_r == -1
                % 完全被困（四周越界），跳出
                break;
            end
        end

        % 防止原地循环（如果下一步是当前像素，则终止，避免死循环）
        if next_r == r && next_c == c
            break;
        end

        % 走向 next_r, next_c
        r = next_r;
        c = next_c;
    end

    % 最后把小区域本身也设为前景（确保没有漏掉）
    bw_final(BW_small) = 1;
end

% 结束计时并显示执行时间
execution_time = toc;
disp(['程序执行时间: ', num2str(execution_time), ' 秒']);

%% 可视化
figure;
imshow(bw_final);
title('基于距离或形状相似性的区域合并结果');

%% 查找二值图像（Binary Image）中对象的边界轮廓,提取连通域（用于后续凸包）
% B: 包含每个连通区域边界坐标的 cell 数组, B{k} 包含第 k 个区域的边界像素坐标 (行, 列)
% L: 标记矩阵（每个连通区域用一个唯一的非零整数标签）
% N: 连通区域的数量
[B, L, N] = bwboundaries(bw_final, 'noholes');      % 'noholes'：只查找外部边界，忽略被物体包围的内部孔洞（Holes）,加速计算。

%% 连通区域可视化
% 连通区域质心计算
stats = regionprops(L, 'Area', 'Centroid');

figure;
imshow(bw_final);
hold on;

% B 是一个 cell 数组，每个 cell B{k} 包含第 k 个区域的边界像素坐标 (行, 列)
for k = 1:N
   boundary = B{k};
   
   % 获取当前区域的面积和质心
   area = stats(k).Area;
   centroid = stats(k).Centroid; % [x, y] 坐标
   
   % 过滤掉不符合面积要求的区域（可选，如果之前 bwareaopen 未完全过滤）
   if area > 100 % 假设阈值为 100 像素
       % 绘制连通边界
       plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
       
       % 提取的 boundary 数组即为您需要的连通的、按序排列的边界像素！
       
       % 示例：计算轮廓周长（像素数）
       perimeter = size(boundary, 1); 
        
       % 显示标签
       text(centroid(1), centroid(2), [num2str(k) "=" num2str(area)], ...
            'Color', 'g', ...              % 标签颜色设为绿色
            'FontSize', 14, ...           % 字体大小
            'FontWeight', 'bold', ...     % 字体加粗
            'HorizontalAlignment', 'center'); % 文本居中对齐
       
       % 可以在此基础上进行其他形状特征计算，例如圆度、最小外接矩形等
   end
end
title('提取的连通边界轮廓');
hold off;

%% 可视化，检查凸包多边形是否在正确位置、大小与方向
figure(Position=[260, 320, 1500, 450]); 
scatter(ptCloud(:,1), ptCloud(:,2), 2, '.'); 
hold on;
for k = 1:length(B)
    boundary = B{k}; % 边界像素坐标转真实坐标
    
    x_coords = x_range(1) + (boundary(:,2) - 1) / (img_width - 1) * diff(x_range);
    y_coords = y_range(1) + (boundary(:,1) - 1) / (img_height - 1) * diff(y_range);

    if numel(x_coords) < 3 continue; end 
    
    ch_idx = convhull(x_coords, y_coords);

    plot(x_coords(ch_idx), y_coords(ch_idx), 'r-', 'LineWidth', 1.5); 

    % 显示区域序号
    cx2 = mean(x_coords);
    cy2 = mean(y_coords);
    text(cx2, cy2, sprintf('%d', k), 'Color','k','FontSize',18,'FontWeight','bold');
end
% axis equal;
title(sprintf('region %d', k));


%% 可视化：真实边界 和 alphashape
figure(Position=[260, 320, 1500, 450]); 
% scatter(ptCloud(:,1), ptCloud(:,2), 2, '.'); 
hold on;

alpha = 0.80; % 缩小比例 (可调)
colors = lines(length(B));  % 给每个区域不同颜色

for k = 1:length(B)
    boundary = B{k}; 
    if size(boundary,1) < 3
        continue;
    end

    % === 1. 像素边界 → 真实坐标 ===
    x_coords = x_range(1) + (boundary(:,2) - 1) / (img_width - 1) * (x_range(2) - x_range(1)); 
    y_coords = y_range(1) + (boundary(:,1) - 1) / (img_height - 1) * (y_range(2) - y_range(1)); 

    % === 2. 构造真实边界 polyshape ===
    p_real = polyshape(x_coords, y_coords);

    % === 6. 可视化所有图层 ===
    % 原始真实边界（细黑线）
    plot(x_coords, y_coords, 'Color', [0 0 0]+0.4, 'LineWidth', 0.8);

    % 凸包（红色） 
    plot(p_real, 'FaceColor', 'none', 'EdgeColor', [1 0 0 0.4], 'LineWidth', 1.2);

    % 显示区域序号
    cx2 = mean(x_coords);
    cy2 = mean(y_coords);
    text(cx2, cy2, sprintf('%d', k), 'Color','k','FontSize',18,'FontWeight','bold');
end

% title('D-plus 多边形交集法（真实边界 ∩ 缩小凸包）可视化');
% axis equal;
hold off;

%% 加载 ClusterCenters
% load([pathname 'maize-60-Cloud-SOR-plants222-optimized-Oct11-clusterCenters.mat']);
load([pathname '目标42.las.segmented.octree11-Cloud-Cloud-GT-clusterCenters.mat']);
% load([pathname '目标4-Cloud-ClassOctree10-GT-labels01-align-clusterCenters.mat']);
disp(size(clusterCenters));


%% 可视化
% ---------- 直接映射 ----------
clusterX_pixel = round((clusterCenters(:,1) - x_range(1)) / (x_range(2) - x_range(1)) * (img_width - 1)) + 1;
clusterY_pixel = round((clusterCenters(:,2) - y_range(1)) / (y_range(2) - y_range(1)) * (img_height - 1)) + 1;

clusterXY = [clusterX_pixel, clusterY_pixel];


figure;
imshow(~binaryImg_KDE);

hold on;
for i = 1:size(clusterCenters, 1)
    scatter(clusterXY(i, 1), clusterXY(i, 2), 50, 'r', 'x', 'LineWidth', 2); % 绘制聚类中心（红色X）
    text(clusterXY(i, 1), clusterXY(i, 2), int2str(i), 'Color', 'red', 'FontSize', 14);
end
hold off;


%% 凸包边界-old,用于保存遮挡集，单个处理进行调试
% targets = {}; % 存储目标点集 
% targets_info = []; % 存储目标的坐标和索引 
% hulls = {}; % 保存凸包坐标用于可视化
% 
% for k = 1:length(B)
%     boundary = B{k}; % 边界像素坐标转真实坐标
%     x_coords = x_range(1) + (boundary(:,2) - 1) / (img_width - 1) * (x_range(2) - x_range(1)); 
%     y_coords = y_range(1) + (boundary(:,1) - 1) / (img_height - 1) * (y_range(2) - y_range(1)); 
%     
%     if numel(x_coords) < 3 continue; end 
%     
%     K = convhull(x_coords, y_coords); % 找出点云中在凸包内的点 
%     inIdx = inpolygon(ptCloud(:,1), ptCloud(:,2), x_coords(K), y_coords(K)); 
%     target_points = ptCloud(inIdx, :); 
%     
%     if size(target_points, 1) < 10 
%         fprintf("%d : % d\n", k, size(target_points, 1)); 
%         continue; % 过滤小区域 
%     end
%     
%     % 每株单独计算边界框 bbox
%     min_x = min(target_points(:,1)); 
%     max_x = max(target_points(:,1)); 
%     min_y = min(target_points(:,2)); 
%     max_y = max(target_points(:,2)); 
%     min_z = min(target_points(:,3)); 
%     max_z = max(target_points(:,3)); 
%     
%     % 获取每个目标的 XY 坐标和索引 
%     xy_coords = mean(target_points(:, 1:2), 1); 
%     % 计算目标的平均 XY 坐标 
%     new_row = [xy_coords, k, min_x, max_x, min_y, max_y, min_z, max_z]; % 存储坐标和目标索引 
%     targets_info = [targets_info; new_row]; 
%     
%     targets{k} = target_points; 
%     
%     % 保存凸包坐标 
%     hulls{end+1} = [x_coords(K), y_coords(K)];
% end

%% 凸包边界
targets = {};       % 存储目标点集
targets_info = [];  % 存储目标的坐标和索引
hulls = {};         % 保存凸包坐标用于可视化

minPts = 10;        % 小于该值被认为是噪点集合
plantIdx = 0;       % 全局单株编号计数器

used_global_idx = false(size(ptCloud,1),1);

for k = 1:length(B)
    boundary = B{k};

    % 坐标映射
    x_norm = (boundary(:,2) - 1) / (img_width - 1);
    y_norm = (boundary(:,1) - 1) / (img_height - 1);

    x_coords = x_range(1) + x_norm * (x_range(2) - x_range(1));
    y_coords = y_range(1) + y_norm * (y_range(2) - y_range(1));

    if numel(x_coords) < 3
        continue;
    end

    ch_idx = convhull(x_coords, y_coords);
    ch_x = x_coords(ch_idx);
    ch_y = y_coords(ch_idx);

    % -----------------------------
    % Step 3：防止凸包重复覆盖
    % -----------------------------
    inIdx = inpolygon(ptCloud(:,1), ptCloud(:,2), ch_x, ch_y);

    % 去掉已经分配给其他 plant 的点
    inIdx = inIdx & ~used_global_idx;

    target_points = ptCloud(inIdx,:);

    if nnz(inIdx) < minPts
        continue;
    end

    % 标记为已使用
    used_global_idx(inIdx) = true;

    % === 根节点判断 & 分裂处理 ===
    hullXY = [ch_x, ch_y];

    inRoot = inpolygon(clusterCenters(:,1), clusterCenters(:,2), ch_x, ch_y);
    numRoots = sum(inRoot);
    rootIdx = find(inRoot);

    if numRoots > 1
        fprintf('Region %d detected occlusion (%d roots)\n', k, numRoots);
        [split_sets, split_hulls] = resolve_occlusion_subsets(target_points, clusterCenters(rootIdx,:));
    else
        split_sets = {target_points};
        split_hulls = {hullXY};
    end

    % === 存储分裂子集 ===
    for s = 1:length(split_sets)
        plantIdx = plantIdx + 1;

        currentPoints = split_sets{s};
        currentHull   = split_hulls{s};

        targets{plantIdx} = currentPoints;

        xy_coords = mean(currentPoints(:,1:2),1);

        bbox = [min(currentPoints(:,1)), max(currentPoints(:,1)), ...
                min(currentPoints(:,2)), max(currentPoints(:,2)), ...
                min(currentPoints(:,3)), max(currentPoints(:,3))];

        row = [xy_coords, plantIdx, bbox];
        targets_info = [targets_info; row];

        hulls{plantIdx} = currentHull;
    end
end

fprintf(">>> DONE. Assigned %d plants.\n", plantIdx);


%% 可视化凸包边界
figure;
imshow(~binaryImg_KDE); hold on;
for i = 1:numel(hulls)
%     disp(i);
    hull = hulls{i};
    hull_X = round((hull(:,1) - x_range(1)) / (x_range(2) - x_range(1)) * (img_width - 1)) + 1;
    hull_Y = round((hull(:,2) - y_range(1)) / (y_range(2) - y_range(1)) * (img_height - 1)) + 1;
    
    new_hull = [hull_X, hull_Y];
    
    plot(new_hull(:,1), new_hull(:,2), 'Color', 'green', 'LineWidth', 3);
end
title(sprintf('Detected Convex Hulls: %d', numel(hulls)));
hold off;

%% 确认 targets{k} 内是否含重复点
totalOriginal = size(ptCloud,1);
totalTargets = sum(cellfun(@(x) size(x,1), targets));

% 检查每一个 k 是否有重复点
for k = 1:length(targets)
    if isempty(targets{k}), continue; end
    [~, ia] = unique(targets{k}, 'rows');
    dup_k = size(targets{k},1) - length(ia);
    if dup_k > 0
        fprintf('第 %d 株存在 %d 个重复点。\n', k, dup_k);
    end
end

fprintf('总点数: %d\n', totalOriginal);
fprintf('targets 总点数: %d (多出 %d)\n', totalTargets, totalTargets-totalOriginal);



%% 保证 targets 内部无重复
for k = 1:length(targets)
    if ~isempty(targets{k})
        targets{k} = unique(targets{k}, 'rows');   % 清理重复点
    end
end



%% “未覆盖点”处理
% 找出所有未被 inpolygon 覆盖的点（被遗漏的点）
% 将每个未覆盖点分配给最近的植株（依据 hull 的最近距离或依据 root cluster center 最近距离）
% 插入到对应的 targets{k} 中

% === Step 1：统计所有被分配的点 ===
assigned = false(size(ptCloud,1),1);

for k = 1:length(targets)
    if isempty(targets{k}), continue; end
    pts = targets{k};
    
    % 找 ptCloud 对应的原始点索引（精确匹配）
    [~,loc] = ismember(pts, ptCloud, 'rows');
    assigned(loc(loc>0)) = true;
end

% === Step 2：未被任何凸包覆盖的点 ===
unassignedIdx = find(~assigned);

% 统计量输出
fprintf('已覆盖点数量：%d\n', length(find(assigned)));
fprintf('未被覆盖点数量：%d\n', length(unassignedIdx));


%% 建立每株的参考中心（用于最近匹配）
% Step 3：每株凸包的质心（比根节点更稳定）
numPlants = length(targets);
centroids = zeros(numPlants, 2);

for k = 1:numPlants
    if isempty(hulls{k}), continue; end
    centroids(k,:) = mean(hulls{k}, 1);    % hullXY 已为 N×2
end

% Step 4：将未覆盖点分配给最近的凸包中心
% 所有未覆盖点的 XY
P = ptCloud(unassignedIdx, 1:2);    % M × 2

% 计算每个点到所有质心的距离:  M×N 距离矩阵
D = pdist2(P, centroids);  

% 为每个点找到最近的植株索引
[~, nearestPlant] = min(D, [], 2);   % M×1

% 按植株批量加入
for k = 1:numPlants
    idx_k = unassignedIdx(nearestPlant == k);  % 属于第 k 株的未覆盖点
    if ~isempty(idx_k)
        % 批量加入，避免重复出现
        targets{k} = [targets{k}; ptCloud(idx_k,:)];
    end
end

finalPointsNum = sum(cellfun(@length, targets));
fprintf('最终点数: %d，与真值相差: %d\n, ', finalPointsNum, targetPoints.Count - finalPointsNum)


%% === 调用排序函数对植株进行分组排序 ===
[groupSorted_targets, groupNums, others] = seg_box(targets_info, groupSize);

%% 生成独立植株点云并进行三维可视化
colors = lines(numel(hulls));

% 创建三维图形
figure;
scatter3(ptCloud(:, 1), ptCloud(:, 2), ptCloud(:, 3), 10, [0.5, 0.5, 0.5], 'filled');
hold on;

% 预分配（可选，提升效率）
targetsNum = size(groupSorted_targets, 1);
allPointsCell = cell(1, targetsNum);
allLabelsCell = cell(1, targetsNum);

for i = 1:targetsNum
    targetIdx = groupSorted_targets(i, 3);

    current_targets = targets{targetIdx};

    try
        K3 = convhull(current_targets(:,1), current_targets(:,2), current_targets(:,3));
        trisurf(K3, current_targets(:,1), current_targets(:,2), current_targets(:,3), ...
            'FaceColor', colors(targetIdx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    catch
        warning('Convex hull failed for cluster %d', targetIdx);
    end

    % 显示标签
    % 计算标签显示位置
    avg_x = min(current_targets(:, 1)) + 0.1;
    avg_y = mean(current_targets(:, 2)) + 0.2;
    avg_z = max(current_targets(:, 3)) + 0.2;

    % 生成标签信息
    text(avg_x, avg_y, avg_z, sprintf('%02d', i), ...
         'Color', 'black', 'FontSize', 14, 'FontWeight', 'normal');

    fprintf('Region %d have %d points.\n', i, size(current_targets, 1));
    
    allPointsCell{i} = current_targets;
    allLabelsCell{i}  = i * ones(size(current_targets, 1), 1);
end

% title([num2str(length(targets)), ' Point Cloud Objects with Bounding Boxes']);
xlabel('X');
ylabel('Y');
zlabel('Z');
% axis equal;
grid on;
hold off;

%% 保存：将每个点的 Intensity 属性保存为相应类别
% 合并所有点和标签
allPoints = vertcat(allPointsCell{:});
allLabels = vertcat(allLabelsCell{:});
% 构建新的点云对象
ptCloud_new = pointCloud(allPoints, 'Intensity', allLabels);

% 写入PLY
[~, name, ~] = fileparts(filename);
pcwrite(ptCloud_new, [pathname name '_wiithSingleLabeles4.ply'], 'PLYFormat', 'binary');

%% 遮挡集判断
regionType = cell(targetsNum,1);

for i = 1:targetsNum
    regionXY = allPointsCell{i};
    if size(regionXY,1) < 3
        regionType{i} = 'noise';
        continue;
    end
    k = convhull(regionXY(:,1), regionXY(:,2));
    hullXY = regionXY(k,:);
    
    % 判断 clusterXY 中哪些点在该区域内
    [in,~] = inpolygon(clusterCenters(:,1), clusterCenters(:,2), hullXY(:,1), hullXY(:,2));
    numRootsInRegion = sum(in);
    
    % 分类
    if numRootsInRegion == 1
        regionType{i} = 'single';
    elseif numRootsInRegion >= 2
        regionType{i} = 'occlusion';
    else
        regionType{i} = 'noise';
    end
end

%% 视觉验证
figure; hold on;
plot(clusterCenters(:,1), clusterCenters(:,2), 'ko', 'MarkerFaceColor','k');
ocNum = 1;
for i = 1:targetsNum
    regionXY = allPointsCell{i};
    % 能把所有点"包"在里面的最小凸多边形
    % k 是一个 索引向量，按逆时针顺序列出构成凸包的那些点在 regionXY 中的行下标
    k = convhull(regionXY(:,1), regionXY(:,2));
    hullXY = regionXY(k,:);
    
    mean_XY = mean(regionXY(:,1:2));

    if strcmp(regionType{i}, 'occlusion')
        plot(hullXY(:,1), hullXY(:,2), 'r-', 'LineWidth',1.5);
        text(mean_XY(1), mean_XY(2), ['Occlucion-' int2str(ocNum)], 'Color', 'black', 'FontSize', 18, 'FontWeight', 'normal');
        ocNum = ocNum + 1;
    elseif strcmp(regionType{i}, 'single')
        plot(hullXY(:,1), hullXY(:,2), 'g-', 'LineWidth',1);
    else
        plot(hullXY(:,1), hullXY(:,2), 'b--');
    end
end
title('单株与遮挡集判定结果');
legend('root centers','occlusion','single','noise');


%% 遮挡集数据保存
[~, name, ~] = fileparts(filename);
mat_file = [pathname name '-occlusion_plants.mat'];
save(mat_file, 'regionType', 'allPointsCell');



%% 辅助函数：Bresenham 直线算法
function line = bresenham(x1, y1, x2, y2)
    dx = abs(x2-x1);
    dy = abs(y2-y1);
    steep = dy > dx;
    if steep
        [x1, y1] = deal(y1, x1);
        [x2, y2] = deal(y2, x2);
        dx = abs(x2-x1);
        dy = abs(y2-y1);
    end
    if x1 > x2
        [x1, x2] = deal(x2, x1);
        [y1, y2] = deal(y2, y1);
    end
    derr = 2 * dy;
    err = 0;
    y = y1;
    ystep = 1;
    if y1 > y2
        ystep = -1;
    end
    line = [];
    for x = x1 : x2
        if steep
            line = [line; y, x];
        else
            line = [line; x, y];
        end
        err = err + derr;
        if err > dx
            y = y + ystep;
            err = err - 2*dx;
        end
    end
end

