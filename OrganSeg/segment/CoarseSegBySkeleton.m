function [ Phi_O, Phi_U, Phi_S] = CoarseSegBySkeleton( points, spls, corresp, sub_skeletons, joints, K, show_results )
% 对输入的点云（points）进行基于骨架的粗分割
% 利用骨架点（sub_skeletons）、点云到骨架点的对应关系（corresp）以及关节点（joints），将点云分为三个部分：
% 主要骨架点集（Phi_S）：与主要骨架（sub_skeletons{end}）对应的点云点。
% 次要骨架点集（Phi_O）：与非主要骨架（sub_skeletons{1:end-1}）对应的点云点，分为多个子集。
%   未分配点集（Phi_U）：不属于任何骨架的点云点。

% 初始化
Phi_O = [];
AllPhi = [];

% 处理次要骨架（Phi_O）
for i = 1:length(sub_skeletons)-1
    indices = sub_skeletons{i};
    Phi_O{i} = [];
    for j = 1:length(indices)-1         % 不计算 indices{end} 对应的 joints 的关联点
        indices1 = find(corresp == indices(j));
        Phi_O{i} = [Phi_O{i}; indices1];
        AllPhi = [AllPhi; indices1];
    end
end

% 计算未分配点集（Phi_U）
AllIndices = 1:length(points);
Phi_U = setdiff(AllIndices, AllPhi)';

% 处理茎秆骨架（Phi_S）
indices = sub_skeletons{end};  % stem skeleton
numPts = length(indices);

% 只标记茎秆骨架上的点是否已被处理（防止同一个 skelIdx 被重复映射）
assignedPts = false(numPts, 1);

Phi_S = [];

% 参数设置
neighborhood_half_width = 3;    % 上下各 4 个点
r_anchor = 0.0115;               % 小半径锚定

for j = 1:numPts
    if assignedPts(j)
        continue;
    end
    
    skelIdx = indices(j);
    
    % 进行局部半径限制映射
    indices1 = mapJointNeighborhoodLimited(j, numPts, indices, corresp, points, spls, neighborhood_half_width, r_anchor);
    
    % 严格过滤非茎秆点
    if ~isempty(indices1)
        indices1 = indices1(ismember(corresp(indices1), indices));
    end
    
    % 标记整个邻域为已处理
    start_pos = max(1, j - neighborhood_half_width);
    end_pos   = min(numPts, j + neighborhood_half_width);
    assignedPts(start_pos:end_pos) = true;
    
    if ~isempty(indices1)
        Phi_S = [Phi_S; indices1];
    end
    
    assignedPts(j) = true;  % 普通点也标记（冗余但安全）
end

if show_results
    figure('Name', 'CoarseSegBySkeleton', 'NumberTitle', 'off');set(gcf,'color','white');movegui('southwest');
    hold on;
    % 先画未分配点
    scatter3(points(Phi_U,1),points(Phi_U,2),points(Phi_U,3), 15, [0 0 0], 'filled');
    color = MyGS.MYCOLOR;
    % 预留图例句柄
    h_leaf = scatter3(NaN,NaN,NaN, 3, [0 1 0], 'filled');
    h_stem = scatter3(NaN,NaN,NaN, 3, color(end, :), 'filled');
    h_unassigned = scatter3(NaN,NaN,NaN, 3, [0 0 0], 'filled');
    for i=1:length(Phi_O)
        indices = Phi_O{i};
        scatter3(points(indices,1),points(indices,2),points(indices,3), 5, color(i,:), 'filled');
    end
    scatter3(points(Phi_S,1),points(Phi_S,2),points(Phi_S,3), 5, color(end, :), 'filled');
    % 图例
    legend([h_leaf h_stem h_unassigned], {'叶片(绿色)','茎秆','未分配'}, 'Location','best');
    hold off;
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view3d rot;
end

end

function ptsIdx = mapJointNeighborhoodLimited(center_pos, numPts, indices, ...
                                              corresp, points, SkelPts, ...
                                              half_width, r_cylinder)
    % MAPJOINTNEIGHBORHOODLIMITED
    %   对茎秆骨架上以 center_pos 为中心的 ±half_width 邻域骨架点，
    %   进行局部半径 r 的点云映射（合并去重）
    %
    % 输入：
    %   center_pos   - 当前 joint 在 indices 中的位置（标量，1~numPts）
    %   numPts       - length(indices)
    %   indices      - 茎秆骨架所有 skeleton 索引向量 (1 x numPts)
    %   corresp, points, SkelPts - 原有映射参数
    %   half_width   - 邻域半宽（例如 4 表示上下各取 4 个点，共 9 个）
    %   r_cylinder   - 局部映射圆柱半径阈值
    %
    % 输出：
    %   ptsIdx       - 满足条件的点云索引（列向量，去重后）

    % Step 1: 获取邻域骨架点位置和坐标
    start_pos = max(1, center_pos - half_width);
    end_pos   = min(numPts, center_pos + half_width);
    neighbor_pos = start_pos : end_pos;
    neighbor_skelIdxs = indices(neighbor_pos);
    local_skel_pts = SkelPts(neighbor_skelIdxs, :);  % [M x 3]

    if size(local_skel_pts, 1) < 2
        % 邻域太小，无法拟合直线，退化为球形距离（安全回退）
        center_pt = local_skel_pts(1, :);
        cand = find(corresp == neighbor_skelIdxs(1));
        if isempty(cand)
            ptsIdx = [];
            return;
        end
        dist = vecnorm(points(cand, :) - center_pt, 2, 2);
        ptsIdx = cand(dist <= r_cylinder);
        return;
    end

    % Step 2: 拟合局部直线（使用 PCA 主方向）
    % 方法：PCA 获取主方向向量 v，选取端点之一作为基点 p0
    centroid = mean(local_skel_pts, 1);
    [~, ~, V] = svd(local_skel_pts - centroid, 'econ');  % V 的第一列是主方向
    direction = V(:, 1)';  % [1 x 3] 主方向向量（单位化）
    direction = direction / norm(direction);

    % 选取邻域中一个端点作为直线上的参考点（更好数值稳定性）
    p0 = local_skel_pts(1, :);  % 或用 centroid 也行

    % Step 3: 收集所有邻域骨架点对应的候选点云点
    all_cand = [];
    for i = 1:length(neighbor_skelIdxs)
        skelIdx = neighbor_skelIdxs(i);
        cand = find(corresp == skelIdx);
        if ~isempty(cand)
            all_cand = [all_cand; cand];
        end
    end
    if isempty(all_cand)
        ptsIdx = [];
        return;
    end
    candidate_pts = points(all_cand, :);  % [K x 3]

    % Step 4: 计算每个候选点到直线的距离
    % 向量：p0 -> pi
    vectors = candidate_pts - p0;  % [K x 3]
    
    % 投影长度（标量）
    proj_length = vectors * direction';  % [K x 1]
    
    % 最近点坐标 = p0 + proj_length * direction
    closest_pts = p0 + proj_length .* direction;  % [K x 3]
    
    % 点到直线的欧氏距离
    dist_to_line = vecnorm(candidate_pts - closest_pts, 2, 2);  % [K x 1]

    % Step 5: 保留圆柱半径内的点
    valid_mask = dist_to_line <= r_cylinder;
    ptsIdx = all_cand(valid_mask);

    % 最终去重（虽然通常不会重复）
    ptsIdx = unique(ptsIdx);
end

