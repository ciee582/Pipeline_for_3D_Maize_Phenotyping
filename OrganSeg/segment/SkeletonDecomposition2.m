function [ organ_subskeleton, skeletonType, joints] = SkeletonDecomposition2(mainskeleton, spls, joints, roots, spls_adj, corresp, features, show_results )
%CLASSIFYROOT  recognize the only stem root vertex from all the root vertices
%output: organ_subskeleton .The last element is the stem sub_skeleton. The rest is the leaf sub_skeleton

%%%%%%%%%%%%% 初始化邻接矩阵 %%%%%%%%%%%%%
s_adj_Num = size(spls, 1);
adjmatrix = zeros(s_adj_Num, s_adj_Num);

% 非连接点权值设为 inf，相邻的骨架点权重为 1
for i = 1:s_adj_Num
    for j = 1:s_adj_Num
        temp = spls_adj(i,j);
        if i == j
            adjmatrix(i,j) = 0;
        elseif temp == 0
            adjmatrix(i,j) = inf;
        else
            adjmatrix(i,j) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%% 提取初步子骨架 %%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------
% 对每个端点 roots，去寻找其通向所有节点 joints 的最短路径
% 取最短路径中的最短者作为叶片骨架片段（对应某个叶片）
% ------------------------------------------------
rootNum = size(roots,1);
jointNum = size(joints,1);

organ_subskeleton = cell(0,1);
sub_index = 1;

for i = 1:rootNum
    r_index = roots(i);
    if ~ismember(mainskeleton, r_index)     % 对不在主骨架中的点不进行处理，因为有可能是孤立点
        continue;
    end

    Dist = inf(jointNum, 1);        % 子骨架距离数组
    Path = cell(jointNum, 1);       % 子骨架路径数组

    for j = 1:jointNum
        c_index = joints(j);

        [t1, t2] = mydijkstra(adjmatrix, r_index, c_index);     % 端点到节点的最短距离 t1 和路径 t2

        if ~isempty(t2)
            Dist(j) = t1;
            Path{j} = t2;
        end
    end

    [~, minDist_idx] = min(Dist);   % 提取最短距离
    minPath = Path{minDist_idx};    % 提取最短路径
    if isempty(minPath)
        continue; % 跳过无法到达的路径
    end

    organ_subskeleton{sub_index} = minPath(:);      % 添加到子骨架元组中
    sub_index = sub_index + 1;
end

% visualize_subskeleton(spls, organ_subskeleton);

%%%%%%%%%%%%% 找出唯一非叶片（茎秆）子骨架 %%%%%%%%%%%%%
% ==== 计算每个子骨架的方向向量、平均 feature 和长度特征 ====
num_organ_sub = length(organ_subskeleton);
skefeature = [];
skelen = [];
d_v = zeros(num_organ_sub, 3);
Z_min_list = zeros(num_organ_sub,1);

for i = 1:num_organ_sub
    indices = organ_subskeleton{i};

    if length(indices) < 1      % 子骨架点数过小
        skefeature = [skefeature; 0];
        d_v(i,:) = [0, 0, 0];
        skelen = [skelen; 0];
        Z_min_list(i) = inf;
        continue;
    end

    Z_min_list(i) = min(spls(indices, 3));

    skelen = [skelen; length(indices) / length(spls)];      % 子骨架点数占骨架整体的比例

    % 平均 feature（例如：点云密度、粗细等）
    pids = ismember(corresp, indices);          % 提取子骨架中所有点在原始点云中对应的数据点集合
    fsum = sum(features(pids));                 % 计算原始点的特征之和
    psum = sum(pids);                           % 子骨架对应原始点云的点数
    mf = fsum / max(psum, 1);  % 避免除0
    skefeature = [skefeature; mf];

    vi_end = indices(end);
    vi_mid = indices(floor(end / 2));
    v_end = spls(vi_end, :);
    v_mid = spls(vi_mid, :);

    dir_vec = v_mid - v_end;        % 中点 和 末端点组成线段的方向向量
    if norm(dir_vec) == 0           % 模长为0，即距离为 0
        d_v(i,:) = [0, 0, 0];
    else
        d_v(i,:) = dir_vec / norm(dir_vec);     % 单位方向向量
    end
end

%点积是一个向量与另一个向量的乘积，其结果是一个标量，表示两个向量之间的相似程度。
% 如果两个向量的点积为正，那么它们的方向相似；
% 如果为负，那么它们的方向相反；
% 如果为零，那么它们垂直。
% 茎秆：点积之和最小（或接近0）的向量，因为茎秆与大多数叶片的方向夹角较大。
% 叶片：点积之和较大的向量，因为叶片之间有较多的正点积。
% d_v2 = d_v * d_v';
% sumV = sum(d_v2, 1);
% ids = find(sumV <= 0);
%
% if length(ids) == 1
%     stem_root =  organ_subskeleton{ids}(1);
%     organ_subskeleton(ids) = [];        % 删除识别到的茎秆子骨架，避免重复计入
% elseif isempty(ids)
%     error('Unable to determine stem root.');
% else
%     % skelen = skelen + skefeature;
%     skelen = skefeature;
%     lens = skelen(ids);
%     % lens = skelen;
%     [~, id] = max(lens);
%     index = ids(id);
% 
%     stem_root = organ_subskeleton{index}(1);
% end


% ==== 进行 stem 判定 ====
% % 主茎方向相似度（与Z轴）
% z_axis_candidates = [0 0 1; 0 0 -1];
% z_scores = abs(d_v * z_axis_candidates');  % Nx2
% dir_similarity = max(z_scores, [], 2);  % 每行取 max（越接近竖直越大）

% ---------- Step 1: 基于连接点计算茎秆主方向 ---------
% 获取连接点
connPts = spls(joints, :);

% 若连接点数量足够，使用 PCA 拟合主方向
if size(connPts,1) >= 3
    [coeff, ~, ~] = pca(connPts);
    stem_axis = coeff(:,1)';  % 第一主成分方向（主茎方向估计）
%     ptMaxZ = connPts(connPts(:,3) == max(connPts(:,3)), :);  % Z最大的点（取第一个）
%     ptMinZ = connPts(connPts(:,3) == min(connPts(:,3)), :);  % Z最小的点（取第一个）
%     stem_axis = ptMaxZ - ptMinZ;
else
    stem_axis = [0 0 1];  % 回退策略：默认竖直方向
end

% ---------- Step 2: 计算每个子骨架方向与 stem_axis 的相似度 ----------
axis_similarity = abs(d_v * stem_axis');  % Nx1，越大越接近茎秆方向

% ==== 综合得分（加权可调）====
% 得分越低越可能是茎秆
score = ...
    0.3 * normalize(skefeature) + ...
    0.2 * normalize(skelen) + ...
    0.5 * normalize(1 - axis_similarity);

% 计算方向向量点积矩阵
dotM = d_v * d_v';
sumV = sum(dotM, 1);
geom_score = sumV ./ num_organ_sub;   % 平均方向一致性（越低越可能为垂直主轴）

% 子骨架Z权重，越靠近最低处（Z小） → 权重越大
Z_min_all = min(Z_min_list(~isinf(Z_min_list)));
Z_max_all = max(Z_min_list(~isinf(Z_min_list)));

% 归一化Z权重（越低越接近1）
Z_weight = (Z_max_all - Z_min_list) ./ (Z_max_all - Z_min_all + 0.001);
Z_weight(isinf(Z_min_list)) = 0;  % 异常子骨架不参与评分


% 综合更新评分
final_score = 0.3*score + 0.7*normalize(geom_score(:)) - Z_weight;

mask = Z_weight > 0.7;          % 只考虑将高度低于 0.7 的子骨架选做茎秆
final_score(~mask) = inf;

[~, stem_idx] = min(final_score);
% ---------- Step 5: 删除茎秆子骨架，避免重复 ----------
if stem_idx <= numel(organ_subskeleton)
    stem_points = organ_subskeleton{stem_idx};
    if isempty(stem_points)
        error('Stem sub-skeleton is empty. Unable to determine stem root.');
    end
else
    error('检测的 stem_idx 超出范围。');
end

% ---------- 自动根节点提取 ----------
% stem_points 局部值作为根节点
if numel(stem_points) >= 3
    pts = spls(stem_points, :);

    % 计算逐段方向向量并单位化
    dirs = diff(pts, 1, 1);
    dirs = dirs ./ vecnorm(dirs, 2, 2);

    % 计算相邻方向夹角（度）
    angles = acosd(sum(dirs(1:end-1,:) .* dirs(2:end,:), 2));

    % 找最大方向突变点（即潜在L型拐点）
    [max_angle, idx_turn] = max(angles);
    angle_thresh = 60;  % 角度阈值，可根据骨架平滑程度调整

    % 检测L型骨架：若整体方向变化明显则更新 local_min_index
    global_angle = acosd(dot(dirs(1,:), dirs(end,:)));  % 整体首尾夹角
    if max_angle > angle_thresh || global_angle > 60
        local_min_index = stem_points(idx_turn + 1);
        fprintf('[L-shape detected] 拐点角度 = %.2f°，更新 local_min_index = %d。\n', ...
                max_angle, local_min_index);
    else
        fprintf('[Normal shape] 未检测到明显拐点（max_angle = %.2f°），使用局部最小Z值点作为根节点。\n', max_angle);
        [~, local_min_idx] = min(spls(stem_points,3));
        local_min_index = stem_points(local_min_idx);
    end
else
    [~, local_min_idx] = min(spls(stem_points,3));
    local_min_index = stem_points(local_min_idx);
end

% 若局部最小点比全局Z略高不多，则使用局部点
% if spls(local_min_index,3) < spls(all_min_idx, 3) + 0.01
%     stem_root = all_min_idx;
% elseif spls(local_min_index,3) < spls(min_joints, 3) + 0.001
%     stem_root = local_min_index;
%     organ_subskeleton(stem_idx) = [];   % 删除茎秆子骨架，避免重复
% else
%     stem_root = min_joints;
% end

% ==== 从确定的主茎起点出发，找到其与 joints 连接的所有最短路径 ====
% ---------- 计算茎秆路径 ----------
Dist = [];
Path = [];

if ismember(local_min_index, joints)
    % 情况1：根节点在 joints 中 —— 按原逻辑计算所有路径
    for j = 1:jointNum
        c_index = joints(j);
        if ~ismember(mainskeleton, c_index)
            continue;
        end
        [t1, t2] = mydijkstra(adjmatrix, local_min_index, c_index);
        Dist(j) = t1;
        Path{j} = t2;
    end

else
    % 情况2：根节点不在 joints 中
    Dist_local = [];
    Path_local = [];

    % Step 1:  找到 joints 最低点到各 joints 的路径    
    % 连接点最小点
    [~, conn_min_idx] = min(spls(joints, 3));
    min_joints = joints(conn_min_idx);
%     [~, all_min_idx] = min(spls(:, 3));
    
    for j = 1:jointNum
        c_index = joints(j);
        if ~ismember(mainskeleton, c_index)
            continue;
        end
        [d_tmp, p_tmp] = mydijkstra(adjmatrix, min_joints, c_index);
        Dist_local(j) = d_tmp;
        Path_local{j} = p_tmp;
    end

    % Step 2: 计算 min_joints → local_min_index 的路径
    [d_root, p_root] = mydijkstra(adjmatrix, local_min_index, min_joints);

    % Step 3: 拼接路径（每个 joints 的路径加上根段）
    for j = 1:jointNum
        if ~ismember(mainskeleton, joints(j))
            continue;
        end
        if isempty(Path_local{j})
            Dist(j) = NaN;
            Path{j} = [];
            continue;
        end
        Dist(j) = Dist_local(j) + d_root;
        Path{j} = [p_root, Path_local{j}(2:end)];  % 从2开始，避免 local_min_index 重复
    end
    
    % Step 4: 更新连接点
    joints(end+1) = local_min_index;

    % Step 5: 更新子器官子集,将原有叶端点至 min_joints 的路径修改为至 local_min_index
    local_idx = find(organ_subskeleton{stem_idx} == local_min_index);
    organ_subskeleton{stem_idx} = organ_subskeleton{stem_idx}(1:local_idx);
end

% ==== 获取相应路径作为茎骨架 ==== 
[~, stem_idx] = max(Dist); % 直接取最大路径作为主骨架 
% [~, sorted_idx] = sort(Dist, 'descend');  % 路径按长度排序（候选列表）
% stem_idx = sorted_idx(2);     % 选择次长路径作为茎秆
pathNodes  = Path{stem_idx};
target_path_dist = Dist(stem_idx) * 0.8;    % 以整个候选路径的 60% 作为路径，避免植株顶部聚集导致的叶片误划分为茎秆

% joints 在该 path 上的交集
pathJoints = intersect(pathNodes , joints, 'stable');

cumDist = 0;
chosenJoint = pathNodes(end); % 默认最后一个 joint，即选择全部路径作为茎秆节点

% 计算主路径上节点间的距离，按照 target_path_dist 进行筛选
for i = 2:length(pathNodes)
    % 当前边长
%     d = norm(spls(pathNodes(i),:) - spls(pathNodes(i-1),:));
    [d, ~] = mydijkstra(adjmatrix, pathNodes(i), pathNodes(i-1));
    cumDist = cumDist + d;

    % 如果是 joint 节点并且累积距离超过目标
    if ismember(pathNodes(i), pathJoints) && cumDist >=  target_path_dist
        chosenJoint = pathNodes(i);
        break;
    end
end

% 茎秆路径 = 从根节点到 chosenJoint
stem_path = pathNodes(1:find(pathNodes==chosenJoint));

% numCandidates = min(5, numel(sorted_idx));  % 最多取前5条长路径来评估
% scores = zeros(numCandidates,1);
% 
% for k = 1:numCandidates
%     idx = sorted_idx(k);
%     path = Path{idx};
% 
%     % ---- 特征计算 ----
%     % 1) 路径长度（已知 Dist）
%     L = Dist(idx);
% 
%     % 2) 路径直度 (端点直线距离 / 实际路径长度)
%     endpts = [spls(path(1), :); spls(path(end), :)];  % 假设 path 是 Nx3 点坐标
%     euclid_len = norm(endpts(2,:) - endpts(1,:));
%     straightness = euclid_len / L;
% 
%     % 3) 与Z轴夹角 (近似用端点连线方向)
%     dir_vec = (endpts(2,:) - endpts(1,:)) / euclid_len;
%     angleZ = acos(abs(dot(dir_vec, [0 0 1]))); % 与竖直的夹角，越小越好
%     vert_score = cos(angleZ);  % cos值越接近1越竖直
% 
%     % 4) 与植株主轴对齐 (用整株点云或所有路径做PCA)
%     % 假设已经有 plantAxis （1x3 向量，植株整体主方向）
%     if exist('plantAxis','var')
%         anglePCA = acos(abs(dot(dir_vec, plantAxis)));
%         pca_score = cos(anglePCA);
%     else
%         pca_score = vert_score; % 没算PCA就先用竖直性代替
%     end
% 
%     % ---- 综合打分 ----
%     % 权重可调: 长度0.25, 直度0.2, 主轴对齐0.3, 垂直0.25
%     scores(k) = 0.10*(L/max(Dist)) + ...
%                 0.30*straightness + ...
%                 0.30*pca_score + ...
%                 0.30*vert_score;
% end
% 
% % 选择得分最高的候选作为茎
% [~, bestIdx] = max(scores);
% stem_idx = sorted_idx(bestIdx);
% stem_path = Path{stem_idx};

if isempty(stem_path)
    error('无法从 stem_root 构建主茎骨架。');
end

% % 主茎路径平滑优化
% for i = 3:length(stem_path)
%     index1 = stem_path(i-2);
%     index2 = stem_path(i-1);
%     index3 = stem_path(i);
%     spl_1 = spls(index1,:);
%     spl_2 = spls(index2,:);
%     spl_3 = spls(index3,:);
%     r = norm(spl_3 - spl_2);
%     dir = (spl_2 - spl_1) / norm(spl_2 - spl_1);
%     spl_4 = spl_2 + dir * r;
%     spls(index3,:) = 0.5 * (spl_4 + spl_3);
% end

% === 关键：剔除原器官中与茎秆重叠的骨架点 ===
isSubset = cellfun(@(x) all(ismember(x, stem_path)), organ_subskeleton);
organ_subskeleton = organ_subskeleton(~isSubset);
organ_subskeleton{end+1} = stem_path(:);

% 输出类型标记（1=叶片，2=主茎）
skeletonType = [ones(1, length(organ_subskeleton))];
skeletonType(end) = 2;

%% 叶片按从下到上排序，主茎放最后一个
numOrgans = numel(organ_subskeleton);
numLeaves = numOrgans - 1;

leaf_indices = 1:numLeaves;
stem_index = numOrgans;

leaf_base_z = zeros(numLeaves, 1);  % 每片叶的连接点Z坐标
leaf_length = zeros(numLeaves, 1);  % 每片叶骨架长度

for i = 1:numLeaves
    leaf_path = organ_subskeleton{i};
    leaf_points = spls(leaf_path, :);

    % 连接点为最后一个骨架点
    leaf_base = leaf_points(end, :);
    leaf_base_z(i) = leaf_base(3);

    % 计算叶片骨架长度（累积欧氏距离）
    if size(leaf_points,1) > 1
        % leaf_length(i) = sum(vecnorm(diff(leaf_points), 2, 2));
        leaf_length(i) = numel(leaf_path);
    else
        leaf_length(i) = 0; % 防止只有一个点的异常情况
    end

%     fprintf('leaf-%d(%d~%d): %.4f %2d\n', i, leaf_path(1), leaf_path(end), leaf_base_z(i), leaf_length(i));
end

%% 多级排序规则：
% 1）按连接点Z坐标降序（高→低）
% 2）若Z相近，则按叶片长度升序（短→长）

% 按Z坐标初步排序
% ascend 正序，从小到达； descend 逆序，从大到小
[~, idx1] = sort(leaf_base_z, 'descend');

% 针对Z相同的情况，进一步按长度排序
sorted_idx = idx1;
for i = 1:length(idx1)
    sameZ = abs(leaf_base_z(idx1) - leaf_base_z(idx1(i))) < 1e-6;
    if sum(sameZ) > 1
        sub_idx = idx1(sameZ);
        if leaf_base_z(idx1(i)) > max(leaf_base_z)*0.8
            [~, sub_order] = sort(leaf_length(sub_idx), 'ascend');      % 顶部叶片长的生长较早
        else
            [~, sub_order] = sort(leaf_length(sub_idx), 'descend');     % 底部叶片短的生长较早
        end
        sorted_idx(sameZ) = sub_idx(sub_order);
    end
end

sorted_leaf_indices = leaf_indices(sorted_idx);

% 最终排序：叶片在前依次排列，主茎最后
sorted_indices = [sorted_leaf_indices, stem_index];

% 重排
organ_subskeleton = organ_subskeleton(sorted_indices);
skeletonType = skeletonType(sorted_indices);

if show_results
    disp(['识别到 ', num2str(length(organ_subskeleton) - 1), ' 片叶子和 1 个主茎']);
    visualize_subskeleton(spls, organ_subskeleton, joints, roots);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visualize_subskeleton(spls, organ_skeleton, joints, roots)
figure('Name','SkeletonDecomposetion','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
hold on;

LineColors = MyGS.MYCOLOR;
organ_num = length(organ_skeleton);
for idx = 1:organ_num
    path = organ_skeleton{idx};

    if idx == organ_num  % 主茎加粗
        plot3(spls(path, 1), spls(path, 2), spls(path, 3), '-', 'Color', 'green', 'LineWidth', 3);
        scatter3(spls(path(1), 1), spls(path(1), 2), spls(path(1), 3), 100, 'go', 'filled');  % 起点
    else
        plot3(spls(path, 1), spls(path, 2), spls(path, 3), '-', 'Color', LineColors(idx, :), 'LineWidth', 1);
    end

    text(spls(path(1), 1), spls(path(1), 2), spls(path(1), 3), int2str(organ_num + 1 - idx), 'Color', 'red', 'FontSize', 18);

%     pause(2);
end

% scatter3(spls(:,1),spls(:,2),spls(:,3), 2, [0,0,0], 'filled');
% for i = 1:length(spls)
%     if ismember(i, joints) || ismember(i, roots)
%         text(spls(i,1),spls(i,2),spls(i,3), num2str(i), 'Color', 'black', 'FontSize', 18);
%     end
% end

hold off;
axis off; axis equal; axis vis3d; view3d rot;
end


function x_norm = normalize(x)
if max(x) == min(x)
    x_norm = zeros(size(x));
else
    x_norm = (x - min(x)) / (max(x) - min(x));
end
end

