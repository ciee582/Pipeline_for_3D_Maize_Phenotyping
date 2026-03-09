function [plant_sets, plant_hulls] = resolve_occlusion_subsets(points, roots)
    %% 1. 沿Z轴分层处理
    % 获取Z轴范围
    zCoords = points(:, 3);
    zMin = min(zCoords);
%     layerHeight = 0.04;     % 每层高度（米）
    layerHeight = 0.08;     % 每层高度（米）
    
    % 每个点的原始层号（可能有跳号）
    origLayerIdx = floor((zCoords - zMin) / layerHeight) + 1;
    
    % 过滤（可选，通常不必要，但稳妥）
    validMask = ~isnan(origLayerIdx) & origLayerIdx>=1;
    pts = points(validMask,:);
    origLayerIdx = origLayerIdx(validMask);
    
    % 把原始层号压缩为 1..K（没有点的原始层号会被移除）
    [uniqueLayers, ~, ic] = unique(origLayerIdx);   % uniqueLayers = 原始层号数组；ic 映射为 1..K
    
    % 按压缩后的组聚合点（不会产生中间空 cell）
    layerPoints = accumarray(ic, (1:numel(ic))', [], @(rows){ pts(rows,:) });
    
    % 对应密度
    layerDensities = cellfun(@(c) size(c,1)/layerHeight, layerPoints);
    
    % 保留映射： layerPoints{i} 对应原始层号 uniqueLayers(i)
    origLayerNums = uniqueLayers;
    
    numLayers = size(layerPoints, 1);

    %% 2.提取每层聚类
    minpts  = 6;         % DBSCAN参数:聚类中心点的最小邻居数量，DBSCAN 的发明者建议 minPts ≥ 2 * 数据维度
    
    allClusters = cell(numLayers,1);
    allCenters  = cell(numLayers,1);
    allHulls    = cell(numLayers,1);
    allIdx      = cell(numLayers,1);

    for i = 1:numLayers
        if any(isnan(layerPoints{i}(:)) | isinf(layerPoints{i}(:)))
            fprintf('Layer %d contains NaN or Inf in layerPoints\n', i);
        end
    
        if size(layerPoints{i}, 1) < minpts
            fprintf('Layer %d has %d points, less than minpts = %d\n', ...
                i, size(layerPoints{i}, 1), minpts);
        end
    
        eps = kdistplot3d(layerPoints{i}, minpts);
        if isnan(eps) || isinf(eps) || eps <= 0
            fprintf('Layer %2d: Invalid eps = %f \n', i, eps);
        end
    
        [clusters, centers, hulls, idxList] = extract_layer_clusters(double(layerPoints{i}), 0.02, minpts);
        allClusters{i} = clusters;
        allCenters{i}  = centers;
        allHulls{i}    = hulls;
        allIdx{i}      = idxList;
    end

    %% 3.建立层间连接相关参数
    % w1 ~ w4 为相应权重，其和等于 1
    params.w1 = 0.20;           % 位置d，设置高则强调两个聚类的距离是否较为接近
    params.w2 = 0.10;           % 方向theta，两个聚类中心连线与Z轴夹角
    params.w3 = 0.30;           % 形状ov，设置高则强调两个聚类的面积是否较多重叠
    params.w4 = 0.10;           % 规模，点数/密度，评价两个聚类的一致性
    params.w5 = 0.30;           % 方向dir_theta向量夹角权重
    
    params.d_max   = 0.26;      % 距离阈值，单位：米, 较大 → 允许更大横向偏移;较小 → 只保留距离很近的连接
    params.dz_max  = layerHeight * 1.0;     % Δz 允许偏差
    params.theta_max = 80;      % 角度阈值（度）
    params.dir_theta_max = 30;
    
    params.d_inlayer_max = 0.5;    % 同层补偿最大距离,同层聚类大于该值进行补偿连接，小于不连接
    
    connections = connect_clusters(allClusters, allCenters, allHulls, params);
    
    bottomLayers = 1:min(2,length(allClusters));
%     bottomLayers = roots;
    [~, ~, ~, pLabels, numPlants, ~] = build_plant_graph_v4(connections, allClusters, allCenters, bottomLayers, 0.7, 2, 'linear');
    fprintf('[build_plant_graph] 最终植株数: %d\n', numPlants);

    %% plantLabels 反向映射到点云，得到每株的点云集合
    plant_sets = assign_points_to_plants(points, pLabels, allClusters, allIdx, numPlants);

    %% 计算凸包信息
    plant_hulls = cell(length(plant_sets), 1);
    for i = 1:length(plant_sets)
        pts = plant_sets{i};

        try
            K = convhull(pts(:,1), pts(:,2));
            plant_hulls{i} = pts(K, 1:2);
        catch
            plant_hulls{i} = [];
        end
    end
end

function eps = kdistplot3d(pts, minPts)
    % 在 k-距离图中寻找"拐点"（elbow point）作为 eps 的值
    % pts    : N×3 double
    % minPts : DBSCAN 的 minPts

    Ksafe = min(minPts, size(pts, 1)-1);   % 安全的K值，防止数据集点数比近邻数小
    if Ksafe < 3    % 只有两个点无法求二阶导数
        eps = 0;
        return
    end
    [~, dist] = knnsearch(pts, pts, 'K', Ksafe+1); % 注意 K=k+1（因为包含自身）
    kDist = dist(:, Ksafe+1); % 第k近邻的距离（跳过自身）
%     idx   = knnsearch(pts, pts, 'K', Ksafe);
%     kDist = sqrt(sum((pts - pts(idx(:,Ksafe),:)).^2, 2));
    kDistSorted = sort(kDist, 'descend');
        
    % 简易肘点:“最大二阶导数”近似
    d1 = diff(kDistSorted);
    [~, elbowIdx] = max(abs(diff(d1)));
    
    eps = kDistSorted(elbowIdx + 1);            % +1 是为了把二阶差分索引还原到原数组坐标
end

function [clusters, centers, hulls, idxList] = extract_layer_clusters(pts, eps, minpts)
% 提取层聚类中心、凸包等信息
% 输入:
%   pts    : 当前层点云 Nx3
%   eps    : DBSCAN 半径
%   minpts : DBSCAN 最小点数
% 输出:
%   clusters : cell，每个聚类的点集
%   centers  : Kx3，每个聚类的几何中心
%   hulls    : cell，每个聚类在 XY 投影的凸包
%   idxList  : cell，每个聚类对应的点索引

labels = dbscan(pts, eps, minpts);
if all(labels == -1)
    fprintf('<== All points are noise\n');
% else
%     fprintf('\n');
end

clusters = {};
centers = [];
hulls = {};
idxList = {};

for c = unique(labels(labels~=-1))'     % 所有聚类的类别编号向量
    idx = find(labels == c);
    cluster = pts(idx,:);
    if isempty(cluster)
        fprintf('Empty cluster detected for label %d\n', c);
    end

    clusters{end+1} = cluster;
    centers(end+1,:) = mean(cluster,1);
    idxList{end+1} = idx;

    % XY 投影凸包
    projXY = cluster(:,1:2);
    if size(projXY,1) > 2
        k = convhull(projXY(:,1), projXY(:,2));     % 凸包索引
        hulls{end+1} = projXY(k,:);                 % 构成凸包的点集（XY 二维）
    else
        hulls{end+1} = projXY; % 点或线情况
    end
end
end


function connections = connect_clusters(allClusters, allCenters, allHulls, params)
% 建立层间连接
% 输入:
%   allClusters : 每层聚类点集
%   allCenters  : 每层聚类中心
%   allHulls    : 每层聚类凸包
%   params      : 超参数 (结构体)
% 输出:
%   connections : cell，存储相邻层连接

numLayers = numel(allClusters);
connections = cell(numLayers-1,1);

for i = 1:numLayers-1
    clusters1 = allClusters{i};
    centers1  = allCenters{i};
    hulls1    = allHulls{i};

    clusters2 = allClusters{i+1};
    centers2  = allCenters{i+1};
    hulls2    = allHulls{i+1};

    % 空层处理
    if isempty(clusters1) || isempty(clusters2) || isempty(centers1) || isempty(centers2)
        connections{i} = struct('src', {}, 'dst', {}, 'score', {}, 'type', {});
        continue;
    end

    % 计算层间距离矩阵
    D = pdist2(centers1(:,1:2), centers2(:,1:2), 'euclidean');               % 使用 pdist2 计算所有聚类中心的欧式距离

    % 预过滤符合距离阈值聚类
    valid = D < params.d_max;

    % 计算方向向量与置信度
    [dir1, conf1] = cellfun(@compute_cluster_direction, clusters1, 'UniformOutput', false);
    [dir2, conf2] = cellfun(@compute_cluster_direction, clusters2, 'UniformOutput', false);

    connList = struct('src', {}, 'dst', {}, 'score', {}, 'type', {});

    % ============ 层间正常连接 ============
    for j = 1:numel(clusters1)
        bestIdx = -1; bestScore = -inf; bestOV = 0;
        for k = 1:numel(clusters2)
%             fprintf("%2d-%2d-%2d\n", i, j, k);

            if ~valid(j,k)
                continue;
            end

            % ---- 调用独立评分函数 ----
            [score, ov] = compute_cluster_connection_score(...
                        clusters1{j}, clusters2{k}, ...
                        centers1(j,:), centers2(k,:), ...
                        hulls1{j}, hulls2{k}, ...
                        dir1{j}, conf1{j}, dir2{k}, conf2{k}, params);

            if score > bestScore
                if ov > 0                           % 情况1：存在重叠，直接更新最佳连接
                    bestScore = score;
                    bestIdx = k;
                    bestOV = ov;
                elseif bestOV <= 0                  % 情况2：无重叠，但当前未找到过重叠连接（bestOV <= 0）
                    bestScore = score;
                    bestIdx = k;
                    bestOV = 0;
                                                    % 情况3：无重叠，且已有更好重叠连接 → 不更新
                end
            end
        end
        
        connList(end+1).src   = j;              % i 层第 j 个聚类索引

        if bestIdx > 0
            connList(end).dst     = bestIdx;        % 连接目标 聚类索引
            connList(end).score   = bestScore;
            connList(end).type    = 'normal';
        else
            connList(end).dst     = 0;        % 下一层聚类索引
            connList(end).score   = 0;
            connList(end).type    = 'invalid';      % 无有效连接
        end
    end

    % ======== 无有效连接聚类同层补偿 ========
    used_src = [connList(strcmp({connList.type}, 'normal')).src];   % 仅统计已经正常连接的源聚类
    invalid_src = setdiff(1:numel(clusters1), used_src);
    
    if ~isempty(invalid_src)
        D_inlayer = pdist2(centers1(:,1:2), centers1(:,1:2));
        D_inlayer(eye(size(D_inlayer))==1) = inf; % 去掉自身
        for j = invalid_src
            [minDist, nearIdx] = min(D_inlayer(j,:));
            if minDist > params.d_inlayer_max && minDist < 1.0
                connList(end+1).src = j;
                connList(end).dst = nearIdx;
                connList(end).score = exp(-minDist/params.d_inlayer_max);
                connList(end).type = 'inlayer';
            end
        end
    end

    
    if all([connList.score])
        connections{i} = connList;
    else
        mask = [connList.score] ~= 0;

        connections{i} = connList(mask);
    end
end
end

function [dir, confidence] = compute_cluster_direction(cluster)
    % 计算更稳健的聚类方向向量及置信度
    % cluster: N×3 点云坐标
    % 输出:
    %   dir: 方向向量（1×3）
    %   confidence: 置信度 [0~1]
    
    n = size(cluster, 1);
    if n == 0
        dir = [];
        confidence = 0;
        return;
    end
    
    % ---------- 情况1：点数太少 (<3) ----------
    if n <= 3
        % 尝试用 PCA 主方向近似
        [coeff, ~, ~] = pca(cluster);
        if isempty(coeff)
            dir = [];
            confidence = 0;
            return;
        end
        dir = coeff(:,1)'; % 第一主成分方向
        dir = dir / norm(dir);
        if dir(3) < 0
            dir = -dir;
        end
        % 置信度与点数成比例
        confidence = min(0.3 + 0.1*n, 0.5); % 小聚类置信度低
        return;
    end
    
    % ---------- 情况2：点数足够 ----------
    % 法向量平均法
    pc = pointCloud(cluster);
    normals = pcnormals(pc, min(10, n-1));
    
    % 法向量分布稳定性
    dir = mean(normals, 1);
    dir = dir / (norm(dir) + 1e-10);
    if dir(3) < 0
        dir = -dir;
    end
    
    % 根据法向量一致性和点数定义置信度
    ang_dispersion = mean(vecnorm(normals - dir, 2, 2));
    confidence = max(0, 1 - ang_dispersion) * min(1, log10(n + 1) / 2);
    confidence = min(confidence, 1);
end

function [score, ov] = compute_cluster_connection_score(clusterA, clusterB, centerA, centerB, hullA, hullB, dirA, confA, dirB, confB, params)
% 计算两个聚类间的匹配评分
%
% 输入：
%   clusterA, clusterB : 聚类点集合
%   centerA, centerB   : 聚类中心 [x y z]
%   hullA, hullB       : 聚类凸包点（Nx3）
%   dirA, dirB         : 聚类方向向量
%   confA, confB       : 聚类方向置信度
%   params             : 权重参数结构体
%
% 输出：
%   score              : 两个聚类间的匹配评分

    % --- 1. 距离特征 ---
    d = norm(centerA(1:2) - centerB(1:2)); % XY距离

    % --- 2. 重叠特征 ---
    pgA = polyshape(hullA(:,1), hullA(:,2), 'Simplify', true);
    pgB = polyshape(hullB(:,1), hullB(:,2), 'Simplify', true);
    ov = polygonOverlap(hullA, hullB, pgA, pgB);

    % --- 3. 垂直夹角特征 ---
    vec = centerB - centerA;
    theta = acosd(abs(vec(3)) / (norm(vec) + 1e-10));

    % --- 4. 聚类方向一致性 ---
    if isempty(dirA) || isempty(dirB) || norm(dirA) == 0 || norm(dirB) == 0
        dir_theta = 90;  % 或 params.dir_theta_max，看你的设计意图
    else
        cosTheta = abs(dot(dirA(:), dirB(:))) / (norm(dirA) * norm(dirB) + 1e-10);  % cosTheta ∈ [0, 1]
        cosTheta = min(max(cosTheta, 0), 1);  % 确保在 [0,1]
        dir_theta = acosd(min(max(cosTheta, -1), 1));                               % dir_theta ∈ [0, 90]
    end

    % --- 5. 方向置信度 ---
    dir_theta_conf = sqrt((confA^2 + confB^2)/2);

    % --- 6. 点数相似度 ---
    sizeSim = min(size(clusterA,1), size(clusterB,1)) / ...
              max(size(clusterA,1), size(clusterB,1));

    % --- 7. 综合评分 ---
    if ov > 0
        score = params.w1 * exp(-d / params.d_max) + ...
                params.w2 * exp(-(theta / params.theta_max)^2) + ...
                params.w3 * ov + ...
                params.w4 * sizeSim + ...
                params.w5 * exp(-(dir_theta / params.dir_theta_max)) * dir_theta_conf;
    else
        score = (params.w1 + params.w3 / 3) * exp(-d / params.d_max) + ...
                (params.w2 + params.w3 / 3) * exp(-(theta / params.theta_max)^2) + ...
                (params.w4 + params.w3 / 3) * sizeSim + ...
                params.w5 * exp(-(dir_theta / params.dir_theta_max)) * dir_theta_conf;
    end
end

function ov = polygonOverlap(poly1, poly2, pg1, pg2)
% 计算两个多边形在 XY 平面上的交并比（IoU）
% 输入：
%   poly1, poly2 : N×2, M×2 多边形顶点坐标
%   pg1, pg2     : 预计算的 polyshape 对象（可选）
% 输出：
%   ov : 交并比

if isempty(poly1) || isempty(poly2) || size(poly1,1) < 3 || size(poly2,1) < 3
    ov = 0; return;
end

if any(isnan(poly1(:)) | isinf(poly1(:))) || any(isnan(poly2(:)) | isinf(poly2(:)))
    warning('polygonOverlap: Invalid input (NaN/Inf).');
    ov = 0; return;
end

% 边界框过滤
bb1 = [min(poly1); max(poly1)];
bb2 = [min(poly2); max(poly2)];
if bb1(2,1) < bb2(1,1) || bb2(2,1) < bb1(1,1) || ...
   bb1(2,2) < bb2(1,2) || bb2(2,2) < bb1(1,2)
    ov = 0; return;
end

try
    if nargin < 3
        pg1 = polyshape(poly1(:,1), poly1(:,2), 'Simplify', true);
        pg2 = polyshape(poly2(:,1), poly2(:,2), 'Simplify', true);
    end
    interPg = intersect(pg1, pg2);
    unionPg = union(pg1, pg2);
    interArea = area(interPg);
    unionArea = area(unionPg);
    minArea = 1e-6;

    if unionArea < minArea
        ov = 0;
    else
        ov = interArea / unionArea;
    end
catch E
    warning(E.identifier, 'polygonOverlap: Error (%s) for poly1 (%d points), poly2 (%d points).', ...
        E.message, size(poly1,1), size(poly2,1));
    ov = 0;
end
end


function [G, nodeList, nodeInfo, plantLabels, numPlants, stemRoots, stemMembers] = build_plant_graph_v4(connections, allClusters, allCenters, bottomLayers, xyTol, maxGapLayer, interpMethod)
% Update from BUILD_PLANT_GRAPH_V3
% 插值补全后基于 stem 重分配（仅建立子图与 stem 最近点对的连接）
%
% 输入：
%   connections    - 层间连接信息，cell 数组，每层包含 src 和 dst
%   allClusters    - 每层的聚类索引，cell 数组
%   allCenters     - 每层聚类的中心坐标，cell 数组
%   bottomLayers   - 底层索引
%   xyTol          - XY 平面容差
%   maxGapLayer    - 最大层间间隙
%   interpMethod   - 插值方法（仅支持 'linear'）
%
% 输出：
%   G             - 图对象
%   nodeList      - 节点列表，[layer, cluster]
%   nodeInfo      - 节点坐标
%   plantLabels   - 植物分组
%   numPlants     - 植物数量
%   stemRoots     - 茎秆根节点
%   stemMembers   - 茎秆成员节点

% 输入验证
assert(numel(allClusters) == numel(allCenters), 'allClusters and allCenters must have same length');
% % assert(all(cellfun(@(x,y) size(x,1)==size(y,1), allClusters, allCenters)), 'Cluster and center sizes must match');
% assert(all(cellfun(@(x, y) size(x, 1) == size(y, 1), allClusters, allCenters)), 'Cluster and center sizes must match');
% assert(strcmpi(interpMethod, 'linear'), 'Only linear interpolation is supported');

numLayers = numel(allClusters);

%% Step 1: 构建节点列表（预分配 + 高效索引）
totalNodes = sum(cellfun(@numel, allClusters));
nodeList = cell(totalNodes, 1);
nodeInfo = zeros(totalNodes, size(allCenters{1}, 2));
maxClusters = max(cellfun(@numel, allClusters));
nodeIdxMap = zeros(numLayers, maxClusters); % 替代 containers.Map
globalIdx = 1;
for L = 1:numLayers
    clusters = allClusters{L};
    for c = 1:numel(clusters)
        nodeList{globalIdx} = [L, c];
        nodeInfo(globalIdx, :) = allCenters{L}(c, :);
        nodeIdxMap(L, c) = globalIdx;
        globalIdx = globalIdx + 1;
    end
end
numNodes = totalNodes;

%% Step 2: 创建空图并添加节点
G = graph();
G = addnode(G, numNodes);

%% Step 3: 添加已有连接（批量）
if ~isempty(connections)
    edges = [];
    for L = 1:length(connections)
        conn = connections{L};
        for k = 1:length(conn)
            i1 = nodeIdxMap(L, conn(k).src);

            if strcmp(conn(k).type, 'inlayer')
                % 同层连接
                i2 = nodeIdxMap(L, conn(k).dst);
            else
                i2 = nodeIdxMap(L+1, conn(k).dst);
            end
            
            if i1 > 0 && i2 > 0
                edges = [edges; i1, i2, norm(nodeInfo(i1,:) - nodeInfo(i2,:))];
            end
        end
    end
    if ~isempty(edges)
        existing = findedge(G, edges(:,1), edges(:,2));
        G = addedge(G, edges(existing==0,1), edges(existing==0,2), edges(existing==0,3));
    end
end

%% Step 4: 识别底层茎秆
% [stemRoots, stemMembers, numStems] = identify_bottom_stems(allCenters, bottomLayers, xyTol);
[stemRoots, stemMembers, numStems] = identify_bottom_stems_old(allClusters, allCenters, bottomLayers, xyTol);

% 转换为节点索引
stemNodeIdxSets = cell(numStems, 1);
stemAnchors = zeros(numStems, 1);
for s = 1:numStems
    mem = stemMembers{s}; % K x 2 [layer, cluster]
    idxs = zeros(size(mem,1), 1);
    valid = mem(:,1) <= numLayers & mem(:,2) <= maxClusters;
    idxs(valid) = nodeIdxMap(sub2ind(size(nodeIdxMap), mem(valid,1), mem(valid,2)));
    stemNodeIdxSets{s} = idxs(idxs > 0);
    
    r = stemRoots(s, :); % [layer, cluster]
    if r(1) <= numLayers && r(2) <= maxClusters
        stemAnchors(s) = nodeIdxMap(r(1), r(2));
    end
end
stemAnchors = stemAnchors(stemAnchors > 0);

%% Step 5: 插值补全跨层连接
bins = conncomp(G); % 缓存连通分量
numSub = max(bins);
currentMaxNodeIdx = numNodes;
if maxGapLayer > 1
    for sg = 1:numSub
        subIdx = find(bins == sg);
        if numel(subIdx) < 2, continue; end
        
        % 按层排序
        layers = cell2mat(nodeList(subIdx));
        subLayers = layers(:,1);
        [sortedLayers, si] = sort(subLayers);
        sortedNodes = subIdx(si);
        
        newEdges = [];
        for p = 1:length(sortedNodes)-1
            layerA = sortedLayers(p);
            layerB = sortedLayers(p+1);
            gap = layerB - layerA;
            if gap > 1 && gap <= maxGapLayer
                idA = sortedNodes(p);
                idB = sortedNodes(p+1);
                ptA = nodeInfo(idA, :);
                ptB = nodeInfo(idB, :);
                
                for g = 1:(gap-1)
                    t = g / gap;
                    ptInterp = (1-t) * ptA + t * ptB; % 线性插值
                    currentMaxNodeIdx = currentMaxNodeIdx + 1;
                    
                    % 添加新节点
                    nodeInfo(currentMaxNodeIdx, :) = ptInterp;
                    nodeList{currentMaxNodeIdx} = [layerA + g, 0]; % 插值节点 cluster=0
                    G = addnode(G, 1);
                    
                    % 收集边
                    newEdges = [newEdges; idA, currentMaxNodeIdx, norm(ptA - ptInterp)];
                    idA = currentMaxNodeIdx;
                end
                newEdges = [newEdges; idA, idB, norm(nodeInfo(idA,:) - nodeInfo(idB,:))];
            end
        end
        % 批量添加边
        if ~isempty(newEdges)
            existing = findedge(G, newEdges(:,1), newEdges(:,2));
            G = addedge(G, newEdges(existing==0,1), newEdges(existing==0,2), newEdges(existing==0,3));
        end
    end
end

%% Step 6: 基于茎秆重新分配子图
plantLabelsByStem = cell(numStems, 1);
for s = 1:numStems
    plantLabelsByStem{s} = zeros(0, 2); % 预分配空矩阵
end

for sg = 1:numSub
    subIdx = find(bins == sg);
    if isempty(subIdx), continue; end
    
    % 检查是否包含茎秆根节点
    presentStemIdx = find(ismember(stemAnchors, subIdx));
    if ~isempty(presentStemIdx)
        % 选择离质心最近的茎秆
        subCenter = mean(nodeInfo(subIdx, :), 1);
        minD = Inf;
        pickStem = presentStemIdx(1);
        for k = 1:length(presentStemIdx)
            sIdx = presentStemIdx(k);
            stemNodes = stemNodeIdxSets{sIdx};
            if isempty(stemNodes), continue; end
            D = pdist2(nodeInfo(stemNodes, :), subCenter, 'euclidean');
            dmin = min(D(:));
            if dmin < minD
                minD = dmin;
                pickStem = sIdx;
            end
        end
        plantLabelsByStem{pickStem} = [plantLabelsByStem{pickStem}; cell2mat(nodeList(subIdx))];
        continue;
    end
    
    % 子图不含茎秆，找最近点对
    bestDist = Inf; bestU = -1; bestV = -1; bestStem = -1;
    for sIdx = 1:numStems
        stemNodes = stemNodeIdxSets{sIdx};
        if isempty(stemNodes), continue; end
        D = pdist2(nodeInfo(subIdx, :), nodeInfo(stemNodes, :), 'euclidean');
        [dmin, pos] = min(D(:));
        [ui, vi] = ind2sub(size(D), pos);
        if dmin < bestDist
            bestDist = dmin;
            bestU = subIdx(ui);
            bestV = stemNodes(vi);
            bestStem = sIdx;
        end
    end
    
    % 添加边并分配子图
    if bestU > 0 && bestV > 0 && bestStem > 0
        if findedge(G, bestU, bestV) == 0
            G = addedge(G, bestU, bestV, bestDist);
        end
        plantLabelsByStem{bestStem} = [plantLabelsByStem{bestStem}; cell2mat(nodeList(subIdx))];
    else
        % 挂接到最近茎秆根节点
        if ~isempty(stemAnchors)
            D = pdist2(nodeInfo(subIdx, :), nodeInfo(stemAnchors, :), 'euclidean');
            [dmin, pos] = min(D(:));
            [ui, vi] = ind2sub(size(D), pos);
            bestU = subIdx(ui);
            bestV = stemAnchors(vi);
            bestStem = vi;
            if bestU > 0 && bestV > 0
                if findedge(G, bestU, bestV) == 0
                    G = addedge(G, bestU, bestV, dmin);
                end
                plantLabelsByStem{bestStem} = [plantLabelsByStem{bestStem}; cell2mat(nodeList(subIdx))];
            end
        end
    end
end

%% Step 7: 整理输出
plantLabels = {};
for s = 1:numStems
    if isempty(plantLabelsByStem{s}), continue; end
    [~, ia] = unique(plantLabelsByStem{s}, 'rows', 'stable');
    plantLabels{end+1} = plantLabelsByStem{s}(ia, :);
end
numPlants = length(plantLabels);

end

function [stemRoots, stemMembers, numRoots] = identify_bottom_stems(allCenters, rootCenters, xyTol)
    % identify_bottom_stems - 基于 rootCenters 查找对应的分层聚类中心作为茎秆根节点
    %
    % 输入：
    %   allClusters - 各层聚类点云 (cell)
    %   allCenters  - 各层聚类中心 (cell)
    %   rootCenters - 全部植株的根节点中心 (N×3)
    %   xyTol       - 匹配容差 (单位与坐标一致)
    %
    % 输出：
    %   stemRoots   - 每个根对应的 [layer, cluster]
    %   stemMembers - 每株对应的成员（此处仅包含根节点）
    %   numStems    - 根节点数量（即植株数量）
    
    % === 收集所有聚类中心及其层索引 ===
    allCentersMat = [];
    layerIdx = [];
    clusterIdx = [];
    for i = 1:numel(allCenters)
        if isempty(allCenters{i})
            continue;
        end
        n = size(allCenters{i}, 1);
        allCentersMat = [allCentersMat; allCenters{i}(:, 1:3)];
        layerIdx = [layerIdx; repmat(i, n, 1)];
        clusterIdx = [clusterIdx; (1:n)'];
    end
    
    numRoots = size(rootCenters, 1);
    stemRoots = zeros(numRoots, 2);
    stemMembers = cell(numRoots, 1);
    
    for r = 1:numRoots
        root = rootCenters(r, 1:2);
        dists = vecnorm(allCentersMat(:, 1:2) - root, 2, 2);
        [minDist, idx] = min(dists);
    
        % 若距离在阈值内则认为匹配成功，否则仍选最近的一个
        if minDist <= xyTol
            stemRoots(r, :) = [layerIdx(idx), clusterIdx(idx)];
        else
            fprintf('⚠ Root #%d 未在 xyTol=%.3f 内匹配，使用最近聚类 (%.3f)\n', r, xyTol, minDist);
            stemRoots(r, :) = [layerIdx(idx), clusterIdx(idx)];
        end
    
        % 可根据需要扩展成员集，目前仅含根节点自身
        stemMembers{r} = [layerIdx(idx), clusterIdx(idx)];
    end

end

function [stemRoots, stemMembers, numStems, stemOrder] = identify_bottom_stems_old(allClusters, allCenters, bottomLayers, xyTol)
%IDENTIFY_BOTTOM_STEMS
% 输入：
%   allClusters : cell{Nlayers,1}
%   allCenters  : cell{Nlayers,1}
%   bottomLayers: 选择用于分析的层号（如 1:5）
%   xyTol       : XY 投影距离阈值
%
% 输出：
%   stemRoots   : numStems x 2 矩阵，每行为 [layerID, clusterID]（每株茎秆最底层聚类）
%   stemMembers : cell(numStems,1)，每个单元是 Kx2 的矩阵，列为 [layerID, clusterID] 表示该茎柱在底层若干层的成员
%   numStems    : 识别到的茎秆数
%
% 说明：
%   - 函数按 layer 从低到高遍历，某茎柱首次出现的聚类即记录为该茎柱的 "根"。
%   - 匹配采用 XY 投影距离阈值 xyTol，可根据数据调整。
%   - 内部通过 PCA + 形状/密度阈值筛选候选聚类 (可按需要调整阈值)。

if nargin < 4, xyTol = 0.03; end

% 判定阈值（可按数据调）
minPoints   = 3;    % 聚类最少点数
pcaCosTol   = 0.01; % 第一主成分与Z轴余弦阈值（越靠近1越竖直）
ratioTol    = 0.7;  % latent(2)/latent(1) 的阈值（越小越柱状）
densityTol  = 2;   % 点密度阈值 (points / xy_area)，根据点云密度调整

clusterFeatures = cell(length(bottomLayers),1);

% Step 1.计算候选特征（使用 allCenters 的中心作为 xyCenter）
for li = 1:length(bottomLayers)
    i = bottomLayers(li);
    clusters = allClusters{i};
    if isempty(clusters)
        clusterFeatures{li} = [];
        continue;
    end
    nC = length(clusters);
    feats = struct('clusterID', cell(nC,1), 'xyCenter', cell(nC,1), 'isCandidate', cell(nC,1));
    for c = 1:nC
        pts = clusters{c};
        feats(c).clusterID = c;
        feats(c).xyCenter = allCenters{i}(c,1:2);
        if size(pts,1) < minPoints
            feats(c).isCandidate = false;
            continue;
        end
        % PCA + density
        try
            [coeff,~,latent] = pca(pts);
        catch
            feats(c).isCandidate = false;
            continue;
        end
        cosTheta = abs(dot(coeff(:,1), [0;0;1]));
        ratio = latent(2)/latent(1);
        xyRange = max(pts(:,1:2)) - min(pts(:,1:2));
        area = max(xyRange(1)*xyRange(2), eps);
        density = size(pts,1)/area;
        feats(c).isCandidate = (cosTheta >= pcaCosTol) && (ratio <= ratioTol) && (density >= densityTol);
    end
    clusterFeatures{li} = feats;
end

% Step 2.层间匹配形成茎柱（按 bottomLayers 的顺序遍历 => 首次出现即为最底层）
stemMembers = {};
stemLastXY = [];
for li = 1:length(bottomLayers)
    feats = clusterFeatures{li};
    if isempty(feats), continue; end
    i = bottomLayers(li); % 全局层号
    for c = 1:length(feats)
        if ~feats(c).isCandidate, continue; end
        xy = feats(c).xyCenter;
        if isempty(stemMembers)
            stemMembers{1} = [i, feats(c).clusterID];
            stemLastXY(1,:) = xy;
            continue;
        end
        dists = sqrt(sum((stemLastXY - xy).^2, 2));
        [mind, idx] = min(dists);
        if mind < xyTol
            stemMembers{idx} = [stemMembers{idx}; i, feats(c).clusterID];
            stemLastXY(idx,:) = xy;
        else
            newIdx = length(stemMembers) + 1;
            stemMembers{newIdx} = [i, feats(c).clusterID];
            stemLastXY(newIdx,:) = xy;
        end
    end
end

% Step 3: 提取最底层聚类
numStems = length(stemMembers);
stemRoots = zeros(numStems,2);
stemCenters = zeros(numStems,2);
for s = 1:numStems
    stemRoots(s,:) = stemMembers{s}(1,:); % 第一行就是最底层聚类
    i = stemRoots(s,1);
    c = stemRoots(s,2);
    stemCenters(s,:) = allCenters{i}(c,1:2); % 从 allCenters 取对应 XY
end

% Step 4: 按 Y->X 排序给茎编号
% ====== 1. 按X方向分组 ======
xTol = 0.5;  % X方向分组阈值，可根据实际点间距调整
[~, sortX] = sort(stemCenters(:,1));     % 按X初步排序
sortedXY = stemCenters(sortX,:);

groups = zeros(numStems,1);
groupID = 1;
groups(sortX(1)) = groupID;
for i = 2:numStems
    if abs(sortedXY(i,1) - sortedXY(i-1,1)) > xTol
        groupID = groupID + 1;  % 新组
    end
    groups(sortX(i)) = groupID;
end

% ====== 2. 组内再按Y排序 ======
stemOrder = zeros(numStems,1);
counter = 1;
for g = 1:groupID
    idx = find(groups == g);
    [~, subSort] = sort(stemCenters(idx,2), 'descend');  % 按Y升序
    for j = 1:length(idx)
        stemOrder(idx(subSort(j))) = counter;
        counter = counter + 1;
    end
end


end

function [plantPointClouds] = assign_points_to_plants(originalPointCloud, plantLabels, allClusters, allIdx, numPlants)
% 把连通子图结果映射回原始点云
%
% 输入:
%   plantLabels : cell(numPlants,1)，每株包含的 [layerID, clusterID]
%   allClusters : cell(numLayers,1)，每层聚类 {clusters}
%   allIdx      : cell(numLayers,1)，每层聚类 {indices}（保存原始点云索引，和 clusters 对应）
%   numPlants   : 植株数
% 输出:
%   plantPointClouds : cell(numPlants,1)，每株点云 (Npx3)
%   pointLabels      : 向量，和原始点云大小一致，每个点的植株 ID

% 获取原始点总数
% N = sum(cellfun(@(x) sum(cellfun(@(c) size(c,1), x)), allClusters));
N = size(originalPointCloud,1);   % ✅ 使用原始点云的总点数

plantPointClouds = cell(numPlants,1);

for p = 1:numPlants
    if isempty(plantLabels{p})
        fprintf('Plant %d: plantLabels empty\n', p);
        continue;
    end

    pts = [];
    idxs = [];
    for k = 1:size(plantLabels{p},1)
        layerID   = plantLabels{p}(k,1);
        clusterID = plantLabels{p}(k,2);
        if clusterID > 0 && clusterID <= length(allClusters{layerID})
            pts  = [pts;  allClusters{layerID}{clusterID}]; %#ok<AGROW>
            idxs = [idxs; allIdx{layerID}{clusterID}];      %#ok<AGROW>
        end
    end
    plantPointClouds{p} = pts;
end

end

