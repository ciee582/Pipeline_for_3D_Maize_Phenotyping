function P = rosa_lineextract(P, RADIUS, bUpdateConnectivity)
%
% we RECOVER CONNECTIVITY after each edge collapse, which is slow, but the
% result sampling is more uniform.
% bUpdateConnectivity: 0 for RECOVER CONNECTIVITY after each edge collapse, which is slow, but the
% result sampling is more uniform.
% bUpdateConnectivity: 1 for update CONNECTIVITY after each edge collapse, which is faster, and the
% result sampling is almost as same as that of the original method (bUpdateConnectivity=0).
%
% @author: Andrea Tagliasacchi
% reformed by jjcao
% @reform-data:     2009-9-3   using 1-ring to construct the connectivity
% graph
% @reform-data:     2010-8-20  decompose the original function into several
% small functions.
SHOW_RESULTS = false;

options.collapse_order = 1;
if (bUpdateConnectivity)
    % 建立 原始点 → 最近骨架点 的映射,每个 P.cpts 点只能属于一个骨架点,离谁最近，就归谁
    [P.spls,P.corresp] = farthest_sampling_by_sphere(P.cpts, P.sample_radius);

    % 如果两个原始点是邻居，那它们对应的骨架点也应该连边，关键代码 A(pc, nc) = 1;
    P.spls_adj = connect_by_inherit_neigh(P.cpts, P.spls, P.corresp, P.rings);

    % Optimization: Post-process corresp to refine assignments in potential branch regions (Direction 3)
%     [P.spls, P.corresp, P.spls_adj] = refine_corresp_in_branches(P.cpts, P.spls, P.corresp, P.rings, P.spls_adj, RADIUS);

    % 按代价（度数 + 距离）不断"去三角 → 得 1D 图”过程
    [P.spls, P.spls_adj,P.corresp] = edge_collapse_update(P.cpts, P.spls, P.corresp, P.rings, P.spls_adj, options);
   
else
    [P.spls,P.corresp] = farthest_sampling_by_sphere(P.cpts, P.sample_radius);
    P.spls_adj = connect_by_inherit_neigh(P.cpts, P.spls, P.corresp, P.rings);
    [P.spls, P.spls_adj,P.corresp] = edge_collapse(P.cpts, P.spls, P.corresp, P.rings, P.spls_adj, options);
end

%%
if (SHOW_RESULTS)
    figure; movegui('northeast');set(gcf,'color','white');hold on;
    plot3( P.spls(:,1), P.spls(:,2), P.spls(:,3), '.r', 'markersize', 5);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
    plot_connectivity(P.spls, P.spls_adj);
    view3d zoom;
end
end


function [spls, corresp, A] = refine_corresp_in_branches(pts, spls, corresp, neigh, A, RADIUS)
% Post-process corresp in potential branch regions (Direction 3: minimal change)
% Identify high-degree spls as branches and reassign corresp using local clustering
% Adjusted: Incorporate topology via distance field (inspired by DFSP) and PCA-based geometric features for stem-leaf distinction, handling wrapping like leaf collars

% Compute degrees (excluding self-loop)
degrees = sum(A, 2) - diag(A);  % 更安全的写法
degrees = degrees(:);           % 确保是列向量

% Find potential branch points: spls with degree > 2
branch_indices = find(degrees > 2);

if isempty(branch_indices)
    return;
end

kdtree_pts = kdtree_build(pts);
kdtree_spls = kdtree_build(spls);

% Assume z is height axis; estimate growth direction (stem axis) as global z for simplicity
% Find approximate stem base: lowest z spl
[~, stem_base_idx] = min(spls(:,3));
stem_base = spls(stem_base_idx,:);

% Precompute distance field: Euclidean distance from stem base (mimicking DFSP base point)
dist_field = zeros(size(pts,1),1);
for i = 1:size(pts,1)
    dist_field(i) = euclidean_distance(pts(i,:), stem_base);
end

% Precompute local geometric features (PCA eigenvalues for linearity/planarity)
geom_features = zeros(size(pts,1),3);  % e1, e2, e3 (sorted descending)
for i = 1:size(pts,1)
    if ~isempty(neigh{i})
        neigh_pts = pts(neigh{i}(1:min(20, length(neigh{i}))), :);
        if size(neigh_pts,1) >= 3
            cov_mat = cov(neigh_pts - pts(i,:));
            evals = sort(eig(cov_mat), 'descend');
            geom_features(i,:) = evals / sum(evals);  % Normalized eigenvalues
        end
    end
end
linearity = geom_features(:,1) - geom_features(:,2);  % High for stems
planarity = geom_features(:,2) - geom_features(:,3);  % High for leaves

% Precompute for spls (average from corresp pts)
spls_linearity = zeros(size(spls,1),1);
spls_planarity = zeros(size(spls,1),1);
for i = 1:size(spls,1)
    corresp_idx = find(corresp == i);
    if ~isempty(corresp_idx)
        spls_linearity(i) = mean(linearity(corresp_idx));
        spls_planarity(i) = mean(planarity(corresp_idx));
    end
end

% Precompute local densities for all spls (points in RADIUS ball / volume approx)
densities = zeros(size(spls,1), 1);
vol_approx = (4/3) * pi * RADIUS^3;  % Sphere volume approximation
for i = 1:size(spls,1)
    [nIdxs, ~] = kdtree_ball_query(kdtree_pts, spls(i,:), RADIUS);
    densities(i) = length(nIdxs) / vol_approx;
end
median_density = median(densities(densities > 0));  % Global median for normalization

for b_idx = branch_indices'
    % 1. 先找当前 spl 对应的原始点
    branch_pts_idx = find(corresp == b_idx);
    if length(branch_pts_idx) < 3
        continue;  % 点太少，不值得重分配
    end

    % 2. Adaptive multiplier: high density (stem-like) -> larger radius to smooth; low density (leaf-like) -> smaller to precise
    local_density = densities(b_idx);
    if local_density > median_density
        multiplier = 1.5;  % Larger for dense stems to reduce invalid branches
    else
        multiplier = 0.5;  % Smaller for sparse leaves, per user test
    end
    local_radius = multiplier * RADIUS;

    % 取以该 spl 为中心、adaptive local_radius 的局部邻域
    [local_idx, ~] = kdtree_ball_query(kdtree_pts, spls(b_idx,:), local_radius);
    if length(local_idx) < 3
        continue;
    end

    local_pts = pts(local_idx, :);
    local_dist_field = dist_field(local_idx);
    local_linearity = linearity(local_idx);
    local_planarity = planarity(local_idx);

    % 3. 确定聚类数目：不超过当前 spl 的度数，也不能超过局部点数
    num_clusters = min(max(2, degrees(b_idx)), size(local_pts,1));

    if num_clusters <= 1
        continue;
    end

    % 4. 进行 k-means 聚类，整合距离场和几何特征以考虑拓扑和包裹结构
    feature_data = [local_pts, local_dist_field * 0.2, ...  % Scale distance field
                    local_linearity * local_radius, local_planarity * local_radius];  % Scale geom features
    try
        [cluster_labels, cluster_centers] = kmeans(feature_data, num_clusters, ...
            'MaxIter', 200, 'Replicates', 3, 'Display', 'off');
        cluster_centers = cluster_centers(:,1:3);  % Extract position centers
    catch ME
        warning('kmeans failed for branch %d: %s', b_idx, ME.message);
        continue;
    end

    % 5. 为每个簇重新分配到最合适的相邻 spl（优先选择与当前 spl 相连的邻居）
    neigh_spls = find(A(b_idx, :) > 0);
    neigh_spls = neigh_spls(neigh_spls ~= b_idx);  % 去掉自己

    for c = 1:num_clusters
        cluster_idx = find(cluster_labels == c);
        if isempty(cluster_idx)
            continue;
        end
        cluster_pt_idx = local_idx(cluster_idx);

        % 计算簇中心到所有邻居 spl 的距离，添加几何相似性以处理包裹（茎线性，叶平面）
        if ~isempty(neigh_spls)
            dist_to_neigh = zeros(length(neigh_spls), 1);
            cluster_lin = mean(local_linearity(cluster_idx));
            cluster_plan = mean(local_planarity(cluster_idx));
            for n = 1:length(neigh_spls)
                xy_dist = euclidean_distance(cluster_centers(c,1:2), spls(neigh_spls(n),1:2));
                lin_diff = abs(cluster_lin - spls_linearity(neigh_spls(n))) * local_radius;
                plan_diff = abs(cluster_plan - spls_planarity(neigh_spls(n))) * local_radius;
                dist_to_neigh(n) = xy_dist + lin_diff + plan_diff;
            end
            [~, min_n] = min(dist_to_neigh);
            new_spl_idx = neigh_spls(min_n);
        else
            % 极端情况：没有邻居，就保持原对应（不应该发生）
            new_spl_idx = b_idx;
        end

        % 重新分配
        corresp(cluster_pt_idx) = new_spl_idx;
    end
end

kdtree_delete(kdtree_pts);
kdtree_delete(kdtree_spls);

% 最后重建邻接矩阵 A
A = connect_by_inherit_neigh(pts, spls, corresp, neigh);
end

function dist = euclidean_distance(p1, p2)
    v= p1 - p2;
    dist = sqrt(dot(v, v));
end
