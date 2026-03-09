function [spls, corresp, spls_adj, joints, segments] = remove_small_cycles_mt(spls, corresp, spls_adj, joints, roots, segments, threshold, show_cycles)
% removing small cycles measured by topological length <= threshold
% Joints will not be deleted or moved; only non-joint points may be merged

G = spls_adj;
cycles = findcycles_undirected(sparse(G));% can be speed up by just using joints, not every nodes.

joint_ids = joints(:,1);

for i = 1:length(cycles)
    cycle = cycles{i};
    if length(cycle) < 3
        continue;
    end

    if show_cycles
        tmp = [cycle, cycle(1)];
        figure('Name', 'Show cycles', 'NumberTitle', 'off', 'color', 'white');
        for j = 1:(length(tmp)-1)
            p1 = spls(tmp(j), :);
            p2 = spls(tmp(j+1), :);
            plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], ...
                  'Color', [1 0 0], 'LineWidth', 1.5, 'LineStyle', '-');
            hold on;
        end
        axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-180,0);view3d rot;
    end

    % Large cycle: break one joint-joint edge
    if length(cycle) > threshold
        % 优先尝试断开几何上距离最远且分布相对均衡的 joint-joint 边
        joint_joint_edges  = [];
        for j = 1:length(cycle)
            u = cycle(j);
            v = cycle(mod(j,length(cycle))+1);
            if (ismember(u, joint_ids) && ismember(v, joint_ids))
                joint_joint_edges(end+1,:) = [u,v];
            end
        end

        % 获取各 joint 外延分支深度
        depth_map = compute_joint_branch_depth(spls_adj, joint_ids, cycle);
    
        scores = inf(size(joint_joint_edges,1),1);
        for k = 1:size(joint_joint_edges,1)
            u = joint_joint_edges(k,1);
            v = joint_joint_edges(k,2);
            d1 = 0; d2 = 0;
            if isKey(depth_map, u), d1 = depth_map(u); end
            if isKey(depth_map, v), d2 = depth_map(v); end
            scores(k) = d1 + d2; % 越小越易断
        end
    
        [~, best_idx] = min(scores);
        if isfinite(scores(best_idx))
            u = joint_joint_edges(best_idx,1);
            v = joint_joint_edges(best_idx,2);
            spls_adj(u,v) = 0;
            spls_adj(v,u) = 0;
            fprintf('选择断开外延最弱的 joint-joint 边: %d-%d\n', u, v);
        else
            warning('所有 joint-joint 边断开后仍未断开环');
        end
        continue;
    end

    % Small cycle
    joint_in_cycle = intersect(cycle, joint_ids);
    non_joint = setdiff(cycle, joint_ids);
    if isempty(non_joint)
        continue;
    end

    % 1. Remove joint–nonjoint edges inside the cycle
    for j = 1:length(joint_in_cycle)
        jid = joint_in_cycle(j);
        neighbors = find(spls_adj(jid,:));
        for n = neighbors
            if ismember(n, non_joint)
                spls_adj(jid,n) = 0;
                spls_adj(n,jid) = 0;
            end
        end
    end

    % 2. Remove non-joint to external edges
    for j = 1:length(non_joint)
        nid = non_joint(j);
        neighbors = find(spls_adj(nid,:));
        for n = neighbors
            if ~ismember(n, cycle)
                spls_adj(nid,n) = 0;
                spls_adj(n,nid) = 0;
            end
        end
    end

    % 3. Reconnect external nodes that lost all edges
    deleted_ids = non_joint;
    rep_candidates = joint_in_cycle;
    if isempty(rep_candidates)
        rep_candidates = non_joint(1);
    end

    for j = 1:length(non_joint)
        nid = non_joint(j);
        neighbors = find(G(nid,:));
        externals = setdiff(neighbors, cycle);
        for e = externals(:)'
            if all(spls_adj(e,:) == 0) % lost all connections
                preserved_coords = spls(rep_candidates,:);
                e_coord = spls(e,:);
                dists = vecnorm(preserved_coords - e_coord, 2, 2);
                [~, idx] = min(dists);
                nearest = rep_candidates(idx);
                spls_adj(e,nearest) = 1;
                spls_adj(nearest,e) = 1;
            end
        end
    end

    % 4. Remove non-joint points from graph, do NOT touch joint nodes
    for j = 1:length(non_joint)
        nid = non_joint(j);
        spls(nid,:) = NaN;
        spls_adj(nid,:) = 0;
        spls_adj(:,nid) = 0;
        segments(nid) = 0;
        corresp(corresp == nid) = joint_in_cycle(1); % map to one joint
    end
end

end


function depth_map = compute_joint_branch_depth(spls_adj, joints, cycle_nodes)
% 返回每个 joint 节点向环外延伸的最大深度
% joints: Nx1 向量或 Nx2，第一列是 joint ID
% cycle_nodes: 构成当前环的所有点

    if size(joints,2) == 2
        joint_ids = joints(:,1);
    else
        joint_ids = joints;
    end

    n = size(spls_adj,1);
    is_in_cycle = false(n,1);
    is_in_cycle(cycle_nodes) = true;

    depth_map = containers.Map('KeyType','double','ValueType','double');

    for i = 1:length(joint_ids)
        joint = joint_ids(i);

        % 向环外搜索
        visited = false(n,1);
        visited(joint) = true;
        max_depth = 0;

        stack = [joint, 0]; % [node, depth]

        while ~isempty(stack)
            node = stack(end,1);
            d    = stack(end,2);
            stack(end,:) = [];

            neighbors = find(spls_adj(node,:) == 1);
            for nb = neighbors
                if visited(nb) || is_in_cycle(nb)
                    continue
                end
                visited(nb) = true;
                max_depth = max(max_depth, d+1);
                stack = [stack; nb, d+1];
            end
        end

        depth_map(joint) = max_depth;
    end
end


function cycles = findcycles_undirected(G)
% 无向图中的环检测，返回所有简单环

    numNodes = size(G, 1);
    visited = false(numNodes, 1);
    cycles = {};

    for startNode = 1:numNodes
        if ~visited(startNode)
            dfs(startNode, -1, []);
        end
    end
    
    function dfs(u, parent, path)
        visited(u) = true;
        path = [path, u];
        neighbors = find(G(u,:));
        for v = neighbors
            if v == parent
                continue;  % 不回头
            end
            if ~ismember(v, path)
                dfs(v, u, path);
            else
                % 找到闭环
                idx = find(path == v, 1);
                cycle = path(idx:end);
                cycle = unique(cycle, 'stable');
                % 防止重复环
                is_duplicate = any(cellfun(@(c) isequal(sort(c), sort(cycle)), cycles));
                if ~is_duplicate
                    cycles{end+1} = cycle;
                end
            end
        end
    end

end

