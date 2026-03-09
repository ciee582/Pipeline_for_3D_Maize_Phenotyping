function [shortest_dist, shortest_path] = mydijkstra(adj_matrix, start_point, end_point)
% MYDIJKSTRA_OPTIMIZED 寻找 i,j 两点最短路径的优化版Dijkstra算法
% 输入：
%   adj_matrix  — 邻接矩阵，adj_matrix(i,j) 指 i 到 j 之间的距离
%                 如果 i 到 j 没有边，则 adj_matrix(i,j) 设为 inf
%                 可以是有向图。
%   start_point — 起始顶点编号
%   end_point   — 目标顶点编号
%
% 输出：
%   shortest_dist  — 最短路的距离。如果不可达，则为 inf。
%   shortest_path  — 最短路的路径。如果不可达，则为空数组 []。
%
% 优化说明：
%   1. 使用优先队列（通过MATLAB的排序功能模拟）来高效查找下一个访问的顶点。
%      这显著提高了算法在稀疏图上的性能。
%   2. 增加了对输入参数的初步校验，提高健壮性。
%   3. 明确了对负权重边的处理警告。

% --- 输入参数校验 ---
ptsNum = size(adj_matrix, 1);

if start_point < 1 || start_point > ptsNum || end_point < 1 || end_point > ptsNum
    error('起始点或结束点超出邻接矩阵的有效范围。');
end

if any(adj_matrix < 0, 'all')
    warning('邻接矩阵包含负权重边。Dijkstra算法不适用于负权重图，结果可能不正确。');
end

% --- 初始化 ---
distance = inf(1, ptsNum);      % 起点到各顶点的当前最短距离
distance(start_point) = 0;      % 起点到自身距离为0

parent = zeros(1, ptsNum);      % 记录最短路径中的前驱顶点

% 优先队列：存储 [距离, 顶点] 对，按照距离从小到大排列。
% 初始时只包含起点。
% 这里我们使用一个N*2的矩阵来模拟优先队列，第一列是距离，第二列是顶点编号。
% 每次取出最小值后，将其标记为已处理。
priority_queue = [0, start_point]; % [distance, vertex_id]

% --- Dijkstra 主循环 ---
while ~isempty(priority_queue)
    % 从优先队列中取出距离最小的顶点
    % sortrows 将按第一列（距离）升序排列，然后取第一个元素
    % 或者更直接地，找到距离最小的行
    [~, min_idx] = min(priority_queue(:, 1));
    
    current_dist = priority_queue(min_idx, 1);
    u = priority_queue(min_idx, 2);
    
    % 将已取出的顶点从队列中移除
    priority_queue(min_idx, :) = []; 

    % 如果当前取出的距离比已知的距离要大，说明这个节点已经被处理过或者有更短的路径，跳过
    % 这是为了处理优先队列中可能存在的重复项（当一个节点的距离被多次更新并加入队列时）
    if current_dist > distance(u)
        continue; 
    end

    % 遍历所有与 u 相邻的顶点 v
    for v = 1:ptsNum
        % 确保 adj_matrix(u,v) 不是 inf (表示有边) 且 u 不是 v 自己 (避免自环)
        if adj_matrix(u, v) ~= inf && adj_matrix(u,v) ~= 0 && u ~= v % 假设0也表示无边
             % 松弛操作：如果通过 u 到 v 的路径更短
            if distance(u) + adj_matrix(u, v) < distance(v)
                distance(v) = distance(u) + adj_matrix(u, v);
                parent(v) = u;
                % 将更新后的顶点及距离加入优先队列
                priority_queue = [priority_queue; distance(v), v];
            end
        end
    end
end

% --- 构建最短路径 ---
shortest_path = [];
shortest_dist = distance(end_point);

if shortest_dist ~= inf % 如果目标点可达
    t = end_point;
    shortest_path = end_point;
    while t ~= start_point
        % 如果在某个环节发现parent(t)为0，说明路径断裂，通常是由于算法逻辑问题或负权重边导致
        if parent(t) == 0
            shortest_path = []; % 清空路径，表示路径无法构建
            break; 
        end
        p = parent(t);
        shortest_path = [p shortest_path];
        t = p;
    end
end

end