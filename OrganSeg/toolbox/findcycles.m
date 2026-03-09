function [cycles]=findcycles(G)
% for directed graph
% 
% example:
% G = sparse([1 2 3 4 5 5 5 6 7 8 8 9],[2 4 1 5 3 6 7 9 8 1 10 2],true,10,10);
% view(biograph(G));
% findcycles(G);

numNodes = size(G,1); 
cycles = cell(0);

for n = 1:numNodes
    % 将图构造成 digraph（有向图），如果是无向图可改为 graph
    DG = digraph(G);

    % 深度优先搜索，从节点 n 开始
    T = dfsearch(DG, n, 'allevents');

    % 建立一个 map 来保存每个节点的前驱节点（相当于 P 向量）
    pred = containers.Map('KeyType', 'double', 'ValueType', 'double');

    for i = 1:size(T,1)
        if strcmp(T{i,1}, 'discover')
            u = T{i,2};
            v = T{i,3};
            if ~isnan(v)
                pred(v) = u;
            end
        end
    end

    % 检查是否存在从某个节点 d 返回 n 的边（即存在环）
    for d = 1:numNodes
        if G(d,n)
            % 存在从 d 到 n 的边
            if isKey(pred, d)
                % 构造从 n 到 d 的路径
                path = d;
                while path(1) ~= n && isKey(pred, path(1))
                    path = [pred(path(1)), path];
                end
                if path(1) == n
                    % 构造闭环 path + n
                    cycles{end+1} = [path, n];
                end
            end
        end
    end

    % 移除 n 的出边（仿照原代码的 G(n,:) = 0）
    G(n,:) = 0;
end
