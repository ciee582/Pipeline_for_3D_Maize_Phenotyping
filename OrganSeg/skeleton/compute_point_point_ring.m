function rings = compute_point_point_ring(pts, k, index)
% 返回的ring是所有点的一环面的各个顶点的索引
% 该索引是按顺序存储的
% pts: n*3 matrix for coordinates where we want compute 1-ring.
% k: k of kNN
% index: index of kNN
%
% example:
% M.rings = compute_point_point_ring(M.verts, 1, M.k_knn);
% M.rings = compute_point_point_ring(M.verts, [1:size(M.verts,1)], M.k_knn);

% 输入验证
if isempty(pts) || size(pts, 2) ~= 3
    error('无效输入：pts 必须是 n*3 矩阵');
end
if k <= 1 || ~isnumeric(k)
    error('k 必须是大于 1 的正整数');
end

% 如果未提供 index，计算 k-NN
npts = size(pts,1);
if nargin < 3 || isempty(index)
%     kdtree = kdtree_build(pts);% kdtree,用来找k近邻
%     index = zeros(npts, k);
%     for i = 1:npts
%         index(i,:)  = kdtree_k_nearest_neighbors(kdtree,pts(i,:),k)';
%          %index(i,:)  = flipud( kdtree_k_nearest_neighbors(kdtree,pts(i,:),k))';
%          % flipud 对矩阵或数组进行上下翻转
%     end

    [index, ~] = knnsearch(pts, pts, 'K', k);
end

rings = cell(npts, 1);

% 处理每个点
for i = 1:npts
    rings{i} = index(i, 2:end); % 默认回退，确保即使后续计算失败，ring{i} 仍然有一个默认值。
    neighbors = pts(index(i,:), :);

    % 将 3D 坐标投影到 2D 平面（由前两个主成分定义）
    [coeff, ~, ~] = pca(neighbors);
    % 投影到局部主平面
    proj2D = neighbors * coeff(:,1:2);
    
    % Delaunay剖分 + robust
    TRI = delaunayn(proj2D);
    
    % 一环邻域提取
    % 找与中心点（索引为 1）关联的所有三角形
    [row, ~] = find(TRI == 1);
    faces = sort(TRI(row,:), 2);
    faces = faces(:, 2:end);
    proj2D = sort(faces(:));

    d = diff([proj2D; max(proj2D)+1]);   % 计算新向量中相邻元素之间的差值。非零元素表示 x 中数值发生变化的位置
    count = diff(find([1;d])); % x 中每个连续出现的唯一数值的频数
    y = [proj2D(find(d)) count];    % x(find(d)) 表示 x 中每个唯一值第一次出现的位置
    n_sorted_index = size(y, 1);
    start = find(count == 1);
    if ~isempty(start)         % 如果有只出现一次的顶点，视为起点，否则任取一个
        want_to_find = y(start(1),1);
    else
        want_to_find = faces(1,1);
        n_sorted_index = n_sorted_index + 1;   % 首尾是封闭的环
    end
    
    j = 0;
    sorted_index = zeros(1, n_sorted_index);
    while j < n_sorted_index
        j = j + 1;
        sorted_index(j) = want_to_find;
        [row, col] = find(faces == want_to_find);
        if ~isempty(col)
            if col(1) == 1
                want_to_find = faces(row(1), 2);
                faces(row(1), 2) = -1;
            else
                want_to_find = faces(row(1),1);
                faces(row(1), 1) = -1;
            end
        end
    end

    % 首位为端点，如果邻域面是封闭的，则首位数字相同，否则不同    
    rings{i} = index(i, sorted_index); 
end
end
