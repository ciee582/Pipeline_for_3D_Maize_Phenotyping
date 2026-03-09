function [cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options)
% point cloud contraction_by_mesh_laplacian
% refer to Skeleton Extraction by Mesh Contraction 08
% mesh could be simplified by MeshLab/Filter/Clustering Decimation, if it
% is too larged.
% 该函数基于拉普拉斯收缩，通过迭代求解线性系统 A \ b，将点云向其骨架方向移动。
% 只调整点的位置，而不直接删除或合并点
% inputs:
%   P.pts
%   P.faces
%   P.npts
%   P.k_knn: k of knn
%   P.rings:
%   options.WL = 3;     初始拉普拉斯约束权重
%   options.WC = 1;     初始吸引权重
%   options.sl:         每次迭代中 WL 的缩放因子
%   options.tc:         终止条件（面积或体积比率）。
%   options.iterate_time = 10;  最大迭代次数
%
% outputs:
%   cpts:               收缩后的点云坐标
%   t                   迭代次数
%
% 输入验证
if isfield(P, 'fileName') && isfield(P, 'pathName')
    [~, name, ~] = fileparts(P.fileName);
    mat_file = fullfile(P.pathName, [name '-laplacianContractionData.mat']);
    
    if exist(mat_file, 'file') == 2
        data = load(mat_file);
        required_vars = {'cpts', 't', 'initWL', 'WC', 'sl'};
        if all(isfield(data, required_vars))
            load(mat_file, required_vars{:});
            
            fprintf('已加载本地拉普拉斯收缩数据：%s\n', mat_file);
            return
        else
            fprintf('警告：%s 缺少必要变量，需重新计算拉普拉斯收缩过程！\n', mat_file);
        end
    else
        fprintf('提示：数据文件不存在，需重新计算拉普拉斯收缩过程！\n');
    end
else
    error('文件路径传输错误！');
end

if ~isfield(P, 'pts') || size(P.pts, 2) ~= 3
    error('P.pts 必须是 n*3 矩阵');
end
if options.USING_POINT_RING && ~isfield(P, 'rings')
    error('P.rings 必须提供当 USING_POINT_RING 为 true');
end

%##########################################################################
%% 参数设置
%##########################################################################
RING_SIZE_TYPE = 1;             % 决定邻域大小的计算方式 1:min, 2:mean, 3:max
Laplace_type = 'conformal';     % 拉普拉斯算子类型（conformal, combinatorial, spring, mvc）

tc = getoptions(options, 'tc', GS.CONTRACT_TERMINATION_CONDITION);
iterate_time = getoptions(options, 'iterate_time', GS.MAX_CONTRACT_NUM); 
initWL = getoptions(options, 'WL', GS.compute_init_laplacian_constraint_weight(P, Laplace_type)); 

% 根据不同类型的离散值，设置初始吸引力权重
if strcmp(Laplace_type, 'mvc')
    WC = getoptions(options, 'WC', 1)*10;
elseif strcmp(Laplace_type,'conformal')
    WC = getoptions(options, 'WC', 1);
else
    WC = getoptions(options, 'WC', 1);
end

% 点的位置约束权重，初始化为 WC 的常数向量
WH = ones(P.npts, 1) * WC; 
% 每次迭代中，拉普拉斯权重（WL）的缩放因子
sl = getoptions(options, 'sl', GS.LAPLACIAN_CONSTRAINT_SCALE);
WL = initWL;    %*sl;


% 初次计算拉普拉斯矩阵
L = compute_laplacian(P, P.pts, Laplace_type, options);

%% 初次收缩
A = [L * WL; sparse(1:P.npts, 1:P.npts, WH)];       % sparse(1:P.npts, 1:P.npts, WH) 对角稀疏矩阵，其中对角线上的元素都是 WC
b = [zeros(P.npts, 3); sparse(1:P.npts, 1:P.npts, WH) * P.pts];
cpts = A \ b;

% 计算初始邻域大小
if options.USING_POINT_RING
    sizes = GS.one_ring_size(P.pts, P.rings, RING_SIZE_TYPE);   % min radius of 1-ring
    size_new = GS.one_ring_size(cpts, P.rings, RING_SIZE_TYPE);
    a(1) = sum(size_new)/sum(sizes);
else 
    ratio_new = area_ratio_1_face_ring(P.pts, cpts, P.faces, P.frings);
%     ratio = ones(size(ratio_new));
    a(1) = mean(ratio_new);%sum(ratio_new)/sum(ratio);
end

%% 迭代收缩
t = 1;
while t < iterate_time
    % 更新拉普拉斯矩阵
    L = compute_laplacian(P, cpts, Laplace_type, options);
    
    % 更新拉普拉斯权重
    WL = min(sl * WL, GS.MAX_LAPLACIAN_CONSTRAINT_WEIGHT);

    % 更新权重
    if options.USING_POINT_RING
        size_new = GS.one_ring_size(cpts, P.rings, RING_SIZE_TYPE);
        WH = WC * (sizes ./ size_new);
        if strcmp(Laplace_type, 'mvc')
            WH = WH * 10;
        end
    else
        ratio_new = area_ratio_1_face_ring(P.pts, cpts, P.faces, P.frings);
        WH = WC * (ratio_new .^ -0.5);
    end
    WH(WH > GS.MAX_POSITION_CONSTRAINT_WEIGHT) = GS.MAX_POSITION_CONSTRAINT_WEIGHT;
    
    % 求解
    A = real([L * WL; sparse(1:P.npts, 1:P.npts, WH)]);     % real 提取复数的实部
    b(P.npts+1:end, :) = sparse(1:P.npts, 1:P.npts, WH) * cpts;
    tmp = A \ b;
    
    % 更新比率
    if options.USING_POINT_RING
        size_new = GS.one_ring_size(tmp, P.rings, RING_SIZE_TYPE);
        a(end+1) = sum(size_new) / sum(sizes);
    else
        ratio_new = area_ratio_1_face_ring(P.pts, tmp, P.faces, P.frings);
        a(end+1) = mean(ratio_new);
    end

    % 检查终止条件
    % 若边界框尺寸超过原始 1.2 倍，退出
    tmpbox = GS.compute_bbox(tmp);
    if sum( (tmpbox(4:6)-tmpbox(1:3)) > ((P.bbox(4:6)-P.bbox(1:3))*1.2) ) > 0 || a(t) - a(end) < tc || isnan(a(end))
        break;
    end

    cpts = tmp;    
    t = t + 1;  
end

% % 保存相关参数至本地
% try
%     save(mat_file, 'cpts', 't', 'initWL', 'WC', 'sl');
% catch e
%     error('无法保存文件 %s：%s', mat_file, e.message);
% end
% fprintf('拉普拉斯收缩计算结果保存至：%s\n', mat_file);

% 可视化（可选）
if getoptions(options, 'SHOW_CONTRACTION_PROGRESS', false)
    visualize_contraction(P, cpts, WH, WL, t);
    figure;
    plot(a); xlabel('迭代次数'); ylabel('体积比率');
    title('收缩过程体积变化');
end

end

function L = compute_laplacian(P, pts, Laplace_type, options)
    if options.USING_POINT_RING
        L = -compute_point_laplacian(pts, Laplace_type, P.rings, options);
    else
        L = -compute_mesh_laplacian(pts, P.faces, Laplace_type, options);
    end
end

function visualize_contraction(P, cpts, WH, WL, t)
    figure('Name', sprintf('迭代 %d 次', t), 'Color', 'white');
    scatter3(P.pts(:,1), P.pts(:,2), P.pts(:,3), 10, WH, 'filled');
    colormap(parula); colorbar;
    hold on;
    scatter3(cpts(:,1), cpts(:,2), cpts(:,3), 10, 'r', 'filled');
    axis off; axis equal; view3d rot;
    title(sprintf('迭代 %d 次，WL = %.2f', t, WL));
    legend('原始点云', '收缩点云');
end
