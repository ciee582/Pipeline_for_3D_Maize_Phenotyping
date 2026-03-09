function [ P_S ] = skeleton_laplacian( P, Parameters )
try
    [~, name, ~] = fileparts(P.fileName);

    outputPath = fullfile(P.pathName, [name, '_skeleton.mat']);
    
    % 若相关数据文件已存在则跳过后续执行过程，直接加载文件中的数据
    if isfile(outputPath)
        load(outputPath, 'spls', 'corresp', 'spls_adj', 'joints', 'roots', 'branches');
        P.spls = spls;            % 骨架点云，代表骨架结构的关键点（n * 3）
        P.corresp = corresp;      % 原始点云（pts）中的每个点与采样后的骨架点（spls）之间的对应关系,一个 n x 1 的数组
        P.spls_adj = spls_adj;    % 骨架点的邻接矩阵，表示骨架点之间的连接关系
        P.joints = joints;        % 骨架拓扑结构中的关键节点，如分支点或连接点
        P.roots = roots;          % 骨架的端点，骨架的末段点或拓扑结构的参考点      
        P.branches = branches;    % 骨架的分支结构，即骨架中从一个关节点到另一个关节点（或端点）的路径段

        P_S = P;

        fprintf('存在文件 %s 对应的骨架缓存数据: %s，已成功加载！\n', P.fileName, outputPath);

        return;
    else
        fprintf('对应的骨架缓存数据不存在，将开始计算。。。。。。');
        
        % 启动并行池
%         pool = gcp('nocreate');
%         if isempty(pool)
%             parpool;  % 自动按CPU核心数开足线程
%         end
end
catch ME
    warning('处理文件 %s 出错: %s', P.fileName, ME.message);
end


%% 预处理与参数初始化
P.npts = size(P.pts, 1);

[P.bbox, P.diameter] = GS.compute_bbox(P.pts);
% 运行时间: 0.001225 秒

P.k_knn = Parameters.KnnNum;

options.USING_POINT_RING = GS.USING_POINT_RING;

% 计算点邻域
if options.USING_POINT_RING
    P.rings = compute_point_point_ring(P.pts, P.k_knn, []);
    % 运行时间: 9.910287 秒
else
    P.frings = compute_vertex_face_ring(P.faces);
    P.rings = compute_vertex_ring(P.faces, P.frings);
end

% 基于拉普拉斯算子的点云收缩算法
[P.cpts, ~, ~, ~, ~] = contraction_by_mesh_laplacian(P, options);
% 把原始 3D 结构向中轴收缩,让枝干、叶片更接近"中心线”,为 骨架 抽取服务
% 运行时间: 89.179449 秒

% 下采样 (控制骨架点数量)
P.sample_radius = P.diameter * Parameters.sampleScale;
P = rosa_lineextract(P, P.sample_radius, 1);
% 运行时间: 6.103151 秒

% 参数准备
t1 = Parameters.t1;     % 内部分支点距离阈值
a1 = Parameters.a1;     % 内部分支角度阈值
t2 = Parameters.t2;     % 去除无关极值阈值
t3 = Parameters.t3;     % 小环路移除阈值

SHOW_JOINTS = false;
SHOW_CYCLES = false;
SHOW_ROOT_JOINT = false;
SHOW_IRRELEVANT_EXTRAMA = false;

%% ====================== 骨架优化流程 ===================================

% 初始骨架提取后的关节点、端点、分支
[joints, roots, segments] = find_joints(P.pts, P.spls, P.corresp, P.spls_adj, SHOW_JOINTS);
% 运行时间: 0.005097 秒

%% 1. 移除拓扑长度过小的环路
[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_small_cycles_mt(P.spls, P.corresp, P.spls_adj, joints, roots, segments, t3, SHOW_CYCLES);
% 运行时间: 0.303719 秒

%% 2. 计算根节点 global distance relative to root_id, "size of skeleton"
[root_id, global_dist] = find_root_node(P.spls, P.spls_adj, joints, SHOW_ROOT_JOINT);
% 运行时间: 0.004178 秒

%% 3. 合并过近关节点 (可选优化) The effect is bad for some cases.
if ~isempty(joints)
    [P.spls, P.corresp, P.spls_adj, joints, root_id] = merge_nearby_joints(P.spls, P.corresp, P.spls_adj, joints, root_id);
end
% 运行时间: 0.001919 秒

%% 4. 去除无效极值
[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_irrelevant_extrama(P.spls, P.corresp, P.spls_adj, joints, root_id, segments, global_dist, t2, SHOW_IRRELEVANT_EXTRAMA);
% P.root_id = root_id;  % The first node is root node if there is no joints.
% 运行时间: 0.002478 秒

%% 5. 重新构图
[P.spls,P.corresp,P.spls_adj, graph] = build_graph(P.spls, P.corresp, P.spls_adj);% may contain cycles
% 运行时间: 0.038704 秒

%% 6. 进一步平滑骨架位置（对孤立点或中间点求均值修正）
if Parameters.smoothSkeleton
    for i = 1:size(P.spls,1)
        verts = P.pts(P.corresp==i, :);
        if size(verts,1) == 1
            P.spls(i,:) = verts;
        else
            P.spls(i,:) = mean(verts);
        end
    end
    for i = 1:size(P.spls_adj,1)
        links = find( P.spls_adj(i,:) == 1 );
        if length(links) == 2
            verts = P.spls(links,:);
            P.spls(i,:) = mean(verts);
        end
    end
end

%% 7. 重新检测关节点、端点、分支点
[joints ,roots, branches] = find_Joints_mt(P.spls, P.spls_adj, P.pts, SHOW_JOINTS);
% 运行时间: 0.245268 秒
P.joints = joints;
P.roots = roots;
P.branches = branches;

P_S = P;

% 将重要结果写入缓存，避免重复计算
if (Parameters.saveSkeleton)
    spls = P_S.spls;
    corresp = P_S.corresp;
    spls_adj = P_S.spls_adj;
    joints = P_S.joints;
    roots = P_S.roots;
    branches = P_S.branches;
    try
        save(outputPath, 'spls', 'corresp', 'spls_adj', 'joints', 'roots', 'branches');
    catch e
        error('无法保存文件 %s：%s', outputPath, e.message);
    end
    fprintf('骨架结果已保存至缓存: %s\n', outputPath);
end

end

