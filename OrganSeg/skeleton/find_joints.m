function [joints, roots, segments] = find_joints(pts, spls, corresp, A, show_joints)
% 识别骨架图中关键节点类型
% joints: 多分支连接点（连接数 > 3）
% roots: 端点或叶子节点（连接数 < 3）
% segments: 每个骨架点所属的“分支段”编号（用负数索引），方便后续分段或分析。

joints = zeros(0,2);                % [index, diameter] [骨架索引, 局部直径]
roots = [];                         % 端点索引
segments = NaN(size(spls, 1), 1);   % 每个点对应的分支段编号（负数为joint）
jidx = 0;                           % 分支编号（负号编码）

% ----- 1. 根节点识别 -----
% 假设根节点为 Z 最低点
[~, root_id] = min(spls(:,3));
roots = root_id;

% ----- 2. 遍历骨架点 -----
for i = 1:length(segments)
    if ~isnan(segments(i)), continue, end 

    if isnan(spls(i, 1))      % 骨架点为空，失效点
        segments(i) = 0;
        continue;
    end

    links = find( A(i,:)==1 );

    % ----- 2a. 复杂关节点 -----
    if length(links) > 3        
        jidx = jidx - 1;
        segments(i) = jidx;

        % 骨架点覆盖的原始点云尺度
        [~, diam] = GS.compute_bbox(pts(corresp==i, :));
        joints(end+1, :) = [i, diam];
    end

    % ----- 2b. 端点识别 -----
    % 原逻辑：length(links)<3
    % 增强：包含根节点 + 根部附近点
    isRootCandidate = (length(links)<3) || (i==root_id);
    
    % 可选：根部附近阈值，例如 Z < rootZ + h_thresh
    h_thresh = 0.1 * (max(spls(:,3)) - min(spls(:,3))); % 10% 高度范围
    if spls(i,3) <= spls(root_id,3) + h_thresh
        isRootCandidate = true;
    end

    if isRootCandidate
        jidx = jidx - 1;
        if ~ismember(i, roots)
            roots = [roots; i];
        end
    end


if show_joints
    figure('Name','Find joints Result','NumberTitle','off');set(gcf,'color','white');movegui('northwest');
    plot_skeleton(spls, A);hold on;
    scatter3(spls(joints(:,1),1),spls(joints(:,1),2),spls(joints(:,1),3),200,'b', 'filled');
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
end

end