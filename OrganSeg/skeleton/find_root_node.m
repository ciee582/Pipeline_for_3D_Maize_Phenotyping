function [root_id, global_dist] = find_root_node(spls, spls_adj, joints, show_root)
    if isempty(joints)
        root_id = 1;
    else
        joint_z = spls(joints(:,1), 3);
        [~, idx] = min(joint_z);
        root_id = joints(idx,1);
    end

    % 构造局部距离矩阵
    local_dist = zeros(size(spls_adj)); 
    for i = 1:size(spls_adj,1)
        for j = 1:size(spls_adj,2)
            if spls_adj(i,j) == 1
                rs = spls(i,:) - spls(j,:);
                local_dist(i,j) = norm(rs);
            end
        end
    end

    % 构建图并计算全局距离
    G = graph(local_dist .* spls_adj);  % 保留结构和权重的图
    global_dist = distances(G, root_id);  % 返回从 root_id 到所有点的最短路径距离

    % 可视化根节点
    if show_root
        figure('Name','Find joints and ROOTS','NumberTitle','off');set(gcf,'color','white');
        plot_skeleton(spls, spls_adj);hold on;
        scatter3(spls(joints(:,1),1),spls(joints(:,1),2),spls(joints(:,1),3),200,'b', 'filled');
        scatter3(spls(root_id,1), spls(root_id,2), spls(root_id,3), 400, 'r', 'filled');
        axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
    end
end
