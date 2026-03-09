function [spls, corresp, spls_adj, joints, segments] = remove_irrelevant_extrama(spls, corresp, spls_adj, joints, root_id, segments, global_dist, threshold, SHOW_IRRELEVANT_EXTRAMA)
    tmp = global_dist(global_dist ~= Inf);
    skel_size = max(tmp) - min(tmp);

    while true
        edges = []; % 存储候选的 [joint, leaf, normalized_distance]
        
        for i = 1:size(joints,1)
            joint_id = joints(i,1);
            neighbors = find(spls_adj(joint_id,:) == 1);
            neighbors(neighbors == joint_id) = [];  % 去除自身
            
            for j = neighbors
                if sum(spls_adj(j,:) == 1) == 1  % 是叶子节点
                    if joint_id == root_id || j == root_id
                        continue  % 跳过 root 节点本身或其连接点，防止误删根部结构
                    end
                    normalized_dist = abs(global_dist(joint_id) - global_dist(j)) / skel_size;
                    edges(end+1,:) = [joint_id, j, normalized_dist];
                end
            end
        end

        if isempty(edges), break; end

        edges = sortrows(edges, 3);  % 按 normalized_dist 升序排序
        if edges(1,3) > threshold, break; end

        % 删除最短的叶子分支
        joint_id = edges(1,1);
        leaf_id = edges(1,2);

        segments(leaf_id) = 0;
        corresp(corresp == leaf_id) = joint_id;
        spls(leaf_id, :) = NaN;

        % 更新邻接矩阵（对称）
        spls_adj(leaf_id, :) = 0;
        spls_adj(:, leaf_id) = 0;
        spls_adj(joint_id, leaf_id) = 0;
        spls_adj(leaf_id, joint_id) = 0;

        % 如果 joint 不再是 joint（连接数过少），移除它
        links = find(spls_adj(joint_id,:) == 1);
        if length(links) < 4
            segments(joint_id) = NaN;
            joints(joints(:,1) == joint_id, :) = [];
        end
    end
    
    % draw joints, draw roots; 
    if SHOW_IRRELEVANT_EXTRAMA
        figure('Name','Remove irrelevant extrama','NumberTitle','off');set(gcf,'color','white');
        movegui('south');hold on;
        plot_skeleton(spls, spls_adj);
        scatter3(spls(joints(:,1),1),spls(joints(:,1),2),spls(joints(:,1),3),200,'b','filled');
        scatter3(spls(root_id,1),spls(root_id,2),spls(root_id,3),30,'r','filled');
        axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
    end

end