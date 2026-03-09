function [organ_subskeleton, skeletonType] = SkeletonDecomposition2(mainskeleton, spls, joints, roots, spls_adj, corresp, features, show_results)
% 输出 organ_subskeleton，其中最后一个是茎秆，其余是叶子
% skeletonType（暂未使用）保留接口

%%%%%%%%%%%%% 初始化 %%%%%%%%%%%%%
rootNum = size(roots,1);
jointNum = size(joints,1);
s_adj_Num = size(spls,1);
adjmatrix = zeros(s_adj_Num, s_adj_Num);

% 构造邻接矩阵
for i = 1:s_adj_Num
    for j = 1:s_adj_Num
        temp = spls_adj(i,j);
        if i == j
            adjmatrix(i,j) = 0;
        elseif temp == 0
            adjmatrix(i,j) = inf;
        else
            adjmatrix(i,j) = 1;
        end
    end
end

%%%%%%%%%%%%% 构建叶片子骨架 %%%%%%%%%%%%%
organ_subskeleton = {};
sindex = 1;

for i = 1:rootNum
    r_index = roots(i);
    if ~ismember(mainskeleton, r_index)
        continue;
    end

    Dist = inf(jointNum, 1);
    Path = cell(jointNum, 1);

    for j = 1:jointNum
        c_index = joints(j);
        if ~ismember(mainskeleton, c_index)
            continue;
        end
        [t1, t2] = mydijkstra(adjmatrix, r_index, c_index);
        if ~isempty(t2)
            Dist(j) = t1;
            Path{j} = t2;
        end
    end

    [~, a] = min(Dist);
    b = Path{a};
    if isempty(b)
        continue; % 跳过无法到达的路径
    end

    organ_subskeleton{sindex} = b(:);
    sindex = sindex + 1;
end

%%%%%%%%%%%%% 找出唯一非叶片（茎秆）子骨架 %%%%%%%%%%%%%
% ==== 计算每个子骨架的方向向量、平均 feature 和长度特征 ====
d_v = zeros(length(organ_subskeleton), 3);
skefeature = [];
skelen = [];

for i = 1:length(organ_subskeleton)
    indices = organ_subskeleton{i};
    
    if length(indices) < 3
        skefeature = [skefeature; -Inf]; 
        d_v(i,:) = [0, 0, 0];
        skelen = [skelen; 0];
        continue;
    end

    skelen = [skelen; length(indices) / length(spls)];
    
    % 平均 feature（例如：点云密度、粗细等）
    pids = ismember(corresp, indices);
    fsum = sum(features(pids));
    psum = sum(pids);
    mf = fsum / max(psum, 1);  % 避免除0
    skefeature = [skefeature; mf];

    vi_end = indices(end);
    vi_mid = indices(floor(end / 2));
    v_end = spls(vi_end, :);
    v_mid = spls(vi_mid, :);
    
    dir_vec = v_mid - v_end;
    if norm(dir_vec) == 0
        d_v(i,:) = [0, 0, 0];
    else
        d_v(i,:) = dir_vec / norm(dir_vec);
    end
end

% ==== 使用综合得分进行 stem 判定 ====
% 组合多个特征（你也可以调整权重）
% ==== 主茎应与 Z 轴近似平行 ====
z_axis_candidates = [0 0 1; 0 0 -1];
z_scores = abs(d_v * z_axis_candidates');  % Nx2
best_z_score = max(z_scores, [], 2);  % 每行取 max

dir_similarity = best_z_score;  % 越接近 1 表示越垂直

% ==== 综合得分（加权可调）====
score = ...
    0.2 * normalize(skefeature) + ...
    0.2 * normalize(skelen) + ...
    0.6 * dir_similarity;  % 越接近垂直，得分越高

[~, index] = max(score);

% ==== 提取主茎根节点 ====
if isempty(organ_subskeleton{index})
    error('Stem sub-skeleton is empty. Unable to determine stem root.');
end

stem_root = organ_subskeleton{index}(1);  % 主茎起点
organ_subskeleton(index) = [];            % 从叶片中移除主茎骨架

%%%%%%%%%%%%% 重新计算主茎骨架并放在最后 %%%%%%%%%%%%%
Dist = inf(jointNum, 1);
Path = cell(jointNum, 1);

for j = 1:jointNum
    c_index = joints(j);
    if ~ismember(mainskeleton, c_index)
        continue;
    end
    [t1, t2] = mydijkstra(adjmatrix, stem_root, c_index);
    if ~isempty(t2)
        Dist(j) = t1;
        Path{j} = t2;
    end
end

[~, a] = max(Dist);
b = Path{a};
if isempty(b)
    error('无法从 stem_root 到任意 joint 构建主茎骨架。');
end

organ_subskeleton{end+1} = b(:);
skeletonType = []; % 可扩展为每个子骨架的类型标识

if show_results
    disp(['成功识别到 ', num2str(length(organ_subskeleton) - 1), ' 片叶子和 1 个主茎']);
end
  
  
%%%%%%%%Optimize the 3D position of stem vertices  
%   for i=3:length(b)
%      index1=b(i-2); 
%      index2=b(i-1);
%      index3=b(i);
%      spl_1=spls(index1,:);
%      spl_2=spls(index2,:);
%      spl_3=spls(index3,:);
%      r=norm(spl_3-spl_2);
%      dir=(spl_2-spl_1)./norm(spl_2-spl_1);
%      spl_4=spl_2+dir*r;
%      spls(index3,:)=0.5*(spl_4+spl_3);
%   end
 
 %%%%%%%%%%%sort leaf %%%%%%%%%%%%%%%%%%%%
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_results
figure('Name','SkeletonDecomposetion','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
   for i=1:length(organ_subskeleton)
       indices=organ_subskeleton{i};
       cr=rand(1,1);
       cg=rand(1,1);
       cb=rand(1,1);
       color=[cr cg cb];
         scatter3(spls(indices,1),spls(indices,2),spls(indices,3),20,color, 'filled');
       hold on;
       if(i==length(organ_subskeleton))
           color=[0 0 0];
              scatter3(spls(indices,1),spls(indices,2),spls(indices,3),20,color, 'filled');
              hold on;
              scatter3(spls(indices(1),1),spls(indices(1),2),spls(indices(1),3),100,color, 'filled');  
              hold on;
     end
    
       
       hold on;
   end
   axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
end

end

