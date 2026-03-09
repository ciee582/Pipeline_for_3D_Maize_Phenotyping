function [ Regions ] = FineSegmentation( points, bottomStem, leafSeeds, unSegment, StopX, sample1, sample2, coeff_ratio, coeff_plane, DebugShow )
%   Inputs:
%             points：   整个植株的点云数据，N×3。
%         bottomStem：   植株茎秆底部的参考点，用于计算与未分割点的距离。
%          leafSeeds：   粗分割阶段的初始叶片区域种子，每个叶片对应一个 cell，包含该叶片的点索引。
%          unSegment：   未被分配到任何叶片的点索引。
%              StopX：   每个叶片的 Z 方向截断高度，未达到该高度的点不被考虑为叶片的候选。
%            sample1：   距离计算中使用的样本点数，供 caldis1 使用。
%            sample2：   PCA法投影距离计算中使用的样本点数，供 caldis2 使用。
%        coeff_ratio：   几何近邻距离阈值
%        coeff_plane：   局部领域平面投影距离贡献比例

    leafNum = length(leafSeeds);
    UnSegPts = points(unSegment, :);
    junctionIdx = [];     % 保存叶-茎交界点索引
    
    % ========= 1. 基于基部距离排序（由远及近，优先处理容易划分的点） =========
    EFF = sum((UnSegPts - bottomStem).^2, 2);
    [~, sortIdx] = sort(EFF, 'descend');
    
    for ii = 1:length(sortIdx)
    
        idx = sortIdx(ii);
        pt  = UnSegPts(idx, :);
    
        %% ========= 2. caldis1：局部空间邻近度 =========
        Dis1 = inf(leafNum, 1);
    
        for j = 1:leafNum
            ids = leafSeeds{j};
    
            % StopX 生长顺序约束
            if pt(3) < StopX(j)
                continue;
            end
    
            leafPts = points(ids, :);
            Dis1(j) = caldis1(leafPts, pt, sample1);
        end
    
        [sortDis1, ord] = sort(Dis1, 'ascend');
        s1 = ord(1); 
        s2 = ord(2);
        d1 = sortDis1(1); 
        d2 = sortDis1(2);
    
        %% ========= 3. 若差距明显，直接归属 =========
        if abs(d2 - d1) / max(d1, eps) > coeff_ratio
            leafSeeds{s1} = [leafSeeds{s1}; unSegment(idx)];
            continue;
        end
    
        %% ========= 4. caldis2（局部平面一致性 + 可信度） =========
        [d1_plane, valid1] = caldis2(points(leafSeeds{s1}, :), pt, sample2);
        [d2_plane, valid2] = caldis2(points(leafSeeds{s2}, :), pt, sample2);
    
        %% ========= 5. 决策逻辑 =========
        % 情况 A：至少一个叶片局部是“可信平面”
        if valid1 || valid2
    
            if ~valid1, d1_plane = inf; end
            if ~valid2, d2_plane = inf; end
    
            score1 = d1 + coeff_plane * (valid1 * d1_plane);
            score2 = d2 + coeff_plane * (valid2 * d2_plane);
    
            if score1 <= score2
                leafSeeds{s1} = [leafSeeds{s1}; unSegment(idx)];
            else
                leafSeeds{s2} = [leafSeeds{s2}; unSegment(idx)];
            end
    
        else
            % 情况 B：两者都不可信 → junction 区
            junctionIdx = [junctionIdx; unSegment(idx)];        %#ok<AGROW>
            continue;
        end
    end
    
    Regions = leafSeeds';
    
    if DebugShow
        fprintf('[FineSegmentation] Junction points: %d\n', length(junctionIdx));
    end

end

%% 投影距离法——辅助函数
% 核心思想：点到主方向构成的叶面距离，如果某点正好落在叶面上，说明是属于该叶片的概率高
% 计算给定点 pt 到一个由其邻域点拟合出的局部平面的距离
function [projd, isValid] = caldis2(points, pt, SampleNum)
    % 点到叶片局部主平面的投影距离
    % 同时返回该平面是否可信    
    isValid = false;
    projd   = inf;
    
    %% ========= 1. 最近邻采样 =========
    vs  = points - pt;
    rad = sqrt(sum(vs .* vs, 2));
    [~, idx] = sort(rad);
    
    num = min(length(idx), SampleNum);
    samplePts = points(idx(1:num), :);
    
    if size(samplePts, 1) < 5
        return;
    end
    
    %% ========= 2. PCA =========
    [coeff, ~, latent] = pca(samplePts);
    
    if length(latent) < 3
        return;
    end
    
    %% ========= 3. 平面可信度判断 =========
    % 叶片应满足：第三特征值显著小
    % latent(1) >= latent(2) >= latent(3)
    planarity = latent(3) / sum(latent);
    
    % 阈值经验范围：0.12 ~ 0.2（交界区会明显偏大），建议 0.15 ~ 0.2
    if planarity > 0.18
        return;   % 不认为是可靠平面
    end
    
    %% ========= 4. 投影距离 =========
    normal = coeff(:, 3)';
    Spt    = mean(samplePts, 1);
    
    pjs = find_ProjCoord_mt(normal, Spt, pt);
    projd = norm(pjs - pt);
    
    isValid = true;
end


% 最近邻平均距离: 局部密度或稀疏度指标
function [min_dis] = caldis1(points, pt, SampleNum)
% 计算点集 points 到指定点 pt 的“加权欧氏距离”，取最近的 SampleNum 个距离做平均并返回
    vs = points-pt;
%     rad = sqrt(sum(vs .* vs, 2));   % 欧氏距离
    rad = abs(vs(:,3)) + 1.5 * abs(vs(:, 1)) + 1.5 * abs(vs(:,2));    % 加权曼哈顿距离
    [sortr, ~] = sort(rad);     % 距离从小到大排序
    num = size(sortr, 1);

    if(num > SampleNum)
        min_dis = mean(sortr(1:SampleNum));
    else
        min_dis = mean(sortr(1:num));
    end
end