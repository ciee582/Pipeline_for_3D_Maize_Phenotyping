function [ StemRegions, UnSegmentRegions] = StemSegment2( points, unknowns, spls, stem_skeleton, beta, PA_StartX, JN)
    %                        [Phi_S2, Phi_U2] = StemSegment2( PA_Pts, Phi_U,   PA_Spls, sub_skeletons{end}, 1, PA_StartX, 2);
    % 将点云分为树干区域（StemRegions）和非树干区域（UnSegmentRegions）
    % unknowns ==> Phi_U : 不属于任何骨架的点云点
    % stem_skeleton ==> sub_skeletons{end} : 茎秆骨架点
    
    % 提取所有未分配点（可能包含：茎秆点及其对应原始点、因中断等未被划分到任何叶片的点）
    unAllocPts = points(unknowns, :);
    unAllocPts_Z = unAllocPts(:, 3);
    unAllocPts_mz = median(unAllocPts_Z);
    
    % 低于总高度 20% 的点被视为候选茎秆点
    bottomIdx = find( unAllocPts_Z < unAllocPts_mz / 20);
    StemPts = unAllocPts(bottomIdx, :);
    
    % 拟合底部点形成的直线，估计茎秆方向，作为主轴
    [center, dir] = fitline(StemPts);
    
    % 计算所有未分配点到拟合直线的距离
    total=[];
    for i = 1:size(unAllocPts, 1)
        pt = unAllocPts(i, :);
        [dis, ~ ] = P2LineDistance(center, dir, pt);
        total = [total; dis];
    end
    
    % 利用底部点到主轴的距离估计半径
    SearchR = beta * median(total(bottomIdx));
    
    %%%%%% 沿茎秆骨架，对骨架点邻域内进行球形搜索 %%%%%%
    % 利用主茎骨架点为中心，以 SearchR 半径找周围点
    % 聚合这些邻域点，组成候选的茎秆点集
    StemRegions = [];
    for i = 1:length(stem_skeleton)
        sp = points(stem_skeleton(i), :);      % 当前骨架点坐标
        ptCloud = pointCloud(unAllocPts);     % 换为 MATLAB 的 pointCloud 对象，用于高效邻域搜索
        [indices, ~] = findNeighborsInRadius(ptCloud, sp, SearchR);     % 找到 ptCloud 中与 sp 距离小于 SearchR 的点的索引
        StemRegions = [StemRegions; unknowns(indices)];
    end
    StemRegions = unique(StemRegions);      % 确保树干点索引不重复（因为不同骨架点的邻域可能重叠）
    
    % 对提取出的茎秆区域进行顶部截断（防止将连接点或部分叶片误识为茎）
    temp = sort(unique(PA_StartX)); % 各子骨架最低点进行由低到高的排序（升序）
    StopZ = temp(end);              % 最高的子骨架末端点作为茎秆的分割上限
    if(length(temp) > JN)
        StopZ = temp(end - JN);     % 如果子骨架多于 JN ，则分割上限往下降低一些
    end
    
    StemZ = points(StemRegions, 3);
    % uindex = find(StemZ < StopZ);      
    StemRegions = StemRegions(StemZ < StopZ);   % 仅保留 Z 值低于 StopZ 的点作为最终 Stem
    UnSegmentRegions = setdiff(unknowns, StemRegions);

end

