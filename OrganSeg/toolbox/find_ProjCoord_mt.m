function [ PCoord ] = find_ProjCoord_mt(n, c, v)
    % 计算点 V 到面 P 的投影点坐标
    % %所有输入都是列向量或行向量，且长度相等（3维）。
    
    % 平面的点法式方程为： n * (x - c) = 0 ===> n * x = n * c
    % 点 v 到平面 P 的投影点 p' 位于从 v 沿着法向量 n 方向出发的射线上，并且 p' 必须满足平面方程
    % 投影点 p' 可表示为： p' = v + t * n, 其中 t 为标量，表示沿着法向量方向的位移
    % 所以关键在于求 t
    
    % n · c = dot(n, c)
    % n · n = dot(n, n) = ||n||^2   % 法向量模长的平方
    
    % 1. 计算 t = (n · c - n · v) / (n · n)
    numerator = dot(n, c) - dot(n, v);
    denominator = dot(n, n);
    
    % 检查法向量是否为零向量，防止除以零
    if denominator == 0
        error('法向量 dir1 不能为零向量。');
    end
    
    t = numerator / denominator;
    
    % 2. 计算投影点 p' = v + t * n
    % v + n .* t 相当于 v + t * n
    PCoord = v + n .* t;

end

