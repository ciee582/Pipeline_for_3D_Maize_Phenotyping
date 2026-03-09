function OBB = find_OBB_mt(points, dir1Num, dir2Num, dir3Num, show_result)
% FIND_OBB_MT 计算三维点云的有向包围盒（OBB）及其体素分割
% 输入：
%   points      - N×3 点云数据
%   dir1Num~3Num- 在主轴方向上的体素划分数
%   show_result - 是否可视化显示OBB结果
% 输出：
%   OBB         - 包含方向、尺寸、体素划分等信息的结构体

%% --- 主方向计算（PCA） ---
[coeff, ~, score] = pca(points);   % coeff: 3x3, 每列为一个主方向
center = mean(points, 1);          % 1x3

n = size(points, 1);
displace = points - center;        % N×3

% 直接用矩阵乘法得到在主轴上的投影
V1 = displace * coeff(:,1);        % N×1
V2 = displace * coeff(:,2);        % N×1
V3 = displace * coeff(:,3);        % N×1

minX = min(V1); maxX = max(V1);
minY = min(V2); maxY = max(V2);
minZ = min(V3); maxZ = max(V3);

%% --- 构造 OBB 基本信息（在 PCA 坐标系下的盒子） ---
OBB.ExtentLen = [maxX - minX; maxY - minY; maxZ - minZ];  % 3x1
OBB.Min       = [minX; minY; minZ];
OBB.Max       = [maxX; maxY; maxZ];

% 主方向向量（世界坐标系）
OBB.dir1 = coeff(:,1)';   % 1x3
OBB.dir2 = coeff(:,2)';   % 1x3
OBB.dir3 = coeff(:,3)';   % 1x3

OBB.center = center;      % 1x3

% 额外保留 PCA 信息（更语义化）
OBB.PrincipalAxes = coeff;   % 3x3
OBB.PCAScore      = score;   % Nx3

%% --- 包围盒八个顶点（世界坐标系） ---
d1 = OBB.dir1;
d2 = OBB.dir2;
d3 = OBB.dir3;

p0 = minX*d1 + minY*d2 + minZ*d3 + center;
p1 = maxX*d1 + minY*d2 + minZ*d3 + center;
p2 = maxX*d1 + maxY*d2 + minZ*d3 + center;
p3 = minX*d1 + maxY*d2 + minZ*d3 + center;
p4 = minX*d1 + minY*d2 + maxZ*d3 + center;
p5 = maxX*d1 + minY*d2 + maxZ*d3 + center;
p6 = maxX*d1 + maxY*d2 + maxZ*d3 + center;
p7 = minX*d1 + maxY*d2 + maxZ*d3 + center;

OBB.p0 = p0; OBB.p1 = p1; OBB.p2 = p2; OBB.p3 = p3;
OBB.p4 = p4; OBB.p5 = p5; OBB.p6 = p6; OBB.p7 = p7;

%% --- 初始化体素分割结构体（在 PCA 坐标系下划分） ---
stepX = OBB.ExtentLen(1) / dir1Num;
stepY = OBB.ExtentLen(2) / dir2Num;
stepZ = OBB.ExtentLen(3) / dir3Num;

template.Min     = [0; 0; 0];
template.Max     = [0; 0; 0];
template.Num     = 0;
template.Indices = [];

OBB.partition = repmat(template, [dir1Num, dir2Num, dir3Num]);

% 填充每个体素单元边界（在 PCA 坐标系下的 min/max）
for i = 1:dir1Num
    for j = 1:dir2Num
        for k = 1:dir3Num
            OBB.partition(i,j,k).Min = OBB.Min + [stepX*(i-1); stepY*(j-1); stepZ*(k-1)];
            OBB.partition(i,j,k).Max = OBB.Min + [stepX*i;       stepY*j;       stepZ*k];
        end
    end
end

%% --- 每个点映射到体素单元中 ---
Vpar1 = min(floor((V1 - OBB.Min(1)) / stepX) + 1, dir1Num);
Vpar2 = min(floor((V2 - OBB.Min(2)) / stepY) + 1, dir2Num);
Vpar3 = min(floor((V3 - OBB.Min(3)) / stepZ) + 1, dir3Num);

for i = 1:n
    I1 = Vpar1(i);
    I2 = Vpar2(i);
    I3 = Vpar3(i);
    OBB.partition(I1, I2, I3).Num = OBB.partition(I1, I2, I3).Num + 1;
    OBB.partition(I1, I2, I3).Indices(end+1) = i;
end

OBB.HasPointPartNum = sum(arrayfun(@(p) p.Num > 0, OBB.partition), 'all');
OBB.dir1_ParNum = dir1Num;
OBB.dir2_ParNum = dir2Num;
OBB.dir3_ParNum = dir3Num;

%% --- 可视化结果（可选） ---
if show_result
    figure; hold on;
    scatter3(points(:,1), points(:,2), points(:,3), 5, [0.3 0.3 0.3], 'filled');
    for i = 1:dir1Num
        for j = 1:dir2Num
            for k = 1:dir3Num
                part = OBB.partition(i,j,k);
                if part.Num > 0
                    color = 'r';
                else
                    color = 'g';
                end
                draw_obb_box(part.Min, part.Max, d1, d2, d3, center, color);
            end
        end
    end
    axis equal; view(3); title('OBB + Voxel Partition');
end
end

function draw_obb_box(minV, maxV, dir1, dir2, dir3, center, color)
% minV, maxV: 3x1 (在 PCA 坐标系下)
% dir1~dir3: 1x3 (世界坐标系)
% center:    1x3

    p0 = minV(1)*dir1 + minV(2)*dir2 + minV(3)*dir3 + center;
    p1 = maxV(1)*dir1 + minV(2)*dir2 + minV(3)*dir3 + center;
    p2 = maxV(1)*dir1 + maxV(2)*dir2 + minV(3)*dir3 + center;
    p3 = minV(1)*dir1 + maxV(2)*dir2 + minV(3)*dir3 + center;
    p4 = minV(1)*dir1 + minV(2)*dir2 + maxV(3)*dir3 + center;
    p5 = maxV(1)*dir1 + minV(2)*dir2 + maxV(3)*dir3 + center;
    p6 = maxV(1)*dir1 + maxV(2)*dir2 + maxV(3)*dir3 + center;
    p7 = minV(1)*dir1 + maxV(2)*dir2 + maxV(3)*dir3 + center;

    edges = [p0; p1; p1; p2; p2; p3; p3; p0;
             p4; p5; p5; p6; p6; p7; p7; p4;
             p0; p4; p1; p5; p2; p6; p3; p7];
    for i = 1:2:size(edges,1)
        plot3(edges(i:i+1,1), edges(i:i+1,2), edges(i:i+1,3), color, 'LineWidth', 1);
    end
end
