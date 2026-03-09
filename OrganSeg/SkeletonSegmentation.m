%% Init_Env
clc;
clear;

addpath('toolbox');
addpath('skeleton');
addpath('segment');
addpath('Visuallization');

%% Load PCD data
% pathname = 'E:\Datasets\3d\20241012\seg\';
% pathname = 'E:\Datasets\3d\20241012\plants240\';
% filename = '08-2.ply';
% filename = '11-3.ply';
% filename = '06-1_GT.txt';
% pathname = 'E:\Datasets\3d\20241012\2-2\';
% filename = '2-2.ply';

% pathname = 'E:\Datasets\3d\20241026\segNew\';
% filename = '09-4a.ply';  % "L"型
% filename = '07-3a.ply';  % "丄"型
% filename = '13-1a.ply';
pathname = 'E:\Datasets\3d\20241026\plants222\';
% filename = '15-4.ply';    % 有下垂断叶
filename = '10-4.ply';
% filename = '03-4_GT.txt';

% pathname = 'E:\Datasets\3d\20241214\plants-Octree104\';
% filename = '13-4.ply';
% filename = '06-3_GT.txt';
% filename = '07-4a-Octree7.ply';
% filename = '15-1.ply';

% pathname = 'E:\Datasets\3d\pheno4d\';
% filename = 'M01_0324_a-Octree8.ply';

ptCloud = pcread([pathname filename]);
points = ptCloud.Location;  % Nx3 点坐标
% points = pointsDataLoad(pathname, filename, 'ply', [1 2 3 4]);  % Nx3 点坐标
% points = pointsDataLoad(pathname, filename, 'txt', [1 2 3 4]);  % Nx3 点坐标

% ptCloud = load([pathname filename]);
% points = ptCloud(:, 1:3);  % Nx3 点坐标

g_P=[];
% g_regions=[];

g_P.pts = double(points);

g_P.pathName = pathname;        % 为中间缓存数据提供相关信息
g_P.fileName = filename;

ShowPoints(points, '');

%% Obtain Skeleton
Parameters.t1 = 0.8;            % 内部分支 距离 阈值
Parameters.a1 = pi*5.0/7.0;     % 内部分支 角度 阈值
Parameters.t2 = 0.1;            % 去除无关极值阈值
Parameters.t3 = 3;              % 移除小环路阈值
Parameters.KnnNum = 10;          % 12 older 11/10 edit
Parameters.sampleScale = 0.02;  % 采样比例，越小骨架包含的点数越多
Parameters.smoothSkeleton = 0;  % 是否对骨架进行平滑处理

Parameters.saveSkeleton = 0;

g_P = skeleton_laplacian(g_P, Parameters);
% g_P.spls        骨架点云，代表骨架结构的关键点（n * 3）
% g_P.corresp     原始点云（pts）中的每个点与采样后的骨架点（spls）之间的对应关系,一个 n x 1 的数组
% g_P.spls_adj    骨架点的邻接矩阵，表示骨架点之间的连接关系
% g_P.joints      骨架拓扑结构中的关键节点，如分支点或连接点
% g_P.roots       骨架的端点，骨架的末段点或拓扑结构的参考点      
% g_P.branches    骨架的分支结构，即骨架中从一个关节点到另一个关节点（或端点）的路径段

%% Show Result
figure('Name','Find joints','NumberTitle','off');
set(gcf,'color','white'); movegui('southwest');
hold on;

scatter3(g_P.spls(:,1),g_P.spls(:,2),g_P.spls(:,3), 2, [0,0,0], 'filled');
% for i = 1:length(g_P.spls)
%     if ismember(i, g_P.joints)
%         text(g_P.spls(i,1),g_P.spls(i,2),g_P.spls(i,3), num2str(i), 'Color', 'black', 'FontSize', 18);
%     end
% end

plot_skeleton(g_P.spls, g_P.spls_adj);

scatter3(g_P.spls(g_P.joints,1),g_P.spls(g_P.joints,2),g_P.spls(g_P.joints,3), 60, [1,0,0], 'filled');

scatter3(g_P.spls(g_P.roots,1),g_P.spls(g_P.roots,2),g_P.spls(g_P.roots,3), 100, [1,0,1], 'filled');

scatter3(g_P.spls(g_P.branches,1),g_P.spls(g_P.branches,2),g_P.spls(g_P.branches,3), 10, [0,1,1], 'filled');

axis auto; axis off; axis equal; camorbit(0, 0, 'camera'); view(-180,0); view3d rot;

%% Decompose Skeleton
% 每个点周围 K 个最近邻点的分布特性
% 最终特征 features(i) 是 e3 / e2，即第三主成分特征值与第二主成分特征值的比值
% 当 e3 / e2 接近 1：e3 和 e2 大小相近，局部点云更接近于平面状（2D 分布），如一个扁平的圆盘或薄片。
% 当 e3 / e2 很小：e3 远小于 e2，局部点云更接近于线状（1D 分布），如一条细长的曲线。
% 当 e3 / e2 介于两者之间：表示点云的形状介于线状和面状之间，可能是一个拉长的椭圆形分布。
g_P.features = computerFeature_mt(g_P.pts(:, 1:3), 12, false);

% 基于广度优先搜索（BFS）的函数，将输入的无向图（通过邻接矩阵 spls_adj 表示）分解为多个连通分量
% 再从多个连通分量中选择节点数量最多的连通分量作为主要骨架
% 主要目的是，去除骨架中的孤立点
g_P.mainskeleton = FindMainSkeleton(g_P.pts, g_P.spls, g_P.joints, g_P.roots, g_P.spls_adj, g_P.corresp);
fprintf('骨架点数: %d\n', size(g_P.spls, 1));
fprintf('去除孤立点的主骨架点数: %d\n', size(g_P.mainskeleton, 1));

%% 器官类型标记（1=叶片，2=主茎）
% sub_skeletons：一个细胞数组（cell array），每个元素是一个骨架子集的点索引列表。
% sub_skeletons{end} 通常是主要骨架，而 sub_skeletons{1:end-1} 是次要骨架（其他连通分量）。
[sub_skeletons, organType, g_P.joints] = SkeletonDecomposition2(g_P.mainskeleton, g_P.spls, g_P.joints, g_P.roots, g_P.spls_adj, g_P.corresp, g_P.features, true);

%% 基于骨架的粗分割
% 主要利用了 corresp 骨架与点云映射关系信息，将点云归属到相应骨架点所属类别中
% 器官骨架点集（Phi_O）：sub_skeletons{1:end-1} 对应的点云点，分为多个子集
% 未分配点集（Phi_U）：不属于任何骨架的点云点
% 茎秆骨架点集（Phi_S）：sub_skeletons{end} 对应的点云点
[Phi_O, Phi_U, Phi_S] = CoarseSegBySkeleton(g_P.pts, g_P.spls, g_P.corresp, sub_skeletons, g_P.joints, 3, true);

% 注意：此时 Phi_U 中，包含 Phi_S，即未分配点中包含茎秆点

%% transform global coordinate axis to plant coordinate axis
[PA_Pts, PA_Spls] = ConstructPlantAxis(g_P.pts, g_P.spls, Phi_S, sub_skeletons);

%
% PA_StopX  对应每个叶片连接点（joint）在主茎骨架往上游走两步的点的 Z 值
% PA_StartX 对应每个叶片连接点（joint）本身的 Z 坐标
[PA_StopX, PA_StartX] = FindLeafStopZ(sub_skeletons, g_P.joints, PA_Spls);

%
% [Phi_S,Phi_U]=StemSegment(PA_Pts,Phi_U,PA_StartX,alpha,JN);
%% StemSegment2 最后一个参数表示茎秆的高度：0-从根节点到最高的连接点；1-最高连接点向下偏移一个连接点；2-向下偏移2个连接点。。。
% 这里可不执行
[Phi_S, Phi_U] = StemSegment2(PA_Pts, Phi_U, PA_Spls, Phi_S, 1.1, PA_StartX, 0);

%% 可视化
figure('Name','Distance constraint','NumberTitle','off');set(gcf,'color','white');
% scatter3(PA_Pts(:,1),PA_Pts(:,2),PA_Pts(:,3),5,[0.0 0.0 0], 'filled');
hold on
scatter3(PA_Pts(Phi_S,1),PA_Pts(Phi_S,2),PA_Pts(Phi_S,3),50,[1 1 0], 'filled');
scatter3(PA_Pts(Phi_U,1),PA_Pts(Phi_U,2),PA_Pts(Phi_U,3),5,[0 0 0], 'filled');
% scatter3(PA_Pts(Phi_O,1),PA_Pts(Phi_O,2),PA_Pts(Phi_O,3),5,[0 1 1], 'filled');
axis on; axis equal; hold off;

%% Phi_O = flip(Phi_O);
Phi_O{end+1} = Phi_S;
PA_StopX(end+1) = -inf;
PA_StopX = flip(PA_StopX);
Phi_O = flip(Phi_O);
sub_skeletons = flip(sub_skeletons);

g_P.PA_Pts = PA_Pts;
g_P.PA_Spls = PA_Spls;
g_P.Phi_U = Phi_U;
g_P.Phi_O = Phi_O;
g_P.sub_skeletons = sub_skeletons;
g_P.PA_StopX = PA_StopX;

% ShowPoints(g_P.PA_Pts, Phi_O);

%% Segmentation Organs
id = g_P.sub_skeletons{1}(1);       % 茎秆最低点,因为进行了 flip 反转操作，所以这里为 sub_skeletons{1} 不是 sub_skeletons{end}
StemBottomPt = g_P.PA_Spls(id, :);
k1 = 20;     % 值越大，叶片吸附的茎秆点越多，典型例子E:\Datasets\3d\20241026\plants222\05-4.ply
% k1 = 3;     % 值越大，叶片吸附的茎秆点越多，典型例子E:\Datasets\3d\20241026\plants222\10-3a.ply
k2 = 5;
Beta = 0;
Sigma = 1;
g_regions = FineSegmentation(g_P.PA_Pts, StemBottomPt, g_P.Phi_O, g_P.Phi_U, g_P.PA_StopX, k1, k2, Beta, Sigma, false);
% 第一个元素为茎秆，后续为叶片

% --------------------
% 可视化
% --------------------
figure; hold on;
numOrgans = length(g_regions);
color = MyGS.MYCOLOR;

visPts = g_P.PA_Pts;

for i = 1:numOrgans
    organPointsIdx = g_regions{i};  % 对应骨架路径
    if isempty(organPointsIdx)
        continue;
    end

    organPoints = visPts(organPointsIdx, :);

%     [~, fName, fExt] = fileparts(g_P.fileName);
%     organFile = fullfile(g_P.pathName, [fName '-' num2str(i) fExt]);
%     pcwrite(pointCloud(organPoints), organFile, 'Encoding','binary');
%     fprintf('植株 %s 的第 %d 个器官，保存至 %s（共 %d 点）\n', fName, i, organFile, size(organPoints, 1));

    basePoint = organPoints(1, :);
    
    if i == 1
        % 主茎序号1，标记在起点，稍作偏移
        text(organPoints(end, 1)+0.05, organPoints(end, 2)-0.01, organPoints(end, 3)-0.2, '1', 'Color', 'black', 'FontSize', 18, 'FontWeight', 'bold');
        scatter3(organPoints(:, 1), organPoints(:, 2), organPoints(:, 3), 10, color(end, :), 'filled');
    else
        % 叶片编号从2开始，标记在骨架中点位置
        scatter3(organPoints(:, 1), organPoints(:, 2), organPoints(:, 3), 10, color(i, :), 'filled');
        text(basePoint(1)-0.05, basePoint(2)-0.01, basePoint(3)+0.1, num2str(i), 'Color', 'black', 'FontSize', 18);
    end

%     pause(1);
end
axis off; axis equal; grid off; view3d rot;
title('植株点云分类结果');

%% Phenotyping
addpath('PhenotypicTrait');

filter = [filename '_traits.txt'];

traitfile = [pathname filter];

PhenotypicTrait(g_P.pts, g_regions, traitfile, false);

disp('Processing done');


%% 辅助函数
function pts = pointsDataLoad(pathname, filename, loadType, validClass)
    if nargin < 4, validClass = [1 2]; end
    filepath = fullfile(pathname, filename);

    if strcmp(loadType, 'txt')
        segData = load(filepath);
        if size(segData, 2) < 4
            error('文件列数不足4列');
        end
        segXYZ = segData(:, 1:3);
        segC   = round(segData(:, 4));
        valid  = ismember(segC, validClass);
        pts    = segXYZ(valid, :);          % 保证 M×3
    else
        ptCloud = pcread(filepath);
        if isempty(ptCloud.Location)
            error('点云文件无有效坐标');
        end
        pts = ptCloud.Location;             % N×3
    end
end