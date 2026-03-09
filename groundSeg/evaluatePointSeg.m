%% Init_Env
clc;
clear;

%% 读取真值和分割结果的PLY点云文件
FileStruct(1) = struct('basePath', 'E:\Datasets\3d\20241012\',...
                    'Ours', 'plants241-withLables.ply', ...
                    'CSF', 'plants240-withLables-CSF2.ply', ...
                    'SMRF', 'plants242-withLables-SMRF.ply', ...
                    'GT', 'maize-60-Cloud-SOR-plants222-optimized-Oct11-Labeled-GT - Cloud.ply', ...
                    'Stage', 'V6');% 真值文件路径;

% segmentedFile = 'plants240-withLables-CSF2.ply';
% segmentedFile = 'plants242-SMRF-withLables.ply';

% groundTruthFile = 'plants241-withLables-Single-GT.ply';
% segmentedFile = 'plants241-withLables_wiithSingleLabeles3.ply';

% folder = [pathname 'plants240\'];
% ext    = '*.txt';
% % pattern = '^\d{1,2}-\d\.ply$';  % 这里 \. 表示字面量点号
% gtPattern = '^\d{1,2}-\d_GT\.txt$';
% segPattern = '^\d{1,2}-\d_XYZC\.txt$';
% 
% gtFiles = findFilesByPattern(folder, ext, gtPattern);
% 
% if isempty(gtFiles)
%     disp('未找到匹配文件');
% end
% 
% segFiles = findFilesByPattern(folder, ext, segPattern);
% 
% if isempty(segFiles)
%     disp('未找到匹配文件');
% end

% pathname = 'E:\Datasets\3d\20241026\';
% groundTruthFile = '目标42.las.segmented-Octree11-2-GT.ply'; % 替换为真值文件路径
% segmentedFile = '目标42.las.plants2415-labels.ply'; % 替换为分割结果文件路径
% segmentedFile = '目标42.las.plants2415-labels-CSF.ply';
% segmentedFile = '目标42.las.plants2415-labels-SMRF.ply';
% segmentedFile = '目标42.las.segmented-Octree11-2-withPlantsLables-0.2.ply';
% segmentedFile = '目标42.las.segmented-Octree11-2-CSF.ply';
% segmentedFile = '目标42.las.segmented-Octree11-2-SMRF.ply';
FileStruct(2) = struct('basePath', 'E:\Datasets\3d\20241026\', ...
                    'Ours', '目标42.las.segmented-Octree11-2-withPlantsLables-0.4.ply', ...
                    'CSF', '目标42.las.plants2415-labels-CSF.ply', ...
                    'SMRF', '目标42.las.segmented-Octree11-2-SMRF.ply', ...
                    'GT', '目标42.las.segmented-Octree11-2-GT.ply', ...
                    'Stage', 'V8');% 真值文件路径;

% groundTruthFile = '目标42.las.plants2415-labels_SingleLabeles_GT.ply';
% segmentedFile = '目标42.las.plants2415-labels_wiithSingleLabeles.ply';

% pathname = 'E:\Datasets\3d\20241214\';
% groundTruthFile = '目标4-Cloud-ClassOctree10-GT-labels01-align.ply'; % 替换为真值文件路径
% segmentedFile = '目标4-Cloud-plants-01-withLables.ply'; % 替换为分割结果文件路径
% segmentedFile = '目标4-Cloud-CSF-withLables.ply';
% segmentedFile = '目标4-Cloud-SMRF-withLables.ply';
FileStruct(3) = struct('basePath', 'E:\Datasets\3d\20241214\', ...
                    'Ours', '目标4-Cloud-plants-01-withLables.ply', ...
                    'CSF', '目标4-Cloud-CSF-withLables.ply', ...
                    'SMRF', '目标4-Cloud-SMRF-withLables.ply', ...
                    'GT', '目标4-Cloud-ClassOctree10-GT-labels01-align.ply', ...
                    'Stage', 'V11');% 真值文件路径;

% groundTruthFile = '目标4-Cloud-plants-01-withLables-SOR-segment-GT-intensity.ply';
% segmentedFile = '目标4-Cloud-plants-01-withLables-SOR-segment-GT_wiithSingleLabeles.ply';

%% 分时期、分方法进行分割效果指标计算
stages = {'V6', 'V8', 'V11'};
methods = {'CSF','SMRF','Ours'}; % 与字段同名

for i = 1:numel(FileStruct)
    resNum = 0;

    pc_gt = pcread(fullfile(FileStruct(i).basePath, FileStruct(i).GT));
    
    for j = 1:numel(methods)
        pc_seg = pcread(fullfile(FileStruct(i).basePath, FileStruct(i).(methods{j})));
        
        Metrics = evaluateSegmentation_Enhanced(pc_gt, pc_seg);
        
        %% 转为 table
        tmpTbl = struct2table(Metrics);

        % 追加文件序号 & 方法名
        tmpTbl.Stage = stages{i};
        tmpTbl.Method = methods{j};
        
        saveTableFile = 'E:\Datasets\3d\result\groundSegResult.csv';
        
        % --- 第一次写入（带表头）---
        if ~isfile(saveTableFile)
            writetable(tmpTbl, saveTableFile, 'Delimiter', ',');  % 自动写入表头
        else
            % --- 后续追加（不写表头）---
            writecell(table2cell(tmpTbl), saveTableFile, 'WriteMode', 'append', 'Delimiter', ',');
        end
    
        resNum = resNum + 1;
    end

    fprintf('===> Stage: %s APPEND %d Result Datasets which are saved to: %s\n', stages{i}, resNum, saveTableFile);
end

%% 单独计算
stages = {'V6', 'V8', 'V11'};
methods = {'CSF','SMRF','Ours'}; % 与字段同名
i = 2;
j = 3;
pc_gt = pcread(fullfile(FileStruct(i).basePath, FileStruct(i).GT));
pc_seg = pcread(fullfile(FileStruct(i).basePath, FileStruct(i).(methods{j})));
% pc_seg = pcread(fullfile(FileStruct(i).basePath, 'maize-60-Cloud-SOR-plants222-optimized-Oct11-withPlantsLables-0.2.ply'));

Metrics = evaluateSegmentation_Enhanced(pc_gt, pc_seg);
disp(Metrics);

%% 转为 table
tmpTbl = struct2table(Metrics);

% 追加文件序号 & 方法名
tmpTbl.Stage = stages{i};
tmpTbl.Method = methods{j};

saveTableFile = 'E:\Datasets\3d\result\groundSegResult.csv';

% --- 第一次写入（带表头）---
if ~isfile(saveTableFile)
    writetable(tmpTbl, saveTableFile, 'Delimiter', ',');  % 自动写入表头
else
    % --- 后续追加（不写表头）---
    writecell(table2cell(tmpTbl), saveTableFile, 'WriteMode', 'append', 'Delimiter', ',');
end

fprintf('===> Stage: %s APPEND 1 Result Datasets which are saved to: %s\n', stages{i}, saveTableFile);

%% 计算体素平均偏差(VD_avg)、点数平均偏差PCD_avg
% metrics = evaluateSegmentation(pc_gt, pc_seg);
voxel_size = 0.0025;
metrics = evaluateSegForMultiClass(pc_gt, pc_seg, voxel_size);

% metrics.Name = 'Ours';
% metrics.Name = 'CSF';
% metrics.Name = 'MSRF';
%% 单株效果评价
sigMetrics = evaluteSinglePlants(pc_gt, pc_seg);

disp(sigMetrics);

%% Visualize
% 创建颜色数组：基于 Intensity 值着色
numPoints = pc_seg.Count;  % 点云总数
colors = zeros(numPoints, 3);  % 初始化颜色数组 (RGB)

% Intensity == 0 的点设为灰色 (背景)
bgColor = [0.5 0.5 0.5];  % 灰色
bgIdx = (pc_seg.Intensity == 0);
colors(bgIdx, :) = repmat(bgColor, sum(bgIdx), 1);

% Intensity == 1 的点设为红色 (前景)
fgColor = [0 1 0];  % 红色
fgIdx = (pc_seg.Intensity == 1);
colors(fgIdx, :) = repmat(fgColor, sum(fgIdx), 1);

points = pc_seg.Location;

figure('Position',[200 150 900 600]);
axNew = gca;
scatter3(points(:, 1), points(:, 2), points(:, 3), 3, colors, "filled");   % 地面（绿色）
% pcshow(pc_seg.Location, colors, 'MarkerSize', 50);  % MarkerSize 可调整点的大小
xlabel('X');
ylabel('Y');
zlabel('Z');
title('基于 Intensity 的点云可视化 (0: 灰色背景, 1: 绿色前景)');
grid off;
% for i = 1:numClusters    
%     scatter(clusterCenters(i, 1), clusterCenters(i, 2), 50, 'b', 'x', 'LineWidth', 2); % 绘制聚类中心（红色X）
%     % 显示标签信息
% %     text(clusterCenters(i, 1), clusterCenters(i, 2), int2str(i), 'Color', 'red', 'FontSize', 14);
% end

% 载入参数并应用到新坐标区
load('myView2.mat','viewParam');   % 结构体 viewParam 回到工作区

axNew.CameraPosition  = viewParam.CameraPosition;
axNew.CameraTarget    = viewParam.CameraTarget;
axNew.CameraUpVector  = viewParam.CameraUpVector;
axNew.XLim = viewParam.XLim;
axNew.YLim = viewParam.YLim;
axNew.ZLim = viewParam.ZLim;


%======================================================================
%% 辅助函数
%======================================================================
function metrics = evaluateSegmentation(points_gt, points_seg)
%==========================================================================
%  evaluateSegmentation
%
%  使用最近邻标签转移（Nearest Neighbor Label Transfer）方法评估
%  分割结果与 GT 之间的一致性。
%
%  输入:
%      points_gt  - pcread 读取的 GT 点云 pointCloud
%      points_seg - pcread 读取的 分割结果 点云 pointCloud
%
%  假设:
%      labels 由 pointCloud.Color(:,1)（R 通道）表示：0 或 1
%      如不是，可根据需要修改 getLabels() 函数
%
%  输出:
%      metrics (struct):
%          TP, FP, FN, TN
%          Precision, Recall, IoU, F1
%
%==========================================================================
    %------------------------%
    % 1. 获取点与标签
    %------------------------%
    pts_gt  = points_gt.Location;
    pts_seg = points_seg.Location;

    labels_gt  = getLabels(points_gt);
    labels_seg = getLabels(points_seg);

    if isempty(pts_gt) || isempty(pts_seg)
        error('GT 或 Seg 点云为空！');
    end

    % 提取两个类别的点
    plant_gt = pts_gt(labels_gt==1, :);
    plant_seg = pts_seg(labels_seg==1, :);
    
    ground_gt = pts_gt(labels_gt==0, :);
    ground_seg = pts_seg(labels_seg==0, :);

    %% =======================
    % 5. Boundary-F1（边界敏感指标）
    % ========================
    voxel_size = 0.011;
    BFScore = computeBoundaryF1(plant_gt, plant_seg, voxel_size);

    %% =======================
    % 6. Volume Deviation (VD) - 体积偏差指标
    % ========================
    shp_gt = alphaShape(double(plant_gt));         % 创建alpha形状（自动处理边界）
    vol_plant_gt = volume(shp_gt);          % 计算体积
    
    shp_seg = alphaShape(double(plant_seg));
    vol_plant_seg = volume(shp_seg);
    
    VD_Plant = abs(vol_plant_seg - vol_plant_gt) / vol_plant_gt;

    %% =======================
    % 7. 点数偏差 PCD（Ground & Plant）
    % ========================
    
    PCD_ground = abs(size(ground_seg,1) - size(ground_gt,1)) / size(ground_gt,1) * 100;
    PCD_plant  = abs(size(plant_seg,1) - size(plant_gt,1)) / size(plant_gt,1) * 100;
    PCD_avg    = (PCD_ground + PCD_plant)/2;

    %------------------------%
    % 5. 返回结果
    %------------------------%
    metrics = struct();
    metrics.BFScore = BFScore;
    metrics.VD_Plant = VD_Plant;
    metrics.PCD_avg = PCD_avg;

end

function metrics = evaluateSegForMultiClass(points_gt, points_seg, voxel_size)
%------------------------
% 1. 获取点与标签
    pts_gt = points_gt.Location;
    pts_seg = points_seg.Location;
    labels_gt = getLabels(points_gt);
    labels_seg = getLabels(points_seg);
if isempty(pts_gt) || isempty(pts_seg)
        error('GT 或 Seg 点云为空！');
end

% 获取唯一的标签
unique_labels = unique(labels_gt);

% 初始化存储每个类别的 VD 和 PCD
VD_values = zeros(length(unique_labels), 1);
PCD_values = zeros(length(unique_labels), 1);

for i = 1:length(unique_labels)    
    class_id = unique_labels(i);
    
    % 提取当前类别的点
    class_gt = pts_gt(labels_gt == class_id, :);
    class_seg = pts_seg(labels_seg == class_id, :);

    % ========================
    % Volume Deviation (VD) via Voxelization
    % ========================
    if isempty(class_gt) && isempty(class_seg)
        VSS = 1;
    elseif isempty(class_gt) || isempty(class_seg)
        VSS = 0;
    else
        vox_gt = voxelize2(class_gt, voxel_size);
        vox_seg = voxelize2(class_seg, voxel_size);
        vol_gt = size(vox_gt, 1);
        vol_seg = size(vox_seg, 1);
        
        if vol_gt == 0 && vol_seg == 0
            VSS = 1;
        elseif vol_gt == 0 || vol_seg == 0
            VSS = 0;
        else
            VSS = 2 * min(vol_seg, vol_gt) / (vol_seg + vol_gt);
        end
    end
    VD_values(i) = VSS;
    
    % =======================
    % Point Count Deviation (PCD)
    num_gt = size(class_gt, 1);
    num_seg = size(class_seg, 1);
    if num_gt > 0
        PCD_values(i) = abs(num_seg - num_gt) / num_gt * 100;
    else
        PCD_values(i) = 0;
    end

    fprintf("Class-%d: VD=%.4f, PCD=%.4f\n", i, VD_values(i), PCD_values(i));
end

% 计算平均值
VD_avg = mean(VD_values);
PCD_avg = mean(PCD_values);

%------------------------%
% 5. 返回结果
%------------------------%
    metrics = struct();
    metrics.VD_avg = VD_avg;
    metrics.PCD_avg = PCD_avg;
end

function Metrics = evaluateSegmentation_Enhanced(pcloud_gt, points_seg)
%% =======================
%  1. 读取点云与标签
% ========================
points_GT = double(pcloud_gt.Location);
points_SEG = double(points_seg.Location);

label_gt  = pcloud_gt.Intensity; 
label_seg = points_seg.Intensity;

%% =======================
%  2. 数量对齐（必要）
% =======================
Mdl = KDTreeSearcher(points_GT);
idx_nn = knnsearch(Mdl, points_SEG);   % 每个 seg 点找最近 GT 点

% 由 GT 最近邻获得标签
labels_gt_nn = label_gt(idx_nn);

%% =======================
% 3. 基础混淆矩阵指标
% =======================
TP = sum(labels_gt_nn==1 & label_seg==1);   % 正确识别为植株的植株点
FP = sum(labels_gt_nn==0 & label_seg==1);   % 误识别为植株的地面点
FN = sum(labels_gt_nn==1 & label_seg==0);   % 误识别为地面的植株点
TN = sum(labels_gt_nn==0 & label_seg==0);   % 正确识别为地面的地面点

Precision = TP / (TP + FP + eps);           % 算法识别植株的纯净度，避免把地面误分为植株
Recall_p  = TP / (TP + FN + eps);           % 算法找到植株全部点的能力
Recall_g  = TN / (TN + FP + eps);           % 地面召回率 (背景区分能力),地面召回率低，意味着有大量的地面点被误分为植株点
IoU_p     = TP / (TP + FP + FN + eps);      % 综合评估植株分割效果
IoU_g     = TN / (TN + FP + FN + eps);
F1        = 2*TP / (2*TP + FP + FN + eps);

%% =======================
% 6. Root-Neighborhood Error (RNE)
%    —— 根部区域（Z < 0.3m）
%  V6: 0-0.2
%  V8: 0.05-0.25
% V11: 0.1-0.4
% =======================
z_th = min(points_GT(:,3)) + 0.15;  % 前 30cm 定义为根部区域

% 找到 GT 根部点索引
root_mask_gt = (points_GT(:,3) >= 0.04) & (points_GT(:,3) <= z_th) ;
num_root_gt = sum(root_mask_gt);

% 处理空集情形
if num_root_gt == 0
    % 没有定义的根部点
    ROI_Prec = 0;
    ROI_Rec  = 0;
    ROI_F1   = 0;
else
    root_pts_gt = points_GT(root_mask_gt, :);         % GT 根部点坐标
    true_labels_root_gt = labels_gt_nn(root_mask_gt); % 这些点的真实标签（0/1）

    % 如果 SEG 中没有任何点，直接返回 0（或按需求设为 NaN）
    if isempty(points_SEG)
        pred_labels_for_gt_root = zeros(size(true_labels_root_gt));
    else
        % 用 SEG 点云构建 KD-tree（对大点云也高效）
        MdlSEG = KDTreeSearcher(points_SEG);

        % 对每个 GT 根部点，找在 SEG 中的最近邻点索引
        idx_near = knnsearch(MdlSEG, root_pts_gt, 'K', 1);

        % 用该最近邻点的分割标签作为预测标签
        pred_labels_for_gt_root = label_seg(idx_near);
    end

    % 计算根部区域的混淆量
    TP_root = sum(pred_labels_for_gt_root == 1 & true_labels_root_gt == 1);
    FP_root = sum(pred_labels_for_gt_root == 1 & true_labels_root_gt == 0);
    FN_root = sum(pred_labels_for_gt_root == 0 & true_labels_root_gt == 1);

    % 计算 Precision / Recall / F1 for root region（稳健处理除0）
    ROI_Prec = TP_root / (TP_root + FP_root + eps);
    ROI_Rec  = TP_root / (TP_root + FN_root + eps);
    ROI_F1   = 2 * ROI_Prec * ROI_Rec / (ROI_Prec + ROI_Rec + eps);
end


%% =======================
% 7. 体积偏差 Volume Deviation 和 点数偏差
% =======================
voxel_size = 0.0020;    % 体素分辨率  0.0025
% 获取唯一的标签
unique_labels = unique(labels_gt_nn);

% 初始化存储每个类别的 VD 和 PCD
VD_values = zeros(length(unique_labels), 1);
PCD_values = zeros(length(unique_labels), 1);

for i = 1:length(unique_labels)
    class_id = unique_labels(i);
    
    % 提取当前类别的点
    class_gt = points_GT(labels_gt_nn == class_id, :);
    class_seg = points_SEG(label_seg == class_id, :);

    % ========================
    % Volume Deviation (VD) via Voxelization
    % ========================
    if isempty(class_gt) && isempty(class_seg)
        VSS = 1;
    elseif isempty(class_gt) || isempty(class_seg)
        VSS = 0;
    else
        vox_gt = voxelize2(class_gt, voxel_size);
        vox_seg = voxelize2(class_seg, voxel_size);
        vol_gt = size(vox_gt, 1);
        vol_seg = size(vox_seg, 1);

%         shp_gt = alphaShape(double(class_gt));  % 创建alpha形状（自动处理边界）
%         vol_gt = volume(shp_gt);                % 计算体积
%         shp_seg = alphaShape(double(class_seg));
%         vol_seg = volume(shp_seg);
        
        if vol_gt == 0 && vol_seg == 0
            VSS = 1;
        elseif vol_gt == 0 || vol_seg == 0
            VSS = 0;
        else
            VSS = 2 * min(vol_seg, vol_gt) / (vol_seg + vol_gt);
        end
    end
    VD_values(i) = VSS;
    
    % =======================
    % Point Count Deviation (PCD)
    num_gt = size(class_gt, 1);
    num_seg = size(class_seg, 1);
    if num_gt > 0
        PCD_values(i) = abs(num_seg - num_gt) / num_gt * 100;
    else
        PCD_values(i) = 0;
    end

    fprintf("Class-%d: VD=%.4f, PCD=%.4f\n", i, VD_values(i), PCD_values(i));
end

% 计算平均值
VD_avg = mean(VD_values);
PCD_avg = mean(PCD_values);

%% =======================
% 9. 输出结构体
% =======================
Metrics = struct( ...
    'TP',TP,'FP',FP,'FN',FN,'TN',TN, 'NPts_P_gt', sum(labels_gt_nn), 'NPts_P_seg', sum(label_seg),...
    'Precision',Precision,'Recall_p',Recall_p,'Recall_g',Recall_g,'IoU_p',IoU_p,'IoU_g',IoU_g,'F1',F1, ...
    'VD_g', VD_values(1), 'VD_p', VD_values(2), 'PCD_g', PCD_values(1), 'PCD_p', PCD_values(2), 'VD_avg', VD_avg, 'PCD_avg', PCD_avg,...
    'RNE_F1',ROI_F1,'RNE_Precision',ROI_Prec,'RNE_Recall',ROI_Rec);

end

function BoundaryF1  = computeBoundaryF1(pts_gt, pts_seg, voxel_size)
%======================================================================
% computeBoundaryF1 计算点云分割边界敏感指标 Boundary-F1
%
% 输入：
%   pts_gt: ground truth
%   pts_seg: 分割结果
%   voxel_size: 体素大小（单位与点云一致）
%
% 输出：
%   BoundaryF1: BF1 分数

    % --------------------------
    % 1. 提取目标点
    % --------------------------
    
    % --------------------------
    % 2. 统一 voxel grid 范围
    % --------------------------
    all_pts = [pts_gt; pts_seg];
    minXYZ = min(all_pts, [], 1);
    maxXYZ = max(all_pts, [], 1);
    
    grid_size = ceil((maxXYZ - minXYZ)/voxel_size) + 1;

    % 3. 体素化
    vox_gt  = voxelize(pts_gt, minXYZ, voxel_size, grid_size);
    vox_seg = voxelize(pts_seg, minXYZ, voxel_size, grid_size);
    
    % --------------------------
    % 4. 边界提取（3D边界）
    % --------------------------
    % 使用 26 邻域提取边界
    Bgt_vox = bwperim(vox_gt, 26);
    Bsg_vox = bwperim(vox_seg, 26);
    
    % --------------------------
    % 5. 边界 voxel 转回空间坐标
    % --------------------------
    [idx_x, idx_y, idx_z] = ind2sub(grid_size, find(Bgt_vox));
    Bgt = [idx_x, idx_y, idx_z] - 0.5;          % voxel center
    Bgt = Bgt * voxel_size + minXYZ;
    
    [idx_x, idx_y, idx_z] = ind2sub(grid_size, find(Bsg_vox));
    Bsg = [idx_x, idx_y, idx_z] - 0.5;
    Bsg = Bsg * voxel_size + minXYZ;
    
    % --------------------------
    % 6. 最近邻距离计算
    % --------------------------
    dist_seg_to_gt = pdist2(Bsg, Bgt, 'euclidean', 'Smallest', 1);
    dist_gt_to_seg = pdist2(Bgt, Bsg, 'euclidean', 'Smallest', 1);
    
    % --------------------------
    % 7. 根据距离阈值计算 Precision/Recall
    % --------------------------
    tau = 1.5 * voxel_size;
    
    Precision = sum(dist_seg_to_gt <= tau) / numel(dist_seg_to_gt);
    Recall    = sum(dist_gt_to_seg <= tau) / numel(dist_gt_to_seg);
    
    % BF1
    if Precision + Recall > 0
        BoundaryF1 = 2 * Precision * Recall / (Precision + Recall);
    else
        BoundaryF1 = 0;
    end

end

% 体素化
function vox = voxelize(points, minXYZ, voxel_size, grid_size)
    idx = floor((points - minXYZ)/voxel_size) + 1;   % 体素索引
    % 防止越界
    idx = max(idx, 1);
    idx(:,1) = min(idx(:,1), grid_size(1));
    idx(:,2) = min(idx(:,2), grid_size(2));
    idx(:,3) = min(idx(:,3), grid_size(3));
    
    lin = sub2ind(grid_size, idx(:,1), idx(:,2), idx(:,3));
    vox = false(grid_size);
    vox(lin) = true;
end

function voxel_grid = voxelize2(points, res)
% 将点云量化到体素网格，并返回唯一的体素坐标（整数索引）
% 输入:
%   points : N×3 点云矩阵
%   res    : 体素分辨率（标量）
% 输出:
%   voxel_grid : M×3 唯一体素索引（已去重）

% 归一化到体素坐标系（避免浮点误差）
voxel_indices = floor(points / res);  % 或使用 round()，但 floor 更稳定

% 去除重复体素（每个体素只计一次）
voxel_grid = unique(voxel_indices, 'rows');
end

% 从点云 Intensity 通道提取标签
function labels = getLabels(ptCloud)
    if isempty(ptCloud.Intensity)
        error('点云 Intensity 字段为空！');
    end

    labels = ptCloud.Intensity;   % 颜色非零视为 1
end

% 单株点云分割效果评价
function sigMetrics = evaluteSinglePlants(GT_pts, Seg_pts)
    % 提取位置和标签
    GT_xyz = GT_pts.Location;
    GT_lab = GT_pts.Intensity;
    Pre_xyz = Seg_pts.Location;
    Pre_lab = Seg_pts.Intensity;

    % 检查点数是否一致
    if size(GT_xyz,1) ~= size(Pre_xyz,1)
        error('GT 和 Seg 点云点数不一致，无法使用排序匹配');
    end

    % 通过 sortrows 排序位置（假设位置浮点值精确相同）
    [~, sortIdx_gt] = sortrows(GT_xyz);
    sorted_GT_lab = GT_lab(sortIdx_gt);
    [~, sortIdx_seg] = sortrows(Pre_xyz);
    sorted_Pre_lab = Pre_lab(sortIdx_seg);

    % 获取唯一的标签
    unique_labels = unique(sorted_GT_lab);
    num_classes = length(unique_labels);

    % 新增：计算每个标签下的点数（使用 histcounts 统计）
    % GT_counts(i) = unique_labels(i) 在 GT 中的点数
    GT_counts = histcounts(sorted_GT_lab, [unique_labels; inf]);
    % Pred_counts(i) = unique_labels(i) 在 Pred 中的点数（注意：Pred 可能有不同标签，但我们基于 GT 的 unique_labels 统计）
    Pred_counts = zeros(num_classes, 1);
    for i = 1:num_classes
        Pred_counts(i) = sum(sorted_Pre_lab == unique_labels(i));
    end

    % 新增：输出每个标签的点数
    disp('每个标签下的点数（基于 GT 的 unique_labels）：');
    for i = 1:num_classes
        fprintf('Label %d: GT 点数 = %d, Pred 点数 = %d\n', unique_labels(i), GT_counts(i), Pred_counts(i));
    end
    
    % 构建逐点的混淆矩阵
    % 类别 0 视作背景，其余是有效类别
    valid_idx = sorted_GT_lab > 0 | sorted_Pre_lab > 0;
    GT_eval = double(sorted_GT_lab(valid_idx));
    Pred_eval = double(sorted_Pre_lab(valid_idx));
    
    % C(i,j) = 真值为 i 类，被预测为 j 类
    C = confusionmat(GT_eval, Pred_eval, 'Order', 1:num_classes);
    
%     TP = diag(C);
%     FP = sum(C,1)' - TP;
%     FN = sum(C,2) - TP;
%     
%     Precision = TP ./ (TP + FP + eps);
%     Recall    = TP ./ (TP + FN + eps);
%     F1        = 2 * (Precision .* Recall) ./ (Precision + Recall + eps);
%     IoU       = TP ./ (TP + FP + FN + eps);

    Precision   = diag(C) ./ (sum(C, 1)' + eps); % 精度
    Recall      = diag(C) ./ (sum(C, 2) + eps);     % 召回率
    F1          = 2 * (Precision .* Recall) ./ (Precision + Recall + eps);
    IoU         = sum(diag(C)) / (sum(C(:)) + eps);

    % 宏平均
    P_macro     = mean(Precision);
    R_macro     = mean(Recall);
    F1_macro    = mean(F1);
    mIoU        = mean(IoU);

    sigMetrics = struct('Precision', P_macro, 'Recall', R_macro, 'F1', F1_macro, 'mIoU', mIoU);
end

function matchedFiles = findFilesByPattern(folderPath, fileExtension, regexPattern)
%FINDFILESBYPATTERN 查找符合指定扩展名和正则表达式的文件
%
% 输入:
%   folderPath    - 要搜索的文件夹路径（字符串）
%   fileExtension - 文件扩展名通配符，例如 '*.ply'
%   regexPattern  - 用于匹配文件名（含扩展名）的正则表达式，例如 '^\d{1,2}-\d\.ply$'
%
% 输出:
%   matchedFiles  - 符合条件的完整文件路径组成的 cell 数组（若无匹配则为空 cell）

    % 输入验证
    if nargin < 3
        error('需要提供三个输入参数：folderPath, fileExtension, regexPattern');
    end

    if ~ischar(folderPath) && ~isstring(folderPath)
        error('folderPath 必须是字符串');
    end
    folderPath = char(folderPath); % 兼容 string 类型

    if ~ischar(fileExtension) && ~isstring(fileExtension)
        error('fileExtension 必须是字符串');
    end
    fileExtension = char(fileExtension);

    if ~ischar(regexPattern) && ~isstring(regexPattern)
        error('regexPattern 必须是字符串');
    end
    regexPattern = char(regexPattern);

    % 获取目录中所有匹配扩展名的文件
    fileList = dir(fullfile(folderPath, fileExtension));

    if isempty(fileList)
        matchedFiles = {};
        return;
    end

    % 筛选符合正则表达式模式的文件
    matchedFiles = {};
    for i = 1:length(fileList)
        if fileList(i).isdir
            continue; % 跳过子目录
        end
        fullName = fileList(i).name; % 包含扩展名
        if ~isempty(regexp(fullName, regexPattern, 'once'))
            matchedFiles{end+1} = fullfile(folderPath, fullName);
        end
    end
end

