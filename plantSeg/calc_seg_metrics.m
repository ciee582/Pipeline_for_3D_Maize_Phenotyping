% calc_seg_metrics.m
%% 计算多个点云分割结果的平均 Precision / Recall / F1 / IoU
% 文件命名规则：
%   GT  :  如 0-1_GT.txt,  1-2_GT.txt,  10-5_GT.txt  （符合 ^\d{1,2}-\d+_GT\.txt$）
%   SEG :  如 0-1_XYZC.txt, 1-2_XYZC.txt ... （对应替换）

clear; clc;

% 1. 设置数据文件夹（修改为你的实际路径）
pathname = 'E:\Datasets\3d\20241012\';
dataDir = [pathname 'plants240\result_knn7'];

% pathname = 'E:\Datasets\3d\20241026\';
% dataDir = [pathname 'plants222\result_knn8'];

% pathname = 'E:\Datasets\3d\20241214\';
% dataDir = [pathname 'plants-Octree104\result_knn10'];

% 2. 查找所有 GT 文件
gtFiles = dir(fullfile(dataDir, '*-*_GT.txt'));   % 匹配 *-*_GT.txt
if isempty(gtFiles)
    error('未在指定目录找到任何 *_GT.txt 文件！');
end

% ================ 主循环：读取并同时保存每对文件的标签 ================
allTrue = [];      % 用于最后的全局检查（可选保留）
allPred = [];

% 新增：保存每一对文件的标签（cell数组）
labelsPerFile = cell(length(gtFiles), 2);   % 第1列存GT，第2列存PRED
validFileCount = 0;

fprintf('开始处理 %d 个 GT 文件...\n', length(gtFiles));

for i = 1:length(gtFiles)
    gtName  = gtFiles(i).name;
%     segName = strrep(gtName, '_GT.txt', '_XYZC.txt');
    segName = strrep(gtName, '_GT.txt', '_GT_XYZC_nobg.txt');
    gtPath  = fullfile(dataDir, gtName);
    segPath = fullfile(dataDir, segName);
    
    if ~isfile(segPath)
        warning('缺失对应分割文件，跳过：%s', segName);
        continue;
    end
    
    % 读取标签
%     gtLabels  = readLabelColumn(gtPath);
%     segLabels = readLabelColumn(segPath);
    [gtLabels, segLabels] = readAndMatchLabels(gtPath, segPath, 1e-6);
    
    if length(gtLabels) ~= length(segLabels)
        warning('点数不匹配，跳过 %s (%d vs %d 点)', gtName, ...
                length(gtLabels), length(segLabels));
        continue;
    end
    
    % 过滤只保留类别 1 和 2
    validIdx = ismember(gtLabels, [1 2 3 4]) & ismember(segLabels, [1 2]);
    gtLabels  = gtLabels(validIdx);
    segLabels = segLabels(validIdx);
    
    if isempty(gtLabels)
        warning('无有效标签，跳过：%s', gtName);
        continue;
    end
    
    validFileCount = validFileCount + 1;
    
    % 保存到 cell 数组，后面用来逐文件统计
    labelsPerFile{validFileCount, 1} = gtLabels;   % GT
    labelsPerFile{validFileCount, 2} = segLabels;  % PRED
    
    % 同时拼接到全局（如果你还想顺便看全局结果）
    allTrue = [allTrue; gtLabels];
    allPred = [allPred; segLabels];
    
    fprintf('(%3d/%d) 已读取 %s (%d 点)\n', i, length(gtFiles), gtName, length(gtLabels));
end

% 裁剪多余的空行
labelsPerFile = labelsPerFile(1:validFileCount, :);

if validFileCount == 0
    error('没有有效的文件对可以用于评估！');
end

fprintf('\n成功加载 %d 个有效文件，开始逐文件计算指标并平均...\n', validFileCount);

% ==================== 逐文件计算指标并取平均 ====================
% 预分配（固定两个类别）
prec_per_plants = zeros(validFileCount, 2);
rec_per_plants  = zeros(validFileCount, 2);
f1_per_plants   = zeros(validFileCount, 2);
iou_per_plants  = zeros(validFileCount, 2);

for i = 1:validFileCount
    gt  = labelsPerFile{i, 1};
    pred = labelsPerFile{i, 2};
    
    for c = 1:2         % 强制计算类别 1 和 2
        TP = sum((gt == c) & (pred == c));
        FP = sum((gt ~= c) & (pred == c));
        FN = sum((gt == c) & (pred ~= c));
        
        prec_per_plants(i,c) = TP / (TP + FP + eps);
        rec_per_plants(i,c)  = TP / (TP + FN + eps);
        f1_per_plants(i,c)   = 2*TP / (2*TP + FP + FN + eps);
        iou_per_plants(i,c)  = TP / (TP + FP + FN + eps);
    end
end

% 每个类别在所有文件上求平均 → 再对类别求平均（即 mPrecision / mIoU）
mPrecision = mean(mean(prec_per_plants, 1));   % 先平均文件，再平均类别
mRecall    = mean(mean(rec_per_plants,  1));
mF1        = mean(mean(f1_per_plants,   1));
mIoU       = mean(mean(iou_per_plants,  1));

% 同时可以输出每个类别的平均值（非常有用）
plantsPrec = mean(prec_per_plants, 1);
plantsRec  = mean(rec_per_plants,  1);
plantsF1   = mean(f1_per_plants,   1);
plantsIoU  = mean(iou_per_plants,  1);

%% ==================== 输出结果 ====================
fprintf('\n============== 逐文件(植株)平均结果 (Per-file Averaging) ==============\n');
fprintf('有效文件数      : %d\n', validFileCount);
fprintf('Mean Precision  : %.4f   (Stem: %.4f,  Leaf: %.4f)\n', mPrecision, plantsPrec);
fprintf('Mean Recall     : %.4f   (Stem: %.4f,  Leaf: %.4f)\n', mRecall,    plantsRec);
fprintf('Mean F1-Score   : %.4f   (Stem: %.4f,  Leaf: %.4f)\n', mF1,        plantsF1);
fprintf('Mean IoU (mIoU) : %.4f   (Stem: %.4f,  Leaf: %.4f)\n', mIoU,       plantsIoU);
fprintf('================================================================\n');

%% Visualize
%% ================== 绘制逐文件指标散点图 ==================
% prec_per_class : [N × 2]   每行一个文件，第1列为类1 Precision，第2列为类2
% rec_per_class  : [N × 2]
% f1_per_class   : [N × 2]
% iou_per_class  : [N × 2]
% validFileCount = N;

N = validFileCount;

% 提取每个类别的指标（方便绘图）
class1_IoU = iou_per_plants(:,1);
class2_IoU = iou_per_plants(:,2);
class1_Prec = prec_per_plants(:,1);
class2_Prec = prec_per_plants(:,2);
class1_Rec  = rec_per_plants(:,1);
class2_Rec  = rec_per_plants(:,2);
class1_F1   = f1_per_plants(:,1);
class2_F1   = f1_per_plants(:,2);

% 计算均值（用于画参考线）
mIoU1 = mean(class1_IoU);
mIoU2 = mean(class2_IoU);
mP1 = mean(class1_Prec);
mP2 = mean(class2_Prec);

% ============================================
%  Precision / Recall / F1 三联箱线图
% =============================================

% 合并三个指标
metrics_stem = [class1_Prec, class1_Rec, class1_F1, class1_IoU];
metrics_leaf = [class2_Prec, class2_Rec, class2_F1, class2_IoU];

% ==== 绘图参数 ====
boxW = 0.5;            % 缩小箱体宽度（默认约 0.5）
jitterW = 0.15;         % 控制散点横向扩散，不依赖数据方差
dotSize = 28;           % 散点大小
alphaDot = 0.55;        % 散点透明度
alphaBox = 0.25;        % 箱线透明度

% ==== 配色（Precision / Recall / F1）====
% C = [0.27 0.54 0.99;    % 蓝色 Precision
%      1.00 0.55 0.00;    % 橙色 Recall
%      0.13 0.70 0.13];   % 绿色 F1
C = [0 0.4470 0.7410;    % 鲜明蓝色
     0.6350 0.0780 0.1840; % 默认红色
     0.27 0.54 0.99;
     0.4660 0.6740 0.1880]; % 鲜明绿色

figure('Color','w', 'Position',[295, 200, 920, 735]);

% ================= STEM =================
subplot(1,2,1)
hold on;
% h = boxplot(metrics_leaf, 'Labels', {'Prec(leaf)','Rec(leaf)','F1(leaf)'}, 'Widths', boxW);
h = boxplot(metrics_stem, 'Colors','k', 'Widths', boxW);
set(h,'LineWidth',1.2);

% --- Swarm scatter ---
% C = [0.27 0.54 0.99; 1.00 0.55 0.00; 0.13 0.70 0.13]; % Precision, Recall, F1
for i = 1:4
    x = i + (rand(size(metrics_stem(:,i)))-0.5)*2*jitterW;
    scatter(x, metrics_stem(:,i), 28, C(i,:), 'filled', 'MarkerFaceAlpha',0.55)
end

% 半透明箱体颜色（RainCloud视觉感受）
boxes = findobj(gca,'Tag','Box');
for i = 1:length(boxes)
    patch(get(boxes(i),'XData'), get(boxes(i),'YData'), ...
        C(end-i+1,:), ...        % 注意：boxplot 顺序是反的，所以反向取色
        'FaceAlpha', alphaBox, 'EdgeColor','none');
end

% --- 设置两行标签 ---
set(gca, 'XTick', 1:4)
set(gca, 'XTickLabel', {'Prec','Rec','F1', 'IoU'})       % 第一行指标
xlabel('Stem','FontWeight','bold')                % 第二行对象名
ylabel('Score');
title('Stem Performance Metrics');
set(gca,'FontSize',14,'LineWidth',1.2);
hold off;

% ================= LEAF =====================
ax2 = subplot(1,2,2);
hold on;

% ---- Boxplot ----
h = boxplot(metrics_leaf, 'Symbol','k+', 'Widths', boxW);
set(h,'LineWidth',1.2);

% 半透明箱体
boxes = findobj(gca,'Tag','Box');
for i = 1:length(boxes)
    patch(get(boxes(i),'XData'), get(boxes(i),'YData'), ...
        C(end-i+1,:), ...
        'FaceAlpha', alphaBox, 'EdgeColor','none');
end

% ---- Swarm scatter（扩散度统一） ----
for i = 1:4
    x = i + (rand(size(metrics_leaf(:,i)))-0.5)*2*jitterW;
    scatter(x, metrics_leaf(:,i), 28, C(i,:), 'filled', 'MarkerFaceAlpha',0.55)
end

set(ax2, 'XTick', 1:4)
set(ax2, 'XTickLabel', {'Prec','Rec','F1', 'IoU'})
xlabel('Leaf','FontWeight','bold')

ymin = min(metrics_leaf(:)) - 0.02;   % 留点下边距
% ymax = max(metrics_leaf(:)) + 0.02;   % 留点上边距
% ylim([0.86 1.0])
ylim([ymin 1.0])

title('Leaf Performance Metrics')
set(gca,'FontSize',14,'LineWidth',1.2);
hold off;

%% 图4：综合子图（最适合论文）
figure('Position',[100 100 1300 900],'Color','w');

% subplot(1,2,1)
% scatter(class1_F1, class1_IoU, 80, 'c^', 'filled', 'MarkerFaceAlpha',0.7); hold on
% scatter(class2_F1, class2_IoU, 80, 'Marker', '*', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
% title('F1-Score vs IoU'); xlabel('F1'); ylabel('IoU'); grid on; xlim([0 1]); ylim([0 1]);
% legend('Stem','Leaf','Location','northwest');
% 
% subplot(1,2,2)
bar([mean([class1_IoU class2_IoU]); mean([class1_Prec class2_Prec]); ...
     mean([class1_Rec  class2_Rec]);  mean([class1_F1  class2_F1]) ]');
% set(gca,'XTickLabel',{'IoU','Precision','Recall','F1'},'FontSize',12);
set(gca, 'XTick', 1:2)
set(gca, 'XTickLabel', {'Prec   Rec    F1    IoU','Prec   Rec    F1    IoU'})
xlabel('Stem                                      Leaf','FontWeight','bold')
ylabel('平均值'); 
% title('各类别平均指标'); 
grid on;

sgtitle(sprintf('逐文件(植株)分割指标分析（共 %d 个有效文件）', N), 'FontSize',14,'FontWeight','bold');
% print('-dpng', '-r300', 'PerFile_Metrics_Summary.png');


%% ================== 对全局标签 allTrue 和 allPred 的可视化 ==================
% 假设 allTrue 和 allPred 已经是列向量（N_total x 1），值只含 1 和 2
fprintf('\n开始绘制全局标签 (allTrue vs allPred) 可视化，总点数：%d\n', length(allTrue));

% 确保是列向量
trueLabels = allTrue(:);
predLabels = allPred(:);

%% 1. 混淆矩阵（Confusion Matrix）—— 最重要的一张图
figure('Position',[100 100 700 600],'Color','w');
cm = confusionmat(trueLabels, predLabels, 'Order',[1 2 3 4]);

cm_reduced = cm(:, 1:2);

% 绘制带数值和百分比的热力图
imagesc(cm_reduced); colormap(flipud(gray)); colorbar;
axis square; title('全局混淆矩阵（所有点）','FontSize',16);
xlabel('预测标签','FontSize',14); 
ylabel('真实标签','FontSize',14);
set(gca,'XTick',1:2,'XTickLabel',{'Stem','Leaf'},...
        'YTick',1:4,'YTickLabel',{'Stem','Leaf', 'Tassel', 'Ear'},'FontSize',12);

% 在格子中标注数量和百分比
[rows, cols] = size(cm_reduced);
for i = 1:rows
    for j = 1:cols
        text(j, i, sprintf('%d\n(%.1f%%)', cm_reduced(i,j), 100*cm_reduced(i,j)/sum(cm_reduced(:))), ...
            'HorizontalAlignment','center', 'FontSize',14, 'FontWeight','bold', ...
            'Color',ifelse(cm_reduced(i,j)>mean(cm_reduced(:)), 'w', 'k'));
    end
end
% print('-dpng','-r300','Global_Confusion_Matrix.png');

%% 2. 各类别点数统计柱状图（真实 vs 预测）
figure('Position',[100 100 800 550],'Color','w');
trueCounts = [sum(trueLabels==1), sum(trueLabels==2)];
predCounts = [sum(predLabels==1), sum(predLabels==2)];

barData = [trueCounts; predCounts]';
bar(barData, 'grouped');
colormap([0 0.45 0.74; 0.85 0.33 0.1]);
grid on; legend('真实标签','预测标签','Location','northwest');
set(gca,'XTickLabel',{'Stem','Leaf'},'FontSize',14);
ylabel('点数','FontSize',14);
title(sprintf('全局点数分布（总点数 %d）', length(trueLabels)),'FontSize',16);

% 在柱子上标数值
for i = 1:2
    text(i-0.15, trueCounts(i)+max(trueCounts)*0.02, num2str(trueCounts(i)), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color','k');
    text(i+0.15, predCounts(i)+max(predCounts)*0.02, num2str(predCounts(i)), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color','k');
end
% print('-dpng','-r300','Global_Class_Distribution.png');

%% 3. 标签一致性散点图（按点索引显示预测 vs 真实）
% 适合快速看出系统性错误（例如某段点云整体错分）
figure('Position',[100 100 1100 500],'Color','w');
hold on; grid on;
idx = 1:length(trueLabels);

% 正确预测的点（灰色底）
correct = trueLabels == predLabels;
plot(idx(correct), trueLabels(correct), '.', 'Color',[0.7 0.7 0.7], 'MarkerSize',3);

% 错误预测的点（红色突出）
wrong = trueLabels ~= predLabels;
plot(idx(wrong), predLabels(wrong), 'r.', 'MarkerSize',6);   % 画预测值
plot(idx(wrong), trueLabels(wrong), 'b.', 'MarkerSize',6);   % 画真实值（可选）

xlabel('点索引（全局顺序）','FontSize',14);
ylabel('标签','FontSize',14);
title('全局标签一致性（灰色=正确，红色=预测错误）','FontSize',16);
set(gca,'YLim',[0.5 2.5],'YTick',[1 2],'YTickLabel',{'Stem','Leaf'});
legend('正确预测','错误预测（显示预测值）','Location','best');
% print('-dpng','-r300','Global_Label_Consistency.png');

%% 4. 错误率热力图（可选：按文件顺序显示每个文件的错误率）
% 如果你还保留了 labelsPerFile，可以画这个（非常有用！）
if exist('labelsPerFile','var') && ~isempty(labelsPerFile)
    fileErrorRate = zeros(validFileCount,1);
    for i = 1:validFileCount
        gt = labelsPerFile{i,1};
        pr = labelsPerFile{i,2};
        fileErrorRate(i) = mean(gt ~= pr);
    end
    
    figure('Position',[100 100 900 400],'Color','w');
    bar(fileErrorRate, 'FaceColor',[0.85 0.33 0.1]);
    grid on;
    xlabel('文件序号','FontSize',14);
    ylabel('错误率','FontSize',14);
    title(sprintf('每个文件的分类错误率（平均 %.2f%%）', 100*mean(fileErrorRate)),'FontSize',16);
    ylim([0 max(0.05, max(fileErrorRate)*1.2)]);
    % 在高错误率的文件上标注文件名（可选）
    thresh = mean(fileErrorRate);  % 错误率大于平均值就标注
    highErrIdx = find(fileErrorRate > thresh);
    for k = highErrIdx'
        text(k, fileErrorRate(k) + 0.01, gtFiles(k).name(1:end-7), ...
            'Rotation',0, 'FontSize',9, 'HorizontalAlignment','center');
    end
%     print('-dpng','-r300','PerFile_ErrorRate.png');
end

fprintf('全局标签可视化完成！已生成4张图（混淆矩阵、分布、一致性、逐文件错误率）。\n');


%% 
% rootName = 'E:\Datasets\3d\20241012\plants240\result_knn7\';
% rootName = 'E:\Datasets\3d\20241026\plants222\result_knn12\';
rootName = 'E:\Datasets\3d\20241214\plants-Octree104\result_knn10\';
pName = '06-3';
gtPath = [rootName, pName, '_GT.txt'];
segPath = [rootName, pName, '_GT_XYZC_nobg.txt'];
% segPath = [rootName, pName, '_XYZC.txt'];

gtData = load(segPath);
pts = gtData(:,1:3);
[gt, pred] = readAndMatchLabels(gtPath, segPath);

validIdx = ismember(gt, [1 2 3 4]) & ismember(pred, [1 2]);
gt  = gt(validIdx);
pred = pred(validIdx);

visualizeStemTPFPFN(pts, gt, pred, 1);


%% 辅助函数：条件颜色（MATLAB 没有内置 ifelse，用这个）
function c = ifelse(condition, trueVal, falseVal)
    if condition
        c = trueVal;
    else
        c = falseVal;
    end
end


%% ==============================
% 基于 XYZ 坐标精确匹配（最准确，推荐！）
function [gtLabels, segLabels] = readAndMatchLabels(gtPath, segPath, tol)
    % tol: 坐标匹配容差，默认 1e-6
    if nargin < 3, tol = 1e-6; end
    
    % 读取完整数据 [X Y Z C]
    gtData  = load(gtPath);   % 假设空格分隔
    segData = load(segPath);
    
    if size(gtData,2) < 4 || size(segData,2) < 4
        error('文件列数不足4列!\n');
    end
    
    gtXYZ  = gtData(:,1:3);
    segXYZ = segData(:,1:3);
    gtC    = round(gtData(:,4));
    segC   = round(segData(:,4));
    
    % 关键：用坐标匹配对应关系（支持轻微浮点误差）
    [~, ia, ib] = intersect(round(gtXYZ/tol), round(segXYZ/tol), 'rows', 'stable');
    
    if length(ia) < 0.9 * min(size(gtData,1), size(segData,1))
        warning('匹配点太少！可能点云严重不一致：%s\n', gtPath);
    end
    
    % 按匹配顺序重新排列标签
    gtLabels  = gtC(ia);
    segLabels = segC(ib);
    
    % 只保留有效类别 1 和 2
%     valid = ismember(gtLabels, [1 2]) & ismember(segLabels, [1 2]);
%     gtLabels  = gtLabels(valid);
%     segLabels = segLabels(valid);
    
    fprintf('匹配成功：%d / %d 点对齐（tol=%.0e）\n', length(gtLabels), size(gtData,1));
end


function visualizeStemTPFPFN(pts, gt, pred, stemClass)

    % 逻辑索引
    isTP = (gt == stemClass) & (pred == stemClass);
    isFP = (gt ~= stemClass) & (pred == stemClass);
    isFN = (gt == stemClass) & (pred ~= stemClass);

    figure; hold on; axis equal; axis off;
    
    % 背景（所有点，淡灰）
    scatter3(pts(:,1), pts(:,2), pts(:,3), ...
             3, [0.8 0.8 0.8], 'filled');

    % TP（绿色）
    scatter3(pts(isTP,1), pts(isTP,2), pts(isTP,3), ...
             12, 'g', 'filled');

    % FP（红色）
    scatter3(pts(isFP,1), pts(isFP,2), pts(isFP,3), ...
             12, 'r', 'filled');

    % FN（蓝色）
    scatter3(pts(isFN,1), pts(isFN,2), pts(isFN,3), ...
             18, 'b', 'filled');

    legend({'All points','TP (Stem)','FP (Stem)','FN (Stem)'}, ...
            'Units','normalized', ...     % 用归一化坐标，随窗口等比
            'Position',[0.70 0.18 0.20 0.16]); % [left bottom width height]

    title('Stem TP / FP / FN Spatial Distribution');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
end

