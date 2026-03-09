%% init_env
clear;
clc;

% 读取点云
% pathname = 'E:\Datasets\3d\20241012\';
% filename = 'maize-60-Cloud-SOR-plants222-optimized-Oct11.ply';

pathname = 'E:\Datasets\3d\20241026\';
filename = '目标42.las.segmented-Octree11-2.ply';

ptCloud = pcread([pathname filename]);

% 运行 SMRF
% [groundIdx, nonGroundIdx, groundPts, nonGroundPts] = smrfGroundFilter2(ptCloud, 'c', 0.8, 's', 0.1, 'w', 10, 'et', 0.5, 'es', 1.5);

%% 通用点云地面滤波 (基于 SMRF)
% 目的：读取 PLY 文件，分割地面点和非地面点

% --- 步骤 1: 清理环境并设置文件路径 ---
% 检查点云是否成功读取
if ptCloud.Count == 0
    error('点云文件读取成功，但点数为零。请检查文件内容。');
end

% --- 步骤 3: 使用 SMRF 算法进行地面分割 ---
% SMRF (Simple Morphological Filter) 是一种常用的地面滤波方法。
% 它的原理是基于形态学开运算和坡度阈值，与PTIN一样是表面拟合的优秀替代。

% **参数解释 (请根据您的玉米点云场景调整)：**
% 1. 'MaxWindowRadius' (最大窗口半径): 
%    - 决定了形态学操作的结构元素大小。应略大于您玉米植株的最大直径。
%    - 对于田间点云，通常需要比默认值（18m）小得多，例如 0.5m ~ 2m。
% 2. 'ElevationThreshold' (高程阈值): 
%    - 用于区分地面和非地面对象的高度阈值。
%    - 决定了地面TIN与植株底部点之间的最大允许高度差。
% 3. 'SlopeThreshold' (坡度阈值): 
%    - 用于判断区域是否为地面的最大允许坡度。

disp('正在执行 SMRF 地面分割...');

% 推荐的参数调整（适用于低矮作物和相对平坦的田间）：
% 请根据您的数据特点调整这些值！
groundPtsIdx = segmentGroundSMRF(ptCloud, ...
    'MaxWindowRadius', 2.0, ...       % 假设植株直径不超过 1.0 米
    'ElevationThreshold', 0.12, ...    % 假设地面点和非地面点的Z坐标差异大于 0.2 米
    'SlopeThreshold', 0.6);           % 允许一定的地面坡度 (0.5 约等于 26.5 度)
% --- 步骤 4: 提取地面点和非地面点 ---
    
% 提取地面点云对象
groundPtCloud = select(ptCloud, groundPtsIdx);

% 提取非地面点云对象 (即植株/作物点)
nonGroundPtCloud = select(ptCloud, ~groundPtsIdx);

disp(['地面点数: ', num2str(groundPtCloud.Count)]);
disp(['非地面点数 (植株): ', num2str(nonGroundPtCloud.Count)]);

%% 保存文件
% 获取位置数据
groundLocations = groundPtCloud.Location;
nonGroundLocations = nonGroundPtCloud.Location;
% 创建 Intensity 数据
groundIntensity = zeros(groundPtCloud.Count, 1, 'uint8');  % Intensity=0 for ground
nonGroundIntensity = ones(nonGroundPtCloud.Count, 1, 'uint8');  % Intensity=1 for non-ground
% 合并位置和 Intensity
mergedLocations = [groundLocations; nonGroundLocations];
mergedIntensity = [groundIntensity; nonGroundIntensity];
% 创建新的 pointCloud 对象，包含 Intensity
mergedPtCloud = pointCloud(mergedLocations, 'Intensity', mergedIntensity);
% 保存为二进制 PLY 文件
% pcwrite(mergedPtCloud, [pathname 'plants242-SMRF-withLables.ply'], 'Encoding', 'binary');

[~, name, ~] = fileparts(filename);
pcwrite(mergedPtCloud, [pathname name '-SMRF2.ply'], 'PLYFormat', 'binary');
disp(['结果已保存至: ', [pathname name '-SMRF2.ply']]);

%% --- 步骤 5: 结果可视化 ---    
pcshow(groundPtCloud.Location, 'MarkerSize', 10, 'VerticalAxis', 'Z', 'VerticalAxisDir', 'up');
title('地面点云 (Ground Points)');
% 可视化颜色可以帮助区分
colorMapGround = [0 0.7 0]; % 绿色
if groundPtCloud.Count > 0
    pcshow(groundPtCloud.Location, repmat(colorMapGround, groundPtCloud.Count, 1));
end

pcshow(nonGroundPtCloud.Location, 'MarkerSize', 10, 'VerticalAxis', 'Z', 'VerticalAxisDir', 'up');
title('非地面点云 / 植株 (Non-Ground/Plant Points)');
% 可视化颜色可以帮助区分
colorMapNonGround = [1 0 0]; % 红色
if nonGroundPtCloud.Count > 0
    pcshow(nonGroundPtCloud.Location, repmat(colorMapNonGround, nonGroundPtCloud.Count, 1));
end

%% 打印点数
fprintf('地面点: %d, 非地面点: %d\n', numel(groundIdx), numel(nonGroundIdx));

%% 保存文件
pcwrite(groundPts, [pathname 'plants243-ground.ply']);

pcwrite(nonGroundPts, [pathname 'plants243-noground.ply']);

%% 主函数
function [groundIdx, nonGroundIdx, groundPts, nonGroundPts] = smrfGroundFilter2(ptCloud, varargin)
%SMRFGROUNDFILTER Simple Morphological Filter (SMRF) ground segmentation
%
%   [groundIdx, nonGroundIdx, groundPts, nonGroundPts] = smrfGroundFilter(ptCloud)
%   [...] = smrfGroundFilter(ptCloud,'c',cellSize,'s',slopeThreshold,'w',wkmax,...)
%
%   INPUT
%       ptCloud          - pointCloud object
%       Name-Value pairs (optional):
%           'c'  cellSize      - grid resolution (m)          [required]
%           's'  slopeThreshold- max terrain slope (m/m)      [required]
%           'w'  wkmax         - max window size (m)          [required]
%           'et' elevationThreshold - max height above ground (m) [optional]
%           'es' elevationScaler - slope-dependent height scaler [optional]
%           'inpaintmethod'    - inpainting method (default 4)
%           'cutnet'           - cut net size (m) [optional]
%
%   OUTPUT
%       groundIdx      - linear indices of ground points in ptCloud.Location
%       nonGroundIdx   - linear indices of non-ground points
%       groundPts      - pointCloud containing only ground points
%       nonGroundPts   - pointCloud containing only non-ground points
%
%   Reference:
%       Pingel, T.J., Clarke, K.C., McBride, W.A., 2013.
%       An improved simple morphological filter for the terrain
%       classification of airborne LIDAR data. ISPRS J. Photogramm.
%       Remote Sens. 77, 21–30.

%% 参数检查
if ~isa(ptCloud,'pointCloud')
    error('First input must be a pointCloud object.');
end

% 解析参数
p = inputParser;
addParameter(p,'c',[],@isscalar);  % cellSize
addParameter(p,'s',[],@isscalar);  % slopeThreshold
addParameter(p,'w',[],@isscalar);  % wkmax
addParameter(p,'et',[],@isscalar); % elevationThreshold
addParameter(p,'es',[],@isscalar); % elevationScaler
addParameter(p,'inpaintmethod',4,@isscalar);
addParameter(p,'cutnet',[],@isscalar);
parse(p,varargin{:});

cellSize = p.Results.c;
slopeThreshold = p.Results.s;
wkmax = p.Results.w;
elevationThreshold = p.Results.et;
elevationScaler = p.Results.es;
inpaintMethod = p.Results.inpaintmethod;
cutNetSize = p.Results.cutnet;

if isempty(cellSize)
    error('Cell size (''c'') must be specified.');
end
if isempty(slopeThreshold)
    error('Slope threshold (''s'') must be specified.');
end
if isempty(wkmax)
    error('Maximum window size (''w'') must be specified.');
end
if ~isempty(elevationThreshold) && isempty(elevationScaler)
    elevationScaler = 0;
end

% 检查依赖函数
if ~exist('inpaint_nans','file')
    error('inpaint_nans.m is required. Download from: https://www.mathworks.com/matlabcentral/fileexchange/4551');
end

%% 网格化
xyz = double(ptCloud.Location);  % 强制 double，避免类型问题
x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);

% 确定网格范围
minX = min(x); maxX = max(x);
minY = min(y); maxY = max(y);
originX = floor(minX/cellSize)*cellSize;
originY = floor(minY/cellSize)*cellSize;
nCols = ceil((maxX - originX)/cellSize);
nRows = ceil((maxY - originY)/cellSize);
xi = originX + (0:nCols-1)*cellSize;
yi = originY + (0:nRows-1)*cellSize;

% 点云 → 网格索引
colIdx = floor((x - originX)/cellSize) + 1;
rowIdx = floor((y - originY)/cellSize) + 1;
gridIdx = (rowIdx-1)*nCols + colIdx;

%% 创建 DSM（最低点）
[gridList,~,ic] = unique(gridIdx);
zMin = accumarray(ic, z, [], @min, NaN);
ZImin = nan(nRows, nCols);
ZImin(gridList) = zMin;

% 填补空网格
if any(isnan(ZImin(:)))
    ZImin = inpaint_nans(ZImin, inpaintMethod);
end

%% 检测低点异常
[isLowOutlierCell] = progressiveFilter(-ZImin, 'c', cellSize, 's', 5, 'w', 1);

%% 分割网格（若提供 cutNetSize）
if ~isempty(cutNetSize)
    [ZInet, isNetCell] = createNet(ZImin, cellSize, cutNetSize);
else
    ZInet = ZImin;
    isNetCell = logical(zeros(size(ZImin)));
end

%% 检测对象
[isObjectCell] = progressiveFilter(ZInet, 'c', cellSize, 's', slopeThreshold, 'w', wkmax);

%% 构建地面模型
ZIpro = ZImin;
ZIpro(isLowOutlierCell | isObjectCell | isNetCell) = NaN;
ZIpro = inpaint_nans(ZIpro, inpaintMethod);
isObjectCell = isLowOutlierCell | isObjectCell | isNetCell;

%% 分类点
groundMask = false(ptCloud.Count,1);
nonGroundMask = false(ptCloud.Count,1);

if ~isempty(elevationThreshold)
    % 计算地面坡度
    [gx, gy] = gradient(ZIpro / cellSize);
    gsurfs = sqrt(gx.^2 + gy.^2); % 坡度
    clear gx gy
    
    % 插值到每个点
    [r, c] = map2pix([originX, maxX, originY, maxY, cellSize], x, y);
    ez = interp2(xi, yi, ZIpro, c, r, 'spline');
    SI = interp2(xi, yi, gsurfs, c, r, 'spline');
    
    % 根据高度阈值分类
    requiredValue = elevationThreshold + (elevationScaler * SI);
    groundMask = abs(ez - z) <= requiredValue;
    nonGroundMask = ~groundMask;
else
    % 若无高度阈值，基于格网分类
    for k = 1:numel(gridList)
        gIdx = gridList(k);
        [r,c] = ind2sub([nRows,nCols], gIdx);
        ptsInCell = (gridIdx == gIdx);
        if ~isObjectCell(r,c)
            % 接近地面高度的点（默认阈值 0.3 m）
            zCell = z(ptsInCell);
            surfaceZ = ZIpro(r,c);
            groundMask(ptsInCell & (abs(zCell - surfaceZ) <= 0.3)) = true;
        end
        nonGroundMask(ptsInCell & ~groundMask(ptsInCell)) = true;
    end
end

%% 输出
groundIdx = find(groundMask);
nonGroundIdx = find(nonGroundMask);
groundPts = select(ptCloud, groundIdx);
nonGroundPts = select(ptCloud, setdiff(1:ptCloud.Count, groundIdx));
end

%% 子函数：createDSM（简化版，内置）
function [ZI, R, isEmptyCell, xi, yi] = createDSM(x, y, z, varargin)
    p = inputParser;
    addParameter(p,'c',[],@isscalar);
    addParameter(p,'xi',[],@isvector);
    addParameter(p,'yi',[],@isvector);
    addParameter(p,'type','min');
    addParameter(p,'inpaintMethod',4,@isscalar);
    parse(p,varargin{:});
    
    cellSize = p.Results.c;
    xi = p.Results.xi;
    yi = p.Results.yi;
    type = p.Results.type;
    inpaintMethod = p.Results.inpaintMethod;
    
    if isempty(xi) || isempty(yi)
        minX = min(x); maxX = max(x);
        minY = min(y); maxY = max(y);
        originX = floor(minX/cellSize)*cellSize;
        originY = floor(minY/cellSize)*cellSize;
        nCols = ceil((maxX - originX)/cellSize);
        nRows = ceil((maxY - originY)/cellSize);
        xi = originX + (0:nCols-1)*cellSize;
        yi = originY + (0:nRows-1)*cellSize;
    else
        cellSize = abs(xi(2) - xi(1));
        nCols = length(xi);
        nRows = length(yi);
        originX = min(xi);
        originY = min(yi);
    end
    
    R = [originX, max(xi), originY, max(yi), cellSize];
    
    colIdx = floor((x - originX)/cellSize) + 1;
    rowIdx = floor((y - originY)/cellSize) + 1;
    gridIdx = (rowIdx-1)*nCols + colIdx;
    
    [gridList,~,ic] = unique(gridIdx);
    if strcmp(type,'min')
        zVal = accumarray(ic, z, [], @min, NaN);
    else
        zVal = accumarray(ic, z, [], @mean, NaN);
    end
    
    ZI = nan(nRows, nCols);
    ZI(gridList) = zVal;
    isEmptyCell = isnan(ZI);
    if any(isEmptyCell(:))
        ZI = inpaint_nans(ZI, inpaintMethod);
    end
end

%% 子函数：createNet（简化版，内置）
function [ZInet, isNetCell] = createNet(ZI, cellSize, cutNetSize)
    [nRows, nCols] = size(ZI);
    netCellsX = ceil(cutNetSize/cellSize);
    netCellsY = netCellsX;
    
    isNetCell = false(nRows, nCols);
    for i = 1:netCellsY:nRows
        for j = 1:netCellsX:nCols
            isNetCell(i:min(i+netCellsY-1,nRows), j:min(j+netCellsX-1,nCols)) = true;
        end
    end
    ZInet = ZI;
    ZInet(isNetCell) = NaN;
    ZInet = inpaint_nans(ZInet, 4);  % 默认 Springs
end

%% 子函数：map2pix（内置）
function [r, c] = map2pix(R, x, y)
    r = (R(4) - y) / R(5) + 1;
    c = (x - R(1)) / R(5) + 1;
end


function [isObjectCell,lastSurface,thisSurface] = progressiveFilter(lastSurface,varargin)
% Define required input parameters
cellSize = [];
slopeThreshold = [];
wkmax = [];

% Define optional input parameters
strelShape = [];


%% Process supplied arguments

i = 1;
while i<=length(varargin)    
    if isstr(varargin{i})
        switchstr = lower(varargin{i});
        switch switchstr
            case 'c'
                cellSize = varargin{i+1};
                i = i + 2;
            case 's'
                slopeThreshold = varargin{i+1};
                i = i + 2;
            case 'w'
                wkmax = varargin{i+1};
                i = i + 2;
            case 'shape'
                strelShape = varargin{i+1};
                i = i + 2;
            otherwise
                i = i + 1;
        end
    else
        i = i + 1;
    end
end


%% Catch some errors

if isempty(cellSize)
    error('Cell size must be specified.');
end
if isempty(wkmax)
    error('Maximum window size must be specified');
end
if isempty(slopeThreshold)
    error('Slope threshold value must be specified.');
end


%% Define some default parameters
if isempty(strelShape)
    strelShape = 'disk';
end


%% Convert wkmax to a vector of window sizes (radii) defined in pixels.  
% If w was supplied as a vector, use those values as the basis; otherwise,
% use 1:1:wkmax

if numel(wkmax)~=1
    wk = ceil(wkmax / cellSize);
else
    wk = 1 : ceil(wkmax / cellSize);
end

% wk = wkmax;
%% Define elevation thresholds based on supplied slope tolerance

eThresh = slopeThreshold * (wk * cellSize);


%% Perform iterative filtering.

isObjectCell = false(size(lastSurface));

for i = 1:length(wk)
    thisSurface = imopen(lastSurface,strel(strelShape,wk(i)));
    isObjectCell = isObjectCell | (lastSurface - thisSurface > eThresh(i));
    lastSurface = thisSurface;
end
end