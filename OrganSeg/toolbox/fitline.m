function [xyz0, direction ] = fitline( lineData)
%FITLINE 此处显示有关此函数的摘要
%   此处显示详细说明
xyz0 = mean(lineData,1);
% xyz0=(lineData(1,:)),

% 协方差矩阵奇异变换，与拟合平面不同的是
% 所得直线的方向实际上与最大奇异值对应的奇异向量相同

centeredLine = bsxfun(@minus, lineData, xyz0);  % lineData 的每一行减去 xyz0
[U, S, V] = svd(centeredLine);      % 奇异值分解（SVD）

% 在几何意义上，V 的列向量表示点集在 3D 空间中的主方向，奇异值 S 的大小反映了点集在这些方向上的分布程度。
direction = V(:,1)';

end

