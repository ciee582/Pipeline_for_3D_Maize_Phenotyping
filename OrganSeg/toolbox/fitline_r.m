function [res] = fitline_r( lineData )
%FITLINE 此处显示有关此函数的摘要
%   此处显示详细说明
xyz0 = mean(lineData, 1);       % 点集的中心点（均值）,拟合直线的参考点
% xyz0=(lineData(1,:));

centeredLine = bsxfun(@minus, lineData, xyz0);      % 中心化, 对 lineData 的每一行减去 xyz0

% 协方差矩阵奇异变换，与拟合平面不同的是
% 所得直线的方向实际上与最大奇异值对应的奇异向量相同
% 对中心化后的点集矩阵 centeredLine 进行奇异值分解（SVD）
% U 是 n * n 的正交矩阵，表示左奇异向量。
% S 是 n * 3 的对角矩阵，包含非负奇异值（按降序排列）。
% V 是 3 * 3 的正交矩阵，表示右奇异向量。
% 在几何意义上，V 的列向量表示点集在 3D 空间中的主方向，奇异值 S 的大小反映了点集在这些方向上的分布程度。
[U, S, V] = svd(centeredLine);

direction = V(:,1)';        % 第一列，对应于最大奇异值, 即拟合直线的主方向

res = [xyz0 direction];

end