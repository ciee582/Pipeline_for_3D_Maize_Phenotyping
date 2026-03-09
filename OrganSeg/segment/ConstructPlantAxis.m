function [ PA_Pts, PA_Spls] = ConstructPlantAxis(points, spls, Phi_S, sub_skeletons)
% 建立一个以植物主干为基准的局部坐标系
%
% 可以达到以下目的：
% 标准化植物姿态： 其主干会沿着某个特定的轴（例如Z轴或Y轴）对齐，达到"扶正"效果，消除了姿态差异对后续过程的影响
% 简化计算： 在一个与植物结构对齐的坐标系中，许多与植物形态相关的计算（如长度、直径、分支角度等）会变得更简单、更直观。
% 特征提取： 在标准化坐标系中，更容易提取出与植物种类、生长状态等相关的几何特征。


% 计算实际茎秆的主方向
skeleton = sub_skeletons{end};
p1 = spls(skeleton(1), :);
p2 = spls(skeleton(end), :);
dir1 = (p2 - p1) ./ norm(p2 - p1);


% 拟合茎秆的主方向和中心点
[center, direction] = fitline(points(Phi_S, :));

% 校正拟合方向 direction 和实际方向 dir1 的一致性
% a = sum(dir1 .* direction);     % 两个向量的点积。如果点积小于0，说明两个方向向量大致相反
% if(a < 0)
%     direction = -1 * direction;
% end

% 构建新的坐标轴（Axis）
Axis = find_AxisByPrincipalDir_mt(points(Phi_S,:), dir1, center, false);

% 将所有点转换到新坐标系
PA_Pts = Transfer_XYZ2AXIS_mt(points, Axis);
PA_Spls = Transfer_XYZ2AXIS_mt(spls, Axis);

end

