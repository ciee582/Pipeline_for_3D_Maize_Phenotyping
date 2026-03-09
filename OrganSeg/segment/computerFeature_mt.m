function [features] = computerFeature_mt(AxisPts, K, SHOW_FEATURES)
% 输出: N×1 统一特征值
% 特征综合了线性度、平面度、球形度

    ptCloud = pointCloud(AxisPts);
    N = size(AxisPts,1);
    features = zeros(N,1);

    for i = 1:N
        pt = AxisPts(i,:);
        [indices, ~] = findNearestNeighbors(ptCloud, pt, K);
        n_pts = AxisPts(indices,:);

        % PCA
        [~, ~, latent] = pca(n_pts);

        % 排序 (确保 λ1 >= λ2 >= λ3)
        eigvals = sort(latent, 'descend');  
        if length(eigvals) < 3 || eigvals(1) == 0
            features(i) = 0;
            continue;
        end

        lambda1 = eigvals(1);
        lambda2 = eigvals(2);
        lambda3 = eigvals(3);

        % 局部几何特征
        linearity   = (lambda1 - lambda2) / lambda1;
        planarity   = (lambda2 - lambda3) / lambda1;
        sphericity  = lambda3 / lambda1;

        % 统一特征值 (线性度主导 + 惩罚项)
        features(i) = linearity - 0.5*planarity - 0.2*sphericity;
    end

if SHOW_FEATURES
    fmin = min(features);
    fmax = max(features);
    
    % 避免除零
    if fmax == fmin
        features_norm = zeros(size(features));
    else
        features_norm = (features - fmin) / (fmax - fmin);
    end
    
    % 转换为 colormap 索引 [1,256]
    idx = floor(features_norm * 255) + 1;  % 保证范围在 1~256

    figure('Name','TensorFeatures','NumberTitle','off');set(gcf,'color','white');movegui('southwest');

    colors= jet(256);
    scatter3(AxisPts(:,1), AxisPts(:, 2), AxisPts(:, 3), 10, colors(idx,:), 'filled');

    colorbar;
    axis off; axis equal; camorbit(0, 0, 'camera'); axis vis3d; view(-90, 0);
end
end

