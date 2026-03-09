function [earid, mainStemRegion, AutoRegions] = IdentifyingOrgans(AutoRegions, Pts, spls, organ_subskeleton)
    maxZ = max(Pts(:,3));
    minZ = min(Pts(:,3));
    Para = zeros(length(AutoRegions)-1, 4);  % 每个区域的 OBB尺寸 + 最高Z

    for i = 2:length(AutoRegions)
        ids = AutoRegions{i};
        organPts = Pts(ids,:);
        oz = max(organPts(:,3));
        if length(ids) > 30
            OBB = find_OBB_mt(organPts, 1, 1, 1, false);
            Para(i-1,:) = [OBB.ExtentLen' oz];
        else
            Para(i-1,:) = [Inf Inf Inf oz];
        end
    end

    % ---------------------
    % 主茎识别：通过方向向上 + 骨架长度来判断
    % ---------------------
    mainScores = zeros(length(organ_subskeleton)-1,1);
    for i = 2:length(organ_subskeleton)
        indices = organ_subskeleton{i};
        if isempty(indices), continue; end
        vi_end = indices(end); vi_mid = indices(1);
        v_end = spls(vi_end,:); v_mid = spls(vi_mid,:);
        dv = (v_end - v_mid) ./ norm(v_end - v_mid);
        cosZ = abs(dv(3));  % 越接近1越垂直
        mainScores(i-1) = cosZ * length(indices);
    end
    [~, mainStemIdx] = max(mainScores);
    mainStemRegion = mainStemIdx + 1;       % AutoRegions 从2开始

    % ---------------------
    % 叶片识别（多片）
    % ---------------------
    skeLen = zeros(length(organ_subskeleton)-1,1);
    cosV = zeros(length(organ_subskeleton)-1,1);
    scores = zeros(length(organ_subskeleton)-1,1);

    for i = 2:length(organ_subskeleton)
        indices = organ_subskeleton{i};
        if isempty(indices), continue; end
        vi_end = indices(2);
        vi_mid = indices(1);
        v_end = spls(vi_end,:); v_mid = spls(vi_mid,:);
        dv = (v_mid - v_end) ./ norm(v_mid - v_end);
        cosV(i-1) = dv(3);
        skeLen(i-1) = length(indices);
    end

    for i = 1:length(scores)
        s = 0;
        if Para(i,4) < (minZ + 2*(maxZ - minZ)/3), s = s + 1; end       % 离地低
        if cosV(i) > 0,                      s = s + 1; end             % 方向朝上
        if skeLen(i) < mean(skeLen),         s = s + 1; end             % 骨架不是最长（避免茎秆）
        if Para(i,3) < 0.1 * (maxZ - minZ),  s = s + 1; end             % 厚度小
        scores(i) = s;
    end

    candidate_ids = find(scores >= 1);  % 宽松策略：得分不低于1即可进入下一步
    candidate_ids = candidate_ids + 1;  % 对应 AutoRegions 的索引（从2开始）

    earid = candidate_ids;

    % ---------------------
    % 如果超过1片叶子，使用结构性PCA特征进一步排序
    % ---------------------
    if length(candidate_ids) > 1
        earNum = zeros(length(candidate_ids),1);
        features = zeros(length(candidate_ids),1);
        for i = 1:length(candidate_ids)
            I1 = AutoRegions{candidate_ids(i)};
            earPts = Pts(I1,:);
            ptCloud = pointCloud(earPts);
            earNum(i) = length(I1);
            if length(I1) < 30, continue; end
            feature = 0;
            for j = 1:size(earPts,1)
                pt = earPts(j,:);
                [indices,~] = findNearestNeighbors(ptCloud, pt, 32);
                n_pts = earPts(indices,:);
                [~,~,coffes] = pca(n_pts,'Algorithm','eig');
%                 e1 = coffes(1)/sum(coffes); 
                e2 = coffes(2)/sum(coffes); 
                e3 = coffes(3)/sum(coffes);
                feature = feature + (e3)/(e2);
            end
            features(i) = feature / length(I1);
        end

        % 若特征低于阈值，则认为不够明显，不作为耳叶
        refined_earid = candidate_ids(features >= 0.18);
        if isempty(refined_earid)
            [~, idx] = sort(earNum, 'descend');
            refined_earid = candidate_ids(idx);
        end
        earid = refined_earid;
    end

    % ---------------------
    % 对 earid 按照起点 Z 坐标排序（从下往上）
    % ---------------------
    leaf_base_z = zeros(length(earid), 1);
    for i = 1:length(earid)
        indices = organ_subskeleton{earid(i)};
        leaf_base_z(i) = mean(spls(indices,3));  % 起点的Z坐标
    end
    [~, sortIdx] = sort(leaf_base_z, 'ascend');  % 升序
    earid = earid(sortIdx);
    
    % ---------------------
    % 将主茎索引放置为第一位
    % ---------------------
    earid = [mainStemRegion; earid];
end