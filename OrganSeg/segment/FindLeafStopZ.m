function [StopZ_upstream, StopZ] = FindLeafStopZ( sub_skeleton, joints, spls )
%{
输入：
    sub_skeleton: cell 数组，每个 cell 是一个分支（骨架）的点索引，第 end 个是主茎的骨架点序列。
    joints: 所有子骨架（通常是叶片）与主干连接的关节点的索引。
    spls: 所有骨架点的坐标（N×3 的矩阵）

输出：
    StopZ: 对应叶片连接点（joint）本身的 Z 坐标。
    StopZ_upstream: 对应每个叶片，在主茎骨架往上游走两步的点的 Z 值（如果可以找到）。
%}

% 主茎（stem）的骨架点
stem_skeleton = sub_skeleton{end};
StopZ_upstream = [];
StopZ = [];

% 提取主茎骨架中属于 joints 的索引点，这些点是主干骨架与各个叶片骨架的连接点。
joint_array = [];
for i = 1:length(stem_skeleton)
    v = stem_skeleton(i);
    if(ismember(v, joints))
        joint_array = [joint_array; v];
    end
end

% 对每个子骨架（叶片）寻找对应连接的 joint
for i = 1:length(sub_skeleton)-1
    skeleton_ = sub_skeleton{i};
    joint_ = skeleton_(end);            % 每个子骨架（叶片）的末尾是它与主茎连接的点
    id = find(joint_array == joint_);   % 当前子骨架末尾点是否与某个连接点重合
    if(isempty(id))         % 如果末尾点不在 joint_array 中，则计算其与所有主干连接点的距离，将最近的点作为匹配点
        jpts = spls(joint_array, :);
        skpts = spls(joint_, :);
        disv = jpts - skpts;
        dis_ = sqrt(disv(:,1) .* disv(:,1) + disv(:,2) .* disv(:,2) + disv(:,3) .* disv(:,3));
        [~, mid] = min(dis_);
        id = mid;
    end

    % 尝试往 joint_array 中往下找低 2 个位置的点的 Z 值作为 StopX。
    % 希望找到“比连接点再往下一点”的主茎点作为参考。
    id2 = id - 1;
    if(id2 > 0)
        StopZ_upstream = [StopZ_upstream; spls(joint_array(id2), 3)];
        %          V=spls(joint_array(id2),:)-spls(joint_array(id),:);
        %          StopV=[StopV;V./norm(V)];
    else
        StopZ_upstream = [StopZ_upstream; -Inf];
        %          id2=stem_skeleton(1);
        %          V=spls(id2,:)-spls(joint_array(id),:);
        %          StopV=[StopV;V./norm(V)];
    end

    id2 = id;

    v1 = spls(joint_array(id2), 3);
    % v2=spls(joint_array(id2+1),1);
    StopZ = [StopZ; v1];

end

end

