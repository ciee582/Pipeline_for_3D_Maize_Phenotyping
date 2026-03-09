function [finalSortedData, numFullGroups, remainingNum] = seg_box(data_info, groupSize)
    % 数据存储在一个名为 data_info 的矩阵中，大小为 Nx3，这里 N 为 60
    % data: (X, Y, label)  XY -- 相应坐标值； label -- 数据的类别标签
    % groupSize = 4
    % 返回值 numFullGroups 代表整除分组数，remainingNum 则表示整除分组的余数/剩余的不足一个分组的目标数
    % 不对 label 进行修改，只对的XY坐标进行升序排序

    N = size(data_info, 1);
    
    % Step 1: 按 Y 升序排序
    sortedData = sortrows(data_info, 1);

    % Step 2: 每 4 个一组，若不能完全划分，则还需要处理剩余数据
    numFullGroups = floor(N / groupSize); % 完整组的数量
    remainingData = sortedData((numFullGroups * groupSize + 1):end, :); % 剩余数据
    
    % 3. 对每个完整组，组内的数据按 X 坐标进行排序
    % 初始化finalSortedData的大小
    finalSortedData = zeros(N, size(data_info, 2)); % 使用zeros来初始化，确保大小正确
    % 完整组组内排序并合并数据
    idx = 1; % 用于追踪finalSortedData中的当前索引
    for i = 1:numFullGroups
        groupIdx = (i-1)*groupSize + (1:groupSize);
        groupData = sortedData(groupIdx, :);
        sortedGroup = sortrows(groupData, 2, 'descend'); % 对组内数据按X排序
        finalSortedData(idx:idx+groupSize-1, :) = sortedGroup; % 将排序后的组复制到finalSortedData中
        idx = idx + groupSize; % 更新索引
    end
    
    % 4. 处理剩余数据
    remainingNum = N - numFullGroups * groupSize; % 剩余数据的数量
    if remainingNum > 0
        % 剩余数据已经是按Y排序的，现在需要按X排序
        remainingDataSorted = sortrows(remainingData, 2, 'descend');
        
        % 将排序后的剩余数据附加到finalSortedData的末尾
        % idx 已达到剩余数据数量的最小值，所以在加上剩余数据量时要减掉 1
        % 例如，总数据 43 个，分10个完整组，最后的 idx=41，剩余3个数据，索引值应是41 42 43，即 41+(1:3)-1
        finalSortedData(idx + (1:remainingNum) -1, :) = remainingDataSorted;
    end
end