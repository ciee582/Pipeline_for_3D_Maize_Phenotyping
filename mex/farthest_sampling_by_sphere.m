function [spls, corresp] = farthest_sampling_by_sphere(pts, RADIUS)
    kdtree = kdtree_build( pts );
    spls = zeros( 0, 3 );
    corresp = zeros( length(pts), 1 );
    mindst = nan( length(pts), 1 ); % mindst(i) is the min distance of pts(i) to the sample piont corresp(i) 
    
    for k=1:length(pts)
        if corresp(k)~=0, continue, end
            
        %--- query all the points for distances
        mindst(k) = inf; % make sure picked first
        
        %--- initialize the priority queue
        while ~all(corresp~=0) %~isempty( find(corresp==0, 1) )
            [~, maxIdx] = max( mindst );
            if mindst(maxIdx) == 0
                break
            end
    
            % query its delta-neighborhood
            [nIdxs, nDsts] = kdtree_ball_query( kdtree,pts(maxIdx,:), RADIUS );%original
            % if maxIdx and all its neighborhood has been marked, skip ahead
            if all( corresp(nIdxs) ~= 0 )
                mindst(maxIdx) = 0; 
                continue;
            end
    
            % create new node and update (closest) distances
            spls(end+1,:) = pts(maxIdx,:); %#ok<AGROW>
            for i=1:length(nIdxs)
                if mindst(nIdxs(i))>nDsts(i) || isnan(mindst(nIdxs(i)))
                   mindst(nIdxs(i)) = nDsts(i);
                   corresp(nIdxs(i)) = size(spls,1);
                end
            end
        end
    end
    kdtree_delete( kdtree );
end
