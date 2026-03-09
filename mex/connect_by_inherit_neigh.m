function A = connect_by_inherit_neigh(pts, spls, corresp, neigh)
    if ~iscell(neigh)
        tmp = cell(size(neigh,1),1);
        for i = 1:size(neigh,1)
            tmp{i}=neigh(i,:);
        end
        neigh = tmp;  clear tmp;
    end
    
    A = zeros( length(spls), length(spls) );
    for pIdx=1:length(pts)    
        ns =neigh{pIdx};
        pc = corresp(pIdx);
        if pc == 0
            warning('some points have no correspondence');continue;
        end
        
        for nIdx=1:length(ns)
            nc = corresp(ns(nIdx));
            if nc == 0      
                warning('some points have no correspondence');continue;
            end
            if nc~=pc
                A(pc,nc) = A(pc,nc) + 1;
                A(nc,pc) = A(nc,pc) + 1;
                break;
            end
        end
    end
    
    %% if there are isolate points, connect it with its nearest neighbors.
    isopts = zeros(1,0);
    for i=1:size(A,1)
        A(i,i) = 1;
        if length( find( A(i,:)>0) ) == 1
            isopts(1, end+1) = i;
            warning('there are isolate points: %d', i);
        end
    end
    if isempty(isopts)
        spls_kdtree = kdtree_build( spls );
        for i=isopts
            neighs = kdtree_k_nearest_neighbors( spls_kdtree, spls(i,:), 2)';
            for j = neighs
                A(i,j) = 1;
                A(j,i) = 1;
            end
        end
        kdtree_delete( spls_kdtree );
    end
end
