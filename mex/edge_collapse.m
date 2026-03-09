function [spls, A, corresp] = edge_collapse(~, spls, corresp, ~, spls_adj, collapse_order)    
    %% --- EDGE COLLAPSE 
    A = spls_adj;
    A(A>0) = 1;
    
    while true    
        %--- RECOVER CONNECTIVITY
        degrees = ones(size(spls,1),1);
        for i=1:size(spls,1)
            ns = find(A(i,:)==1);
            degrees(i) = length(ns)-1;
        end
        tricount = 0;
        skeds = zeros(0,4);% idx1, idx2, average degree of two end points, distance
        for i=1:length(spls)
            ns = find(A(i,:)==1);
            ns = ns( ns>i );%touch every triangle only once
            lns = length(ns);
            for j=1:lns
                for k=j+1:lns
                    if A(ns(j),ns(k)) == 1
                        tricount = tricount+1;
                        skeds(end+1,1:3) = [i,ns(j), 0.5*(lns+degrees(ns(j))) ];      
                        skeds(end,4) = euclidean_distance(spls(i,:), spls(ns(j),:) );
                        skeds(end+1,1:3) = [ns(j),ns(k), 0.5*(degrees(ns(j)) +degrees(ns(k)) )];
                        skeds(end,4) = euclidean_distance( spls(ns(j),:), spls(ns(k),:) );
                        skeds(end+1,1:3) = [i,ns(k), 0.5*(lns+degrees(ns(k))) ]; 
                        skeds(end,4) = euclidean_distance( spls(ns(k),:), spls(i,:) );
                    end
                end
            end
        end
        
        %--- STOP CONDITION
        % no more triangles? then the structure is 1D
        if tricount == 0, break, end;
        
        %--- DECIMATION STEP + UPDATES
        % collapse the edge with minimum cost, remove the second vertex
        if collapse_order == 1 % cost is degree + distance
            mind = min( skeds(:,3) );
            tmpIdx = find(skeds(:,3)==mind);
            tmpSkeds = skeds(tmpIdx,4);
    
            [~, idx] = min( tmpSkeds );
            edge = skeds(tmpIdx(idx),1:2);
        else % cost is distance
            [~, idx] = min( skeds(:,4) );
            edge = skeds(idx,1:2);
        end
        
        % update the location
        spls( edge(2),: ) = mean( spls( edge,: ) );
        spls( edge(1),: ) = NaN;
        % update the A matrix
        for k=1:size(A,1)
            if A(edge(1),k) == 1 
                A(edge(2),k)=1; 
                A(k,edge(2))=1; 
            end
        end
        % remove the row
        A(edge(1),:) = 0;
        A(:,edge(1)) = 0;
        % update the correspondents
        corresp(corresp==edge(1) ) = edge(2);
    end
end

function dist = euclidean_distance(p1, p2)
    v=p1-p2;
    dist = sqrt(dot(v,v));
end
