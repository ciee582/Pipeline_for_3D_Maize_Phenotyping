function [spls, A, corresp] = edge_collapse_update(~, spls, corresp, ~, spls_adj, collapse_order)
    %% recover the set of edges on triangles & count triangles
    A = spls_adj;
    A(A>0) = 1;
    
    degrees = ones(size(spls,1),1);
    for i=1:length(spls)
        ns = find(A(i,:)==1);
        degrees(i) = length(ns)-1;
    end
    
    tricount = 0;
    skeds = zeros(0,4);
    for i=1:length(spls)
        ns = find(A(i,:)==1);
        ns = ns( ns>i );
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
        
    %% --- EDGE COLLAPSE 
    while true
        if size(skeds,1) == 0, break, end
        
        %--- DECIMATION STEP + UPDATES
        % collapse the edge with minimum cost, remove the second vertex
        if collapse_order == 1 % cost is degree + distance
            mind = min( skeds(:,3) );
            tmpIdx = find(skeds(:,3)==mind);
            tmpSkeds = skeds(tmpIdx,4);
    
            [~, idx] = min( tmpSkeds );
            edge = skeds(tmpIdx(idx),1:2);
            skeds(tmpIdx(idx),:)=[];
        else % cost is distance
            [~, idx] = min( skeds(:,4) );
            edge = skeds(idx,1:2);
            skeds(idx,:)=[];
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
        
        %% 1) remove skeds connect edge(1) and neighbor of edge(2), called 12; 
        tmpIdx = skeds( skeds(:,1)==edge(2), 2);
        tmpIdx = [tmpIdx; skeds( skeds(:,2)==edge(2), 1)];
        
        [rows,cols] = find(skeds(:,1:2)==edge(1));
        toBeRemoved =  zeros(0,1);   
        for i = 1:length(rows)
            col = 1 + mod(cols(i),2);
            if ismember( skeds(rows(i), col), tmpIdx )%remove
                toBeRemoved(end+1) = rows(i);
            else
                skeds(rows(i), cols(i)) = edge(2);
            end
        end
        if ~isempty(toBeRemoved)
            skeds(toBeRemoved,: ) = [];
        end
        
        %% 2) remove skeds which contain edge(2) and nolonger a edge of a triangle and
        ns = find( A(edge(2),:)==1 );
        ns = ns( ns~=edge(2) );
        lns = length(ns);
        % triangles contain edges(2) include both the two kinds of tmpEdges 
        tmpEdges = zeros(0,2); % edges contain edge(2)
        tmpEdges1 = zeros(0,2); % edges not contain edge(2)
        for j=1:lns
            for k=j+1:lns
                if A(ns(j),ns(k)) == 1
                    tmpEdges(end+1,:) = [edge(2),ns(j)];  
                    tmpEdges(end+1,:) = [edge(2),ns(k)];
                    tmpEdges1(end+1,:) = [ns(j),ns(k)];                
                end
            end
        end
        
        % 2.2) remove all edges do not belong to tmpEdges
        [rows,cols] = find(skeds(:,1:2)==edge(2));
        toBeRemoved =  zeros(0,1);
        tobedel = zeros(0,1); 
        for j = 1:length(rows)
            col = 1 + mod(cols(j),2);
            tmp = find( tmpEdges(:,2) == skeds(rows(j),col) );
            if tmp % is an edge of tmpEdge, then remove it from tmpEdge    
                for k = 1:length(tmp)
                    if ~ismember (tmp(k), tobedel)
                        tobedel = [tobedel; tmp(k)];
                    end
                end
            else
                toBeRemoved(end+1) = rows(j);
            end
        end
        if ~isempty(toBeRemoved)
            skeds( toBeRemoved,: ) = [];
        end
        if ~isempty(tobedel)
            tmpEdges(tobedel,:) = [];
        end
        
        % 2.3) add triangle edges new formed
        tmpEdges = [tmpEdges; tmpEdges1];
        for j = 1:size(tmpEdges, 1)
            tedge = tmpEdges(j,:);
            [rows,cols] = find(skeds(:,1:2)==tedge(1));
            bin = false;
            for k = 1:length(rows)
                col = 1 + mod(cols(k),2);
                if skeds(rows(k),col)==tedge(2) % the edge already in skeds
                    bin = true;
                    break;
                end
            end
            if ~bin % add edge to skeds
                ns = find(A(tedge(1),:)==1);
                degrees(tedge(1)) = (length(ns)-1)*0.5;
                ns = find(A(tedge(2),:)==1);
                degrees(tedge(2)) = (length(ns)-1)*0.5;
                skeds(end+1,1:2) = tedge;
                skeds(end,3) = 0.5*(degrees(tedge(1))+degrees(tedge(2)));
                skeds(end,4) = euclidean_distance(spls(tedge(1),:), spls(tedge(2),:) );             
            end
        end
        
        %% 3) update distance and degree of edges contain edge(2)
        [rows,cols] = find(skeds(:,1:2)==edge(2)); 
        ns = find(A(edge(2),:)==1);
        degrees(edge(2)) = (length(ns)-1)*0.5;
        
        for j = 1:length(rows)
            col = 1 + mod(cols(j),2);   
            k = skeds(rows(j),col);
            ns = find(A(k,:)==1);
            degrees(k) = (length(ns)-1)*0.5;
            
            skeds(rows(j),3) = 0.5*(degrees(edge(2))+degrees(k));
            skeds(rows(j),4) = euclidean_distance(spls(edge(2),:), spls(k,:) ); 
        end
    end
end

%% euclidean_distance
function dist = euclidean_distance(p1, p2)
    v = p1 - p2;
    dist = sqrt(dot(v, v));
end