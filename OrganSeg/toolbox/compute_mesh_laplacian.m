function L = compute_mesh_laplacian(vertex,face,type,options)
options.null = 0;
normalize = getoptions(options, 'normalize', 0);
symmetrize = getoptions(options, 'symmetrize', 1);

W = compute_mesh_weight(vertex,face,type,options);
n = size(W,1);
if symmetrize==1 && normalize==0
    L = diag(sum(W,2)) - W;
elseif symmetrize==1 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1/2)) * W * diag(sum(W,2).^(-1/2));
elseif symmetrize==0 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1)) * W;
else
    error('Does not work with symmetrize=0 and normalize=0');    
end

for i = 1:n
    tmp = abs(sum(L(i,i)));
    if tmp>10000
        L(i,:) = L(i,:)*10000/tmp;
    end
end
