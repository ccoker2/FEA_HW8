function [KGlobal] = kAssemble(nElems,KElement)



if size(KElement,2)==3
%Assemble the global stiffness matrix
m = nElems;
n = 2*m+1;

% center diagonal
ind = (1:1:n);
    ind2 = (3:2:n-2);
    ind = [ind,ind2];
    
    clear ind2
    
    ind3 = repmat([1,2,3],1,m);
    ind = sort(ind);
    dK = diag(KElement);
    D = sparse(ind,1,dK(ind3));
    A = spdiags(D,0,n,n);
    clear ind ind2 ind3 D
    B = sparse(2:n  ,1:n-1,KElement(1,2)*ones(1,n-1),n,n);
    B = B+B';


% overlapping contributions from each element
%
% For 1st-order Lagrange
if size(KElement,1) == 2  
%     F = sparse(2:n-1,2:n-1,KElement(1,1)*ones(1,n-2),n,n);
    
%
% For 2nd-order Lagrange
elseif size(KElement,1) == 3              
%     ind = nonzeros((2:n-1).*(mod(2:n-1,3)==0));
%     F = sparse(ind,ind,KElement(1,1)*ones(1,length(ind)),n,n);

    ind = (1:1:n-2)';
    diag1 = KElement(3,1) * (ones(n-2,1).*mod(ind,2));
    C = spdiags(diag1,-2,n,n);
    C = C+C';


end

KGlobal = A+B+C;
elseif size(KElement)==2
    m = nElems;
    n = m+1;
%     KElement
    % top left of element matrix --> center diagonal.
    D = sparse(1:n-1, 1:n-1, KElement(1,1)*ones(1,n-1),n,n);
%     full(D)
    
    % lower diagonal (its transpose being the upper diagonal).
    E = sparse(2:n  , 1:n-1, KElement(2,1)*ones(1,n-1),n,n);
%     full(E)
    
    % upper diagonal
    G = sparse(2:n  ,1:n-1,KElement(1,2)*ones(1,n-1),n,n)';
%     full(G)
    
    % bottom left of element matrix --> center diagonal
    F = sparse(2:n,2:n,KElement(2,2)*ones(1,n-1),n,n);
%     full(F)
    
    % add the four components 
    KGlobal = E+D+F+G;
%     full(KGlobal)
else
    error('Only 1st and second-order Lagrange shape functions have been implemented.')
end

