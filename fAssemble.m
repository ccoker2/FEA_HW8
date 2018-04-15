function [FGlobal] = fAssemble(nElems,FElement)
%Assemble the global force vector for static problems
m = nElems;
% n = 2*m+1;
n = m+1;
% FElement
reshape(FElement,length(FElement),1);

if length(FElement)==2
    
    % first element of force vector
    FG1 = sparse(1:n-1, 1, FElement(1)*ones(1,n-1)  ,n,1);
%     full(FG1)
    
    % second element of force vector
    FG2 =  sparse(2:n , 1, FElement(2)*ones(1,n-1),n,1);
%     full(FG2)

    FGlobal = FG1 + FG2;
%     full(FGlobal)
    
elseif length(FElement)==3
    ind = (1:1:n);
    ind2 = (3:2:n-2);
    ind = [ind,ind2];
    
    ind3 = repmat([1,2,3],1,m);
    ind = sort(ind);
   
    FGlobal = sparse(ind,1,FElement(ind3));
    
    % non-overlapping contributions
    
else
    error('Error: Only 1st and 2nd-order Lagrange shape functions have been implemented in this code')
end
end

