function idx = fFindClosestPosition(X,V)
% see friend function 'fFindClosest'

[~,idx] = min(abs(X-V));
Y = X(idx);
end
